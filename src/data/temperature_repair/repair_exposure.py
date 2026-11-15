#!/usr/bin/env python3
"""Deterministic, fail-closed exposure repair. This script fits no models."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

REPOSITORY = Path(__file__).resolve().parents[3]
REPAIR_DATA = REPOSITORY / "src" / "data" / "temperature_repair"
PANEL_DEFAULT = REPOSITORY / "artifacts" / "processed_inputs" / "panel_composite_humidity.parquet"
OUTPUT_DEFAULT = REPOSITORY / "artifacts" / "temperature_repair" / "panel_composite_humidity_exposure_repaired.parquet"
REPORT_DEFAULT = REPOSITORY / "artifacts" / "temperature_repair" / "input_validation_primary.json"
MANIFEST_DEFAULT = REPAIR_DATA / "assembly_manifest.json"
WEIGHTS_DEFAULT = REPAIR_DATA / "final_production_exposure_weights.csv"

EXPECTED = {
    "manifest": ("bae568a58996400c7b62d858445acbb0f597c391300e3b8f13bdfeb6c45e6cf5", 61826),
    "weights": ("2c5311c83064aa9ea4845f941163ad508d7a91305807247290549b26eeda51ae", 6935134),
    "primary_weather": ("19bc7332ee461afe927f410e552b42ff6956028d0e53175c47a633b12d3db294", 296424626),
    "mixed_contributions": ("c6571c6a619a820b6125eba11a747bf8bc0aa4eeaae29f38b46203334bba0622", 651759371),
    "old_panel": ("8791be152c870ab691754e930eeade1b2c82a48d591fc10db963612a61ef2031", None),
}
TOL = 1e-9


class GateError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check(ok: bool, name: str, detail: str, checks: list[dict]) -> None:
    checks.append({"gate": name, "passed": bool(ok), "detail": detail})
    if not ok:
        raise GateError(f"{name}: {detail}")


def finite_close(a, b, *, rtol=1e-10, atol=TOL) -> bool:
    aa = np.asarray(a, dtype="float64")
    bb = np.asarray(b, dtype="float64")
    return aa.shape == bb.shape and np.isfinite(aa).all() and np.isfinite(bb).all() and bool(np.allclose(aa, bb, rtol=rtol, atol=atol))


def verify_inputs(checks: list[dict], input_paths: dict[str, Path]) -> dict:
    hashes = {}
    for name, path in input_paths.items():
        expected_hash, expected_bytes = EXPECTED[name]
        check(path.is_file(), f"exists:{name}", str(path), checks)
        size = path.stat().st_size
        check(expected_bytes is None or size == expected_bytes, f"bytes:{name}", f"observed={size}; expected={expected_bytes}", checks)
        actual = sha256(path)
        hashes[name] = actual
        check(actual == expected_hash, f"sha256:{name}", f"observed={actual}; expected={expected_hash}", checks)
    return hashes


def resolve_path(path: Path) -> Path:
    return (REPOSITORY / path).resolve() if not path.is_absolute() else path.resolve()


def require_columns(table: pa.Table, expected: set[str], name: str, checks: list[dict]) -> None:
    missing = sorted(expected - set(table.column_names))
    check(not missing, f"schema:{name}", f"missing={missing}", checks)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--primary-weather", type=Path, required=True, help="Private, hash-verified primary region-day weather contribution parquet.")
    ap.add_argument("--mixed-contributions", type=Path, required=True, help="Private, hash-verified mixed region-day contribution parquet.")
    ap.add_argument("--old-panel", type=Path, default=PANEL_DEFAULT, help="Original panel, normally rebuilt from packaged parts.")
    ap.add_argument("--assembly-manifest", type=Path, default=MANIFEST_DEFAULT, help="Packaged fixed-weight assembly manifest.")
    ap.add_argument("--weights", type=Path, default=WEIGHTS_DEFAULT, help="Packaged fixed production grid weights.")
    ap.add_argument("--output", type=Path, default=OUTPUT_DEFAULT, help="New corrected panel path; existing outputs are never overwritten.")
    ap.add_argument("--report", type=Path, default=REPORT_DEFAULT, help="Validation report path.")
    ap.add_argument("--validate-only", action="store_true", help="Run all gates but do not write a corrected panel.")
    args = ap.parse_args()
    panel_path = resolve_path(args.old_panel)
    manifest_path = resolve_path(args.assembly_manifest)
    weights_path = resolve_path(args.weights)
    primary_path = resolve_path(args.primary_weather)
    mixed_path = resolve_path(args.mixed_contributions)
    output_path = resolve_path(args.output)
    report_path = resolve_path(args.report)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    input_paths = {
        "manifest": manifest_path,
        "weights": weights_path,
        "primary_weather": primary_path,
        "mixed_contributions": mixed_path,
        "old_panel": panel_path,
    }
    report: dict = {
        "schema": "exposure-input-validation-v1",
        "decision": "complete_fixed_weight_sum_land_plus_era5_single_level",
        "models_fitted": False,
        "input_paths": {},
        "output_path": None,
        "checks": [],
    }
    checks = report["checks"]
    try:
        hashes = verify_inputs(checks, input_paths)
        report["input_paths"] = {name: {"path": str(path), "sha256": hashes[name], "bytes": path.stat().st_size} for name, path in input_paths.items()}
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        check(manifest.get("cross_product_metric_assembly") == "add fixed-weight contributions; no denominator change or renormalization", "assembly_contract", str(manifest.get("cross_product_metric_assembly")), checks)
        check(manifest.get("missing_values_filled") is False and manifest.get("main_models_fitted") is False, "source_scope", "source manifest must report no missing-value filling and no fitted models", checks)

        metric_specs = {
            "tmean": ("temperature_mean_c", "land_temperature_mean_c_weighted_contribution", "era5_temperature_mean_c_weighted_contribution", "mixed_temperature_mean_c_weighted_contribution"),
            "tmax": ("temperature_max_c", "land_temperature_max_c_weighted_contribution", "era5_temperature_max_c_weighted_contribution", "mixed_temperature_max_c_weighted_contribution"),
            "rh": ("relative_humidity_mean_pct", "land_relative_humidity_mean_pct_weighted_contribution", "era5_relative_humidity_mean_pct_weighted_contribution", "mixed_relative_humidity_mean_pct_weighted_contribution"),
        }
        pcols = {"health_region_2023", "date"} | {s[0] for s in metric_specs.values()}
        mcols = {
            "health_region_2023", "local_date", "full_candidate_weight_total", "land_weight_total", "era5_weight_total",
            "land_instantaneous_eligible_weight", "era5_instantaneous_eligible_weight",
        } | {c for spec in metric_specs.values() for c in spec[1:]}
        ptab = pq.read_table(primary_path, columns=sorted(pcols))
        mtab = pq.read_table(mixed_path, columns=sorted(mcols))
        require_columns(ptab, pcols, "primary", checks)
        require_columns(mtab, mcols, "mixed", checks)
        p = ptab.to_pandas()
        m = mtab.to_pandas()
        p["date"] = pd.to_datetime(p["date"]).dt.normalize()
        m["date"] = pd.to_datetime(m.pop("local_date")).dt.normalize()
        for frame in (p, m):
            frame["health_region_2023"] = frame["health_region_2023"].astype(str)
        key = ["health_region_2023", "date"]
        check(len(p) == 3_451_500 and p[key].drop_duplicates().shape[0] == len(p), "primary_region_day_unique", f"rows={len(p)}; keys={p[key].drop_duplicates().shape[0]}", checks)
        check(len(m) == 3_451_500 and m[key].drop_duplicates().shape[0] == len(m), "mixed_region_day_unique", f"rows={len(m)}; keys={m[key].drop_duplicates().shape[0]}", checks)
        check(p["health_region_2023"].nunique() == 450 and p.date.min() == pd.Timestamp("2003-01-01") and p.date.max() == pd.Timestamp("2023-12-31"), "source_dates_regions", "expected 450 regions and daily 2003-01-01 through 2023-12-31", checks)
        check(set(map(tuple, p[key].to_numpy())) == set(map(tuple, m[key].to_numpy())), "source_key_complement", "primary and mixed keys must be identical", checks)
        x = p.merge(m, on=key, validate="one_to_one", suffixes=("_primary", "_mixed"))

        w = pd.read_csv(weights_path)
        wcols = {"health_region_2023", "latitude", "longitude", "source_product", "population_weight"}
        check(wcols.issubset(w.columns), "weight_schema", f"observed={list(w.columns)}", checks)
        w["health_region_2023"] = w["health_region_2023"].astype(str)
        check(not w.duplicated(["health_region_2023", "latitude", "longitude"]).any(), "unique_region_grid_coordinates", "duplicate region/grid coordinates absent; this does not independently establish population-support accounting", checks)
        check(np.isfinite(w.population_weight.to_numpy(dtype=float)).all() and (w.population_weight >= 0).all(), "finite_nonnegative_fixed_weights", "all supplied grid weights must be finite and nonnegative", checks)
        check(set(w["source_product"].dropna()) == {"reanalysis-era5-land", "reanalysis-era5-single-levels"}, "weight_sources", str(sorted(w["source_product"].dropna().unique())), checks)
        wsum = w.groupby("health_region_2023").population_weight.sum()
        check(len(wsum) == 450 and finite_close(wsum.to_numpy(), np.ones(450)), "fixed_weight_mass_one", f"regions={len(wsum)}; min={wsum.min()}; max={wsum.max()}", checks)
        source_mass = w.pivot_table(index="health_region_2023", columns="source_product", values="population_weight", aggfunc="sum", fill_value=0)
        source_mass = source_mass.reindex(sorted(x.health_region_2023.unique()))
        land_mass = source_mass["reanalysis-era5-land"].to_numpy()
        era_mass = source_mass["reanalysis-era5-single-levels"].to_numpy()
        static_mass = pd.DataFrame({"health_region_2023": source_mass.index, "land_static": land_mass, "era5_static": era_mass})
        x = x.merge(static_mass, on="health_region_2023", validate="many_to_one")
        check(finite_close(x.land_weight_total, x.land_static) and finite_close(x.era5_weight_total, x.era5_static), "daily_fixed_source_composition", "daily fixed source masses must match the source weight artifact", checks)
        check(finite_close(x.full_candidate_weight_total, np.ones(len(x))) and finite_close(x.land_weight_total + x.era5_weight_total, x.full_candidate_weight_total), "daily_fixed_weight_mass_one", "LAND and ERA5 single-level fixed masses must be disjoint complements totaling one", checks)
        mass_columns = ["full_candidate_weight_total", "land_weight_total", "era5_weight_total", "land_instantaneous_eligible_weight", "era5_instantaneous_eligible_weight"]
        mass_values = x[mass_columns].to_numpy(dtype=float)
        check(np.isfinite(mass_values).all() and (mass_values >= -TOL).all(), "finite_nonnegative_source_masses", "all fixed and recorded instantaneous source masses must be finite and nonnegative", checks)
        check(finite_close(x.land_instantaneous_eligible_weight, x.land_weight_total) and finite_close(x.era5_instantaneous_eligible_weight, x.era5_weight_total), "source_specific_eligible_mass_complete", "each recorded instantaneous source mass must equal its own fixed source mass", checks)
        eligible = x.land_instantaneous_eligible_weight + x.era5_instantaneous_eligible_weight
        check(np.isfinite(eligible).all() and finite_close(eligible, x.full_candidate_weight_total), "instantaneous_eligible_mass_complete", "source-specific instantaneous eligible weights must cover the full fixed mass; generic weight_coverage is not used", checks)
        check((x.land_instantaneous_eligible_weight <= x.land_weight_total + TOL).all() and (x.era5_instantaneous_eligible_weight <= x.era5_weight_total + TOL).all(), "eligible_support_within_source_mass", "eligible support cannot exceed assigned source weight", checks)

        for label, (primary_col, land_col, era_col, mixed_col) in metric_specs.items():
            vals = x[[primary_col, land_col, era_col, mixed_col]].to_numpy(dtype="float64")
            check(np.isfinite(vals).all(), f"finite:{label}", "all primary and source contribution values must be finite; no missing term may be filled with zero", checks)
            component_sum = x[land_col] + x[era_col]
            check(finite_close(component_sum, x[mixed_col]), f"component_sum:{label}", "mixed contribution must equal LAND plus ERA5 single-level contributions", checks)
            check(finite_close(component_sum, x[primary_col]), f"primary_sum:{label}", "primary weather field must equal the complete fixed-weight sum", checks)
        check(x["temperature_mean_c"].between(-100, 70).all() and x["temperature_max_c"].between(-100, 70).all(), "temperature_units_celsius", "temperature fields named in Celsius must be in plausible Celsius bounds", checks)
        check(x["relative_humidity_mean_pct"].between(0, 100).all(), "rh_units_percent", "relative humidity must be in [0,100] percent", checks)

        old = pq.read_table(panel_path)
        olddf = old.select(["hrid", "date", "tmean", "tmax", "rh", "t_p90", "t_p10", "heat", "cold", "smoke", "s_p75", "composite"]).to_pandas()
        olddf["_key_hr"] = olddf.hrid.astype(str)
        olddf["_key_date"] = pd.to_datetime(olddf.date).dt.normalize()
        xf = x[(x.date >= pd.Timestamp("2008-01-01")) & (x.date <= pd.Timestamp("2023-12-31"))].copy()
        xf = xf.rename(columns={"health_region_2023": "_key_hr", "date": "_key_date"})
        old_values = olddf[["_key_hr", "_key_date", "tmean", "tmax", "rh"]].rename(columns={"tmean": "tmean_old", "tmax": "tmax_old", "rh": "rh_old"})
        combined = old_values.merge(
            xf[["_key_hr", "_key_date", *[s[0] for s in metric_specs.values()], *[s[1] for s in metric_specs.values()], "land_instantaneous_eligible_weight"]],
            on=["_key_hr", "_key_date"], how="left", validate="one_to_one", indicator=True,
        )
        source_cols = [s[0] for s in metric_specs.values()] + [s[1] for s in metric_specs.values()] + ["land_instantaneous_eligible_weight"]
        source_complete = combined[source_cols].notna().to_numpy().all()
        check(len(old) == 2_629_800 and (combined["_merge"] == "both").all() and source_complete, "panel_source_key_completeness", f"panel_rows={len(old)}; matched={(combined['_merge'] == 'both').sum()}; source_fields_complete={source_complete}", checks)
        check(olddf["_key_hr"].nunique() == 450 and olddf["_key_date"].min() == pd.Timestamp("2008-01-01") and olddf["_key_date"].max() == pd.Timestamp("2023-12-31"), "panel_dates_regions", "expected 450 regions and daily 2008-01-01 through 2023-12-31", checks)
        # The old panel fields are direct LAND weighted contributions. LAND-normalized means
        # are retained only as a separate sensitivity and never used as the complete target.
        landden = combined["land_instantaneous_eligible_weight"]
        check((landden > 0).all(), "land_sensitivity_support", "LAND-normalized values require positive metric support", checks)
        for oldname, landcol in [
            ("tmean", "land_temperature_mean_c_weighted_contribution"),
            ("tmax", "land_temperature_max_c_weighted_contribution"),
        ]:
            check(finite_close(combined[f"{oldname}_old"], combined[landcol]), f"old_land_provenance:{oldname}", "old panel value must equal the retained LAND weighted contribution", checks)
        rhmask = combined["_key_date"].dt.year.between(2008, 2022)
        check(combined.loc[rhmask, "rh_old"].notna().all() and combined.loc[~rhmask, "rh_old"].isna().all(), "rh_source_window_mask", "preserve original RH-builder window: 2008–2022, leaving 2023 missing", checks)
        check(finite_close(combined.loc[rhmask, "rh_old"], combined.loc[rhmask, "land_relative_humidity_mean_pct_weighted_contribution"]), "old_land_provenance:rh", "old RH must equal the retained LAND weighted contribution in its source window", checks)

        # Confirm the retained threshold implementation before recalculating it.
        base = olddf[olddf["_key_date"].dt.year <= 2020]
        p90 = base.groupby("_key_hr").tmean.quantile(.90)
        p10 = base.groupby("_key_hr").tmean.quantile(.10)
        stored90 = olddf["_key_hr"].map(p90)
        stored10 = olddf["_key_hr"].map(p10)
        check(finite_close(olddf.t_p90, stored90) and finite_close(olddf.t_p10, stored10), "original_threshold_rule", "stored thresholds must match per-region pandas quantiles in the 2008–2020 reference window", checks)
        check((olddf.heat.to_numpy() == (olddf.tmean > olddf.t_p90).astype("int8").to_numpy()).all(), "original_heat_rule", "heat must use strict greater-than P90", checks)
        check((olddf.cold.to_numpy() == (olddf.tmean < olddf.t_p10).astype("int8").to_numpy()).all(), "original_cold_rule", "cold must use strict less-than P10", checks)
        check((olddf.composite.to_numpy() == (olddf.heat.astype(bool) & olddf.smoke.astype(bool)).astype("int8").to_numpy()).all(), "original_compound_rule", "composite must equal heat AND smoke", checks)

        report["support"] = {
            "static_land_mass_min": float(source_mass["reanalysis-era5-land"].min()),
            "static_land_mass_max": float(source_mass["reanalysis-era5-land"].max()),
            "static_era5_single_level_mass_min": float(source_mass["reanalysis-era5-single-levels"].min()),
            "static_era5_single_level_mass_max": float(source_mass["reanalysis-era5-single-levels"].max()),
            "land_eligible_mass_min": float(x.land_instantaneous_eligible_weight.min()),
            "land_eligible_mass_max": float(x.land_instantaneous_eligible_weight.max()),
            "era5_eligible_mass_min": float(x.era5_instantaneous_eligible_weight.min()),
            "era5_eligible_mass_max": float(x.era5_instantaneous_eligible_weight.max()),
            "eligibility_fields_used": ["land_instantaneous_eligible_weight", "era5_instantaneous_eligible_weight"],
            "coverage_fallback_used": False,
        }
        report["metric_checks"] = {name: {"primary_field": spec[0], "land_contribution": spec[1], "era5_single_level_contribution": spec[2], "mixed_contribution": spec[3], "units": "percent RH" if name == "rh" else "degrees Celsius"} for name, spec in metric_specs.items()}
        report["verification_scope"] = {
            "target": "mean over supplied fixed health-region/grid support, not an independently reconstructed population-weighted exposure",
            "validated": "input identity, recorded source mass and contribution assembly contract",
            "rh_shared_hourly_grid_mask_independently_verified": False,
            "rh_use": "descriptive only; preserve source limitation, no RH-adjusted inferential model admitted by this gate",
            "population_raster_production_lineage_verified": False,
        }
        report["land_normalized_sensitivity"] = {
            "available": True,
            "definition": "LAND contribution divided by LAND instantaneous eligible mass",
            "exported_as_separate_panel": False,
            "note": "Sensitivity panel export is deferred for root review; primary corrected input remains the complete fixed-weight sum.",
        }
        if args.validate_only:
            report["status"] = "validated_no_output_written"
        else:
            check(not output_path.exists(), "new_output_path", "refusing to overwrite an existing corrected input", checks)
            # Replace only exposure fields and directly derived temperature categories/thresholds.
            patched = old
            dates = olddf["_key_date"]
            tmean_new = olddf.tmean.copy()
            tmax_new = olddf.tmax.copy()
            rh_new = olddf.rh.copy()
            tmean_new.loc[:] = combined["temperature_mean_c"].to_numpy()
            tmax_new.loc[:] = combined["temperature_max_c"].to_numpy()
            rh_new.loc[rhmask] = combined.loc[rhmask, "relative_humidity_mean_pct"].to_numpy()
            newdf = olddf.copy()
            newdf["tmean"] = tmean_new
            newdf["tmax"] = tmax_new
            newdf["rh"] = rh_new
            ref = newdf[dates.dt.year <= 2020]
            newp90 = ref.groupby("_key_hr").tmean.quantile(.90)
            newp10 = ref.groupby("_key_hr").tmean.quantile(.10)
            newdf["t_p90"] = newdf["_key_hr"].map(newp90)
            newdf["t_p10"] = newdf["_key_hr"].map(newp10)
            newdf["heat"] = (newdf.tmean > newdf.t_p90).astype("int8")
            newdf["cold"] = (newdf.tmean < newdf.t_p10).astype("int8")
            newdf["composite"] = (newdf.heat.astype(bool) & newdf.smoke.astype(bool)).astype("int8")
            # Keep 2023 RH exactly as stored, including its missing mask.
            updates = {"tmean": newdf.tmean, "tmax": newdf.tmax, "rh": newdf.rh, "t_p90": newdf.t_p90, "t_p10": newdf.t_p10, "heat": newdf.heat, "cold": newdf.cold, "composite": newdf.composite}
            report["changed_row_counts"] = {
                name: int((~np.isclose(olddf[name].to_numpy(dtype="float64"), values.to_numpy(dtype="float64"), rtol=0, atol=0, equal_nan=True)).sum())
                for name, values in updates.items()
            }
            report["heat_and_compound_counts"] = {
                "before": {name: int(olddf[name].sum()) for name in ("heat", "cold", "composite")},
                "after": {name: int(newdf[name].sum()) for name in ("heat", "cold", "composite")},
            }
            for name, values in updates.items():
                i = patched.schema.get_field_index(name)
                arr = pa.array(values, type=patched.schema.field(i).type, from_pandas=True)
                patched = patched.set_column(i, name, arr)
            tmp = output_path.with_suffix(".parquet.partial")
            pq.write_table(patched, tmp, compression="zstd")
            os.replace(tmp, output_path)
            reread = pq.read_table(output_path)
            changed = set(updates)
            untouched = [name for name in old.column_names if name not in changed]
            untouched_ok = all(old[name].combine_chunks().equals(reread[name].combine_chunks()) for name in untouched)
            check(untouched_ok, "preserve_unrelated_fields", "all outcomes, dates, fire/smoke, and unrelated columns must equal the unchanged source panel", checks)
            report["output_path"] = str(output_path)
            report["output_sha256"] = sha256(output_path)
            report["unchanged_columns_verified"] = len(untouched)
            report["changed_fields"] = sorted(changed)
            report["status"] = "repaired_and_validated"
        report["row_counts"] = {"primary_region_days": len(p), "panel_rows": len(old)}
    except Exception as exc:
        report["status"] = "blocked"
        report["error"] = f"{type(exc).__name__}: {exc}"
    tmp_report = report_path.with_suffix(".json.partial")
    tmp_report.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp_report, report_path)
    print(json.dumps({"status": report["status"], "report": str(report_path), "output": report.get("output_path"), "error": report.get("error")}, ensure_ascii=False))
    return 0 if report["status"] in {"validated_no_output_written", "repaired_and_validated"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
