from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[4]
BASE = Path(r"F:\analysis\climate_biliary_20261102\work")
WORK = ROOT / "work" / "revision_20261113" / "descriptive"
OUT = ROOT / "outputs" / "revision_20261113" / "descriptive"
UP = ROOT / "outputs" / "upgrade_20261112"
SERVICE = ROOT / "outputs" / "revision_20261113" / "service"
CLINICAL = ROOT / "work" / "upgrade_20261112" / "clinical_code_day_aggregate.parquet"
PANEL = BASE / "panel_composite_humidity.parquet"
AGE_PANEL = BASE / "panel_age_year_2008_2023.parquet"
MX_SOURCE = BASE / "mexico_n20_proced.parquet"
MX_ID_SOURCE = BASE / "mexico_n20_id_proced.parquet"
MX_MUNI_DAY = BASE / "mexico_n20_municipality_day.parquet"
MX_MUNI_PANEL = BASE / "mexico_n20_municipality_panel.parquet"
MX_STATE_PANEL = BASE / "mexico_composite_state_day_panel.parquet"
MX_STATE_TEMP = BASE / "mexico_state_day_temperature.parquet"

for directory in (WORK, OUT):
    directory.mkdir(parents=True, exist_ok=True)


def write_csv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")


def int_or_none(value: object) -> int | None:
    if pd.isna(value):
        return None
    return int(value)


def float_or_none(value: object) -> float | None:
    if pd.isna(value):
        return None
    return float(value)


def exposure_stats(frame: pd.DataFrame, country: str, panel_name: str,
                   year: int | None, variables: dict[str, str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    period = str(year) if year is not None else "all_years"
    for variable, unit in variables.items():
        if variable not in frame:
            continue
        values = pd.to_numeric(frame[variable], errors="coerce").to_numpy(dtype="float64")
        valid = values[np.isfinite(values)]
        if not len(valid):
            rows.append({"country": country, "panel": panel_name, "period": period,
                         "variable": variable, "unit": unit, "area_day_n": len(frame),
                         "nonmissing_n": 0, "missing_n": len(frame)})
            continue
        q = np.percentile(valid, [1, 5, 25, 50, 75, 90, 95, 99])
        rows.append({
            "country": country, "panel": panel_name, "period": period,
            "variable": variable, "unit": unit, "area_day_n": int(len(frame)),
            "nonmissing_n": int(len(valid)), "missing_n": int(len(frame) - len(valid)),
            "mean": float(valid.mean()), "sd": float(valid.std(ddof=1)) if len(valid) > 1 else 0.0,
            "min": float(valid.min()), "p01": float(q[0]), "p05": float(q[1]),
            "p25": float(q[2]), "median": float(q[3]), "p75": float(q[4]),
            "p90": float(q[5]), "p95": float(q[6]), "p99": float(q[7]),
            "max": float(valid.max()),
        })
    return rows


def add_category(frame: pd.DataFrame, heat_col: str = "heat",
                 smoke_col: str = "smoke") -> pd.Series:
    heat = pd.to_numeric(frame[heat_col], errors="coerce")
    smoke = pd.to_numeric(frame[smoke_col], errors="coerce")
    result = pd.Series("missing", index=frame.index, dtype="object")
    known = heat.isin([0, 1]) & smoke.isin([0, 1])
    result.loc[known & heat.eq(0) & smoke.eq(0)] = "neither"
    result.loc[known & heat.eq(1) & smoke.eq(0)] = "heat_only"
    result.loc[known & heat.eq(0) & smoke.eq(1)] = "smoke_only"
    result.loc[known & heat.eq(1) & smoke.eq(1)] = "both"
    return result


def category_support(frame: pd.DataFrame, country: str, area_col: str,
                     outcomes: list[str], panel_name: str,
                     valid_mask: pd.Series | None = None) -> pd.DataFrame:
    work = frame.copy()
    work["exposure_category"] = add_category(work)
    if valid_mask is not None:
        work = work.loc[valid_mask].copy()
    work = work.loc[work.exposure_category.ne("missing")]
    grouped = work.groupby("exposure_category", observed=False)
    rows = []
    for category in ["neither", "heat_only", "smoke_only", "both"]:
        if category not in grouped.groups:
            rows.append({"country": country, "panel": panel_name, "exposure_category": category,
                         "area_day_n": 0, "areas_with_exposure_n": 0})
            continue
        g = grouped.get_group(category)
        row: dict[str, object] = {
            "country": country, "panel": panel_name, "exposure_category": category,
            "area_day_n": int(len(g)),
            "areas_with_exposure_n": int(g[area_col].nunique(dropna=True)),
        }
        for outcome in outcomes:
            row[f"{outcome}_events"] = int(pd.to_numeric(g[outcome], errors="coerce").fillna(0).sum())
            row[f"{outcome}_positive_area_days"] = int(pd.to_numeric(g[outcome], errors="coerce").fillna(0).gt(0).sum())
        rows.append(row)
    result = pd.DataFrame(rows)
    result["area_day_pct"] = 100 * result.area_day_n / max(1, int(result.area_day_n.sum()))
    for outcome in outcomes:
        event_col = f"{outcome}_events"
        result[f"{outcome}_event_pct"] = 100 * result[event_col] / max(1, int(result[event_col].sum()))
    return result


def annual_diagnosis(frame: pd.DataFrame, country: str, date_col: str,
                     diagnoses: dict[str, str], date_group: pd.Series,
                     area_day_n: pd.Series | None = None) -> pd.DataFrame:
    rows = []
    years = sorted(pd.to_numeric(date_group, errors="coerce").dropna().astype(int).unique())
    for year in years:
        sel = date_group.eq(year)
        total = sum(int(pd.to_numeric(frame.loc[sel, col], errors="coerce").fillna(0).sum())
                    for col in diagnoses)
        for col, label in diagnoses.items():
            n = int(pd.to_numeric(frame.loc[sel, col], errors="coerce").fillna(0).sum())
            rows.append({"country": country, "year": int(year), "diagnosis": label,
                         "events": n, "share_of_listed_diagnoses_pct": 100 * n / total if total else np.nan,
                         "area_day_n": int(area_day_n.loc[sel].sum()) if area_day_n is not None else None})
    return pd.DataFrame(rows)


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    panel_cols = [
        "hrid", "date", "year", "pop", "K80", "K81", "K85", "N20_N23",
        "K80_sev", "K81_sev", "K85_sev", "N20_N23_sev",
        "K80_death", "K81_death", "K85_death", "N20_N23_death",
        "K80_icu", "K81_icu", "K85_icu", "N20_N23_icu",
        "K80_bd", "K81_bd", "K85_bd", "N20_N23_bd",
        "tmean", "tmax", "fire_pm25", "rh", "humidex", "heat", "smoke",
        "composite", "cold", "t_p90", "s_p75", "t_p10",
    ]
    panel = pd.read_parquet(PANEL, columns=panel_cols)
    panel["date"] = pd.to_datetime(panel.date).dt.normalize()
    panel["year"] = panel.date.dt.year.astype(int)
    panel = panel.loc[panel.year.between(2008, 2022)].copy()

    # Main Brazil panel: counts are hospitalization admissions aggregated to health-region/day.
    diagnosis_cols = {"K80": "K80", "K81": "K81", "K85": "K85", "N20_N23": "N20–N23"}
    annual = panel.groupby("year", as_index=False).agg(
        area_day_n=("date", "size"), health_regions=("hrid", "nunique"),
        **{f"{c}_events": (c, "sum") for c in diagnosis_cols},
    )
    write_csv(annual, "brazil_primary_annual.csv")
    diag_long = annual_diagnosis(
        panel, "Brazil", "date", diagnosis_cols, panel.year,
        pd.Series(1, index=panel.index),
    )
    write_csv(diag_long, "brazil_primary_annual_diagnosis_long.csv")

    # Population is area/year metadata repeated across days; deduplicate before summing.
    pop_cells = panel[["hrid", "year", "pop"]].drop_duplicates()
    pop_variation = pop_cells.groupby(["hrid", "year"]).pop.nunique(dropna=False)
    pop_year = pop_cells.groupby("year", as_index=False).agg(
        areas_with_population=("hrid", lambda x: int(x[pop_cells.loc[x.index, "pop"].notna()].nunique())),
        area_population_sum=("pop", "sum"), area_population_mean=("pop", "mean"),
        area_population_median=("pop", "median"), area_population_min=("pop", "min"),
        area_population_max=("pop", "max"), population_missing_area_years=("pop", lambda x: int(x.isna().sum())),
    )
    pop_year["regions_with_nonunique_population_within_year"] = int((pop_variation > 1).sum())
    write_csv(pop_year, "brazil_population_by_year.csv")

    severity_rows = []
    severity_vars = [c for c in panel.columns if c.endswith(("_sev", "_death", "_icu", "_bd"))]
    for year, g in panel.groupby("year", sort=True):
        for outcome in ["K80", "K81", "K85", "N20_N23"]:
            row = {"country": "Brazil", "year": int(year), "diagnosis": "N20–N23" if outcome == "N20_N23" else outcome}
            for suffix in ["_sev", "_death", "_icu", "_bd"]:
                col = outcome + suffix
                if col in severity_vars:
                    row[col] = float(pd.to_numeric(g[col], errors="coerce").fillna(0).sum())
            severity_rows.append(row)
    write_csv(pd.DataFrame(severity_rows), "brazil_severity_annual.csv")

    # Exposure distributions: area-day support; event support is tabulated separately.
    br_vars = {"tmean": "°C", "tmax": "°C", "fire_pm25": "µg/m³", "rh": "%", "humidex": "index"}
    exposure_rows = exposure_stats(panel, "Brazil", "450 health regions", None, br_vars)
    for year, g in panel.groupby("year", sort=True):
        exposure_rows.extend(exposure_stats(g, "Brazil", "450 health regions", int(year), br_vars))
    write_csv(pd.DataFrame(exposure_rows), "brazil_exposure_distribution.csv")
    br_support = category_support(panel, "Brazil", "hrid", ["N20_N23", "K80", "K81", "K85"], "450 health regions")
    write_csv(br_support, "brazil_exposure_category_support.csv")

    # Existing Brazilian age table supports age groups, but no sex field.
    age = pd.read_parquet(AGE_PANEL, columns=["year", "age_group", "K80", "K81", "K85", "N20_N23"])
    age["year"] = pd.to_numeric(age.year, errors="coerce")
    age = age.loc[age.year.between(2008, 2022)].copy()
    age_rows = []
    for outcome in ["K80", "K81", "K85", "N20_N23"]:
        g = age.groupby("age_group", as_index=False)[outcome].sum().rename(columns={outcome: "events"})
        denom = int(g.events.sum())
        g["country"] = "Brazil"
        g["period"] = "2008–2022"
        g["diagnosis"] = "N20–N23" if outcome == "N20_N23" else outcome
        g["share_within_diagnosis_pct"] = 100 * g.events / denom if denom else np.nan
        age_rows.append(g[["country", "period", "diagnosis", "age_group", "events", "share_within_diagnosis_pct"]])
    write_csv(pd.concat(age_rows, ignore_index=True), "brazil_age_group_distribution_cache.csv")

    # Brazil clinical aggregate: retain principal diagnosis and CAR_INT separately.
    clinical = pd.read_parquet(CLINICAL, columns=["IDENT", "diagnosis", "CAR_INT", "hrid", "date", "admissions"])
    clinical["date"] = pd.to_datetime(clinical.date).dt.normalize()
    clinical["year"] = clinical.date.dt.year
    clinical["diagnosis"] = clinical.diagnosis.astype(str).str.strip().str.upper()
    clinical["CAR_INT"] = clinical.CAR_INT.astype(str).str.strip().str.zfill(2)
    clinical = clinical.loc[
        clinical.year.between(2018, 2022) & clinical.IDENT.astype(str).eq("1")
        & clinical.diagnosis.isin(["N20", "N21", "N22", "N23"])
        & clinical.CAR_INT.isin(["01", "02"])
    ].copy()
    type_label = {"01": "elective", "02": "urgent"}
    clinical["admission_type"] = clinical.CAR_INT.map(type_label)
    clinical_annual = clinical.groupby(["year", "diagnosis", "CAR_INT", "admission_type"], as_index=False).agg(
        aggregate_cells=("admissions", "size"), admissions=("admissions", "sum"),
        health_regions=("hrid", "nunique"),
    )
    clinical_annual["diagnosis_scope"] = "strict diagnosis code"
    write_csv(clinical_annual, "brazil_clinical_annual_diagnosis_admission_type.csv")
    clinical_overall = clinical.groupby(["diagnosis", "CAR_INT", "admission_type"], as_index=False).agg(
        aggregate_cells=("admissions", "size"), admissions=("admissions", "sum"),
        health_regions=("hrid", "nunique"),
    )
    broad = clinical_overall.groupby(["CAR_INT", "admission_type"], as_index=False).agg(
        aggregate_cells=("aggregate_cells", "sum"), admissions=("admissions", "sum"),
        health_regions=("health_regions", "max"),
    )
    broad["diagnosis"] = "N20–N23 broad"
    strict = clinical_overall.loc[clinical_overall.diagnosis.eq("N20")].copy()
    strict["diagnosis"] = "N20 strict"
    clinical_summary = pd.concat([broad, strict], ignore_index=True)
    totals = clinical_summary.groupby("diagnosis").admissions.sum().to_dict()
    if int(totals.get("N20–N23 broad", -1)) != 431374:
        raise ValueError("Brazil broad 2018–2022 admission total does not match the saved source audit")
    if int(totals.get("N20 strict", -1)) != 321201:
        raise ValueError("Brazil strict N20 2018–2022 admission total does not match the saved source audit")
    write_csv(clinical_summary, "brazil_clinical_overall_diagnosis_admission_type.csv")

    # Mexico baseline cache has one de-identified row per admission, without exporting any person IDs.
    mx = pd.read_parquet(MX_SOURCE, columns=["INGRE", "ENTIDAD", "MUNIC", "AFECPRIN", "state", "muni", "date"])
    mx["date"] = pd.to_datetime(mx.date).dt.normalize()
    mx["year"] = mx.date.dt.year
    mx["diagnosis"] = mx.AFECPRIN.astype(str).str.strip().str.upper().str[:3]
    mx["state_num"] = pd.to_numeric(mx.ENTIDAD, errors="coerce")
    mx["valid_state"] = mx.state_num.between(1, 32)
    mx_target = mx.loc[mx.year.between(2018, 2022) & mx.diagnosis.isin(["N20", "N21", "N22", "N23"])].copy()
    mx_valid = mx_target.loc[mx_target.valid_state].copy()
    mx_annual = mx_valid.groupby(["year", "diagnosis"], as_index=False).size().rename(columns={"size": "admissions"})
    mx_annual["country"] = "Mexico"
    mx_annual["period"] = "2018–2022 baseline source"
    totals_by_year = mx_valid.groupby("year").size()
    mx_annual["annual_all_diagnoses"] = mx_annual.year.map(totals_by_year).astype(int)
    mx_annual["share_within_year_pct"] = 100 * mx_annual.admissions / mx_annual.annual_all_diagnoses
    write_csv(mx_annual[["country", "period", "year", "diagnosis", "admissions", "share_within_year_pct"]],
              "mexico_annual_diagnosis.csv")

    # Mexico's complete-case municipality branch is derived from the cached municipality-day and climate panels.
    md = pd.read_parquet(MX_MUNI_DAY, columns=["ENTIDAD", "key", "date"])
    md["date"] = pd.to_datetime(md.date).dt.normalize()
    md["state_num"] = pd.to_numeric(md.ENTIDAD, errors="coerce")
    md["valid_state"] = md.state_num.between(1, 32)
    md_day = md.groupby(["key", "date"], as_index=False).agg(
        source_admissions=("key", "size"), valid_state_admissions=("valid_state", "sum"),
        invalid_state_admissions=("valid_state", lambda x: int((~x).sum())),
    )
    mxc = pd.read_parquet(MX_MUNI_PANEL, columns=[
        "key", "date", "N20", "tmean", "fire_pm25", "rh", "heat", "smoke", "composite",
    ])
    mxc["date"] = pd.to_datetime(mxc.date).dt.normalize()
    mxc["year"] = mxc.date.dt.year
    mxc["exposure_category"] = add_category(mxc)
    mxc["climate_complete"] = mxc[["tmean", "fire_pm25", "heat", "smoke"]].notna().all(axis=1)
    mxc["climate_rh_complete"] = mxc[["tmean", "fire_pm25", "rh", "heat", "smoke"]].notna().all(axis=1)
    mxc = mxc.loc[mxc.year.between(2018, 2022)].copy()
    event_days = md_day.merge(mxc, on=["key", "date"], how="left", validate="one_to_one")
    event_days["N20"] = pd.to_numeric(event_days.N20, errors="coerce").fillna(0)
    recon = int((event_days.source_admissions != event_days.N20).sum())
    valid_state_cases = int(event_days.valid_state_admissions.sum())
    invalid_state_cases = int(event_days.invalid_state_admissions.sum())
    valid_complete = int(event_days.loc[event_days.climate_complete.fillna(False), "valid_state_admissions"].sum())
    valid_rh_complete = int(event_days.loc[event_days.climate_rh_complete.fillna(False), "valid_state_admissions"].sum())
    no_climate = int(valid_state_cases - valid_complete)
    municipality_n = int(event_days.loc[event_days.climate_complete.fillna(False) & event_days.valid_state_admissions.gt(0), "key"].nunique())
    mx_panel_total = int(mxc.N20.sum())
    mx_climate_total_allgeo = int(mxc.loc[mxc.climate_complete, "N20"].sum())
    mx_rh_total_allgeo = int(mxc.loc[mxc.climate_rh_complete, "N20"].sum())

    mx_exposure_rows = []
    valid_keys = set(md.loc[md.valid_state, "key"].dropna().astype(str).unique())
    mx_exp_valid = mxc.loc[mxc.key.astype(str).isin(valid_keys)].copy()
    mx_exp_complete = mx_exp_valid.loc[mx_exp_valid.climate_complete].copy()
    mx_exp_rh_complete = mx_exp_valid.loc[mx_exp_valid.climate_rh_complete].copy()
    mx_vars = {"tmean": "°C", "fire_pm25": "µg/m³", "rh": "%"}
    mx_exposure_rows.extend(exposure_stats(mx_exp_valid, "Mexico", "valid-state municipalities", None, mx_vars))
    for year, g in mx_exp_valid.groupby("year", sort=True):
        mx_exposure_rows.extend(exposure_stats(g, "Mexico", "valid-state municipalities", int(year), mx_vars))
    write_csv(pd.DataFrame(mx_exposure_rows), "mexico_municipality_exposure_distribution.csv")
    mx_support = category_support(mx_exp_valid, "Mexico", "key", ["N20"], "valid-state municipalities",
                                  valid_mask=mx_exp_valid.climate_complete)
    write_csv(mx_support, "mexico_municipality_exposure_category_support.csv")

    # The state-day panel has 32 × 1,826 fully populated dates and is distinct from municipality exposure support.
    mx_state = pd.read_parquet(MX_STATE_PANEL, columns=["state", "date", "N20", "tmean", "rh", "fire_pm25", "heat", "smoke", "composite"])
    mx_state["date"] = pd.to_datetime(mx_state.date).dt.normalize()
    mx_state["year"] = mx_state.date.dt.year
    mx_state = mx_state.loc[mx_state.year.between(2018, 2022)].copy()
    mx_state_support = category_support(mx_state, "Mexico", "state", ["N20"], "32-state daily panel")
    write_csv(mx_state_support, "mexico_state_exposure_category_support.csv")
    mx_state_year = mx_state.groupby("year", as_index=False).agg(
        state_day_n=("date", "size"), states=("state", "nunique"), N20_N23_events=("N20", "sum"),
    )
    write_csv(mx_state_year, "mexico_state_panel_annual.csv")

    # Use the previously audited 2023-source extension only as a separate boundary-flow row.
    boundary = read_json(UP / "mexico_boundary_audit.json")
    reconciliation = boundary["daily_count_reconciliation"]
    base_counts = reconciliation["base_counts_2018_2022_source_years"]
    additions = reconciliation["accepted_source_2023_additions"]
    augmented = reconciliation["augmented_counts_source_years_2018_2023"]
    if int(base_counts["broad_N20_N23"]) != int(mx_valid.shape[0]):
        raise ValueError("Mexico baseline valid-state cohort does not match saved source audit")
    if int(additions["broad_N20_N23"]) != 33 or int(augmented["broad_N20_N23"]) != int(base_counts["broad_N20_N23"]) + 33:
        raise ValueError("Mexico boundary source extension does not match saved source audit")

    # Reference sample flow: existing fitted outputs only; no association models are rerun here.
    service_flow = pd.read_csv(SERVICE / "sampleflow.csv")
    service_audit = pd.read_csv(SERVICE / "sample_audit.csv")
    service_events = pd.read_csv(SERVICE / "model_events.csv")
    write_csv(service_flow, "service_sampleflow_reference.csv")
    write_csv(service_audit, "service_sample_audit_reference.csv")
    write_csv(service_events, "service_model_events_reference.csv")

    # Fit-row attrition is reported separately from singleton observations. The calendar has
    # >=4 observations per stratum in the full window and >=3 after lag-0-7 restriction.
    br18 = panel.loc[panel.year.between(2018, 2022)].copy()
    br18["stratum"] = (br18.hrid.astype(str) + "_" + br18.date.dt.strftime("%Y%m") + "_" + br18.date.dt.dayofweek.astype(str))
    br_lag = br18.loc[br18.date.ge(pd.Timestamp("2018-01-08"))].copy()
    def size_summary(frame: pd.DataFrame) -> tuple[int, int, int]:
        sizes = frame.groupby("stratum").size()
        return int(sizes.size), int(sizes.min()), int((sizes == 1).sum())
    full_strata, full_min, full_singletons = size_summary(br18)
    lag_strata, lag_min, lag_singletons = size_summary(br_lag)
    model_rows = []
    for _, r in service_audit.iterrows():
        input_strata = int(r.strata_input)
        retained_strata = int(r.strata_model)
        if input_strata == full_strata:
            min_per_stratum, singleton_n = full_min, full_singletons
        elif input_strata == lag_strata:
            min_per_stratum, singleton_n = lag_min, lag_singletons
        else:
            min_per_stratum, singleton_n = None, None
        model_rows.append({
            "model_id": r.model_id, "model_family": r.model_family,
            "input_rows": int(r.input_rows), "retained_rows": int(r.model_rows),
            "rows_dropped_by_fit": int(r.rows_dropped_by_fit),
            "input_events": float(r.input_events), "retained_events": float(r.model_events),
            "input_base_strata": input_strata, "retained_base_strata": retained_strata,
            "base_strata_without_retained_rows": input_strata - retained_strata,
            "minimum_input_observations_per_base_stratum": min_per_stratum,
            "single_observation_base_strata": singleton_n,
            "note": "Saved fit audit; dropped rows are not labelled as singletons. Zero-event FE-cell counts are separately available in acute_repair_model_summary.csv for that analysis family.",
        })
    write_csv(pd.DataFrame(model_rows), "service_model_sample_flow.csv")

    acute_summary = pd.read_csv(UP / "acute_repair_model_summary.csv")
    acute_summary["single_observation_input_strata"] = 0
    acute_summary["singleton_interpretation"] = "The complete 2018–2022 health-region×year-month×weekday calendar has 4–5 rows per stratum; lag-0–7 input has ≥3. Single-observation strata are distinct from all-zero FE strata."
    write_csv(acute_summary, "acute_repair_conditional_poisson_sample_flow.csv")

    # Overall sample flow with original and exposure-complete denominators kept separate.
    mx_id_rows = int(pq.ParquetFile(MX_ID_SOURCE).metadata.num_rows)
    main_events = int(panel.N20_N23.sum())
    br_full_rows = int(len(br18))
    br_lag_rows = int(len(br_lag))
    sample_flow_rows = [
        {"country": "Brazil", "stage": "primary_panel_full_window", "period": "2008–2022", "unit": "health_region-day", "source_or_input_n": len(panel), "valid_geography_n": panel.hrid.nunique(), "exposure_complete_n": int(panel[["tmean", "fire_pm25", "rh"]].notna().all(axis=1).sum()), "lag0_7_complete_n": None, "events_broad_or_named": main_events, "events_urgent": None, "events_elective": None, "events_after_fit": None, "rows_after_fit": None, "note": "N20–N23 aggregate events in 450×calendar-day panel; not unique persons."},
        {"country": "Brazil", "stage": "clinical_aggregate_broad", "period": "2018–2022", "unit": "admission", "source_or_input_n": int(clinical.admissions.sum()), "valid_geography_n": int(clinical.loc[clinical.hrid.notna(), "admissions"].sum()), "exposure_complete_n": None, "lag0_7_complete_n": None, "events_broad_or_named": 431374, "events_urgent": 315504, "events_elective": 115870, "events_after_fit": None, "rows_after_fit": None, "note": "N20–N23, IDENT=1; CAR_INT 02 urgent and 01 elective; sum of aggregate admissions."},
        {"country": "Brazil", "stage": "clinical_aggregate_strict_N20", "period": "2018–2022", "unit": "admission", "source_or_input_n": 321201, "valid_geography_n": 321201, "exposure_complete_n": None, "lag0_7_complete_n": None, "events_broad_or_named": 321201, "events_urgent": 226402, "events_elective": 94799, "events_after_fit": None, "rows_after_fit": None, "note": "Strict diagnosis N20 only; separate from broad N20–N23."},
        {"country": "Brazil", "stage": "service_full_window_area_day", "period": "2018–2022", "unit": "health_region-day", "source_or_input_n": br_full_rows, "valid_geography_n": 450, "exposure_complete_n": br_full_rows, "lag0_7_complete_n": br_lag_rows, "events_broad_or_named": 226402, "events_urgent": 226402, "events_elective": None, "events_after_fit": 226402, "rows_after_fit": 366584, "note": "Saved service sample flow: strict urgent N20. Same-day full-window input then separate lag-0–7 complete input."},
        {"country": "Brazil", "stage": "SP_any_drainage_full_window", "period": "2018–2022", "unit": "procedure-linked admission", "source_or_input_n": 821700, "valid_geography_n": 450, "exposure_complete_n": 821700, "lag0_7_complete_n": 818550, "events_broad_or_named": 84764, "events_urgent": None, "events_elective": None, "events_after_fit": 84572, "rows_after_fit": 180241, "note": "84,764 full-window procedure events; 84,572 events in the lag-0–7 complete input/fit sample. See separate sample-flow stages."},
        {"country": "Brazil", "stage": "SP_any_removal_full_window", "period": "2018–2022", "unit": "procedure-linked admission", "source_or_input_n": 821700, "valid_geography_n": 450, "exposure_complete_n": 821700, "lag0_7_complete_n": 818550, "events_broad_or_named": 71833, "events_urgent": None, "events_elective": None, "events_after_fit": 71684, "rows_after_fit": 158873, "note": "71,833 full-window procedure events; continuous lag-0–7 model has 71,684 events and 158,873 fitted rows; categorical full-window model has 71,833 events and 159,591 fitted rows."},
        {"country": "Mexico", "stage": "baseline_raw_cache", "period": "2018–2022", "unit": "admission", "source_or_input_n": mx_id_rows, "valid_geography_n": int(mx_valid.shape[0]), "exposure_complete_n": valid_complete, "lag0_7_complete_n": None, "events_broad_or_named": int(mx_valid.shape[0]), "events_urgent": None, "events_elective": None, "events_after_fit": None, "rows_after_fit": None, "note": f"Raw cache {mx_id_rows}; invalid/out-of-range state {invalid_state_cases}; valid-state admissions {valid_state_cases}; complete climate exposure among valid-state admissions {valid_complete}."},
        {"country": "Mexico", "stage": "municipality_exposure_branch", "period": "2018–2022", "unit": "valid-state admission", "source_or_input_n": valid_state_cases, "valid_geography_n": valid_state_cases, "exposure_complete_n": valid_complete, "lag0_7_complete_n": None, "events_broad_or_named": valid_state_cases, "events_urgent": None, "events_elective": None, "events_after_fit": None, "rows_after_fit": None, "note": f"Valid-state admissions with complete temperature/fire categories: {valid_complete} across {municipality_n} municipalities; {no_climate} valid-state admissions lack complete municipality-date climate exposure."},
        {"country": "Mexico", "stage": "source_2023_boundary_additions", "period": "admission date in 2022", "unit": "admission", "source_or_input_n": int(additions["broad_N20_N23"]), "valid_geography_n": int(additions["broad_N20_N23"]), "exposure_complete_n": None, "lag0_7_complete_n": None, "events_broad_or_named": int(additions["broad_N20_N23"]), "events_urgent": None, "events_elective": None, "events_after_fit": None, "rows_after_fit": None, "note": "Separate source-boundary extension; 33 broad N20–N23, including 29 strict N20; no municipality-level exposure branch."},
    ]
    flow = pd.DataFrame(sample_flow_rows)
    write_csv(flow, "sample_flow_summary.csv")

    # Operational definitions from inspected builders/audits; no IDs or row-level data are written.
    definitions = [
        {"domain": "Brazil primary outcome", "period": "2008–2022", "definition": "Primary diagnosis DIAG_PRINC three-character prefix; N20, N21, N22, N23 are the broad N20–N23 outcome. K80, K81 and K85 are separate diagnosis groups.", "source_or_field": "panel_composite_humidity.parquet; build_panel_health_region.py", "limits": "Counts are hospitalization admissions summed to 450 health-region-days, not unique persons."},
        {"domain": "Brazil urgent/elective", "period": "2018–2022", "definition": "Clinical aggregate restricted to IDENT=1; CAR_INT=02 urgent and CAR_INT=01 elective. Broad N20–N23 and strict N20 are reported separately.", "source_or_field": "clinical_code_day_aggregate.parquet; repair_acute_phenotypes.py", "limits": "Urgent/elective are admission categories; procedure codes are not used to infer them."},
        {"domain": "Mexico primary outcome", "period": "2018–2022", "definition": "Principal diagnosis AFECPRIN prefix N20–N23 for broad outcome; prefix N20 for strict outcome; eligible geographic code is ENTIDAD 1–32.", "source_or_field": "mexico_n20_proced.parquet; mexico_boundary_audit.json", "limits": "Municipality exposure branch is a subset of the valid-state admission cohort."},
        {"domain": "Mexico urgent/elective", "period": "2018–2022", "definition": "Not available as a field in the retained Mexico admission/procedure caches; no urgent/elective classification is inferred from procedures.", "source_or_field": "mexico_n20_proced.parquet schema", "limits": "A separate Mexico admission-type description is unavailable in this cache."},
        {"domain": "Brazil age cache", "period": "2008–2022", "definition": "Existing aggregated age groups: 0–14, 15–44, 45–64 and 65+; produced from DATASUS age coding in build_age_counts.py.", "source_or_field": "panel_age_year_2008_2023.parquet; build_age_counts.py", "limits": "This legacy cache is K80-specific and has no sex field; the separate source audit reports age and raw SEXO source-code counts."},
        {"domain": "Brazil source age and sex audit", "period": "2008–2022", "definition": "Filtered RD aggregates use DIAG_PRINC, DT_INTER, MUNIC_RES, SEXO, COD_IDADE and IDADE. COD_IDADE=2/3 sub-year units are grouped under 0–14; code 4 uses IDADE in years; code 5 is grouped under 65+; code 0 is retained as unknown. Age and sex are grouped by year and diagnosis.", "source_or_field": "aggregate_case_demographics.py; DATASUS RD fields SEXO, COD_IDADE, IDADE", "limits": "Source counts reconcile exactly to the cached annual diagnosis panel. SEXO values remain source codes; no person or admission IDs were exported."},
        {"domain": "Mexico source age and sex audit", "period": "2018–2022", "definition": "Valid-state N20–N23 Egresos admissions are aggregated by SEXO, CVEEDAD and age category. Official SAEH mappings are CVEEDAD 0=hours, 1=days, 2=months, 3=years and 9=unspecified; age bands are assigned from EDAD only when CVEEDAD=3. Official SEXO labels are 1=male, 2=female, 3=intersex, 8=ignored and 9=unspecified.", "source_or_field": "aggregate_case_demographics.py; Egresos fields ENTIDAD, AFECPRIN, SEXO, EDAD, CVEEDAD", "limits": "Counts reconcile exactly to the valid-state baseline cache. CVEEDAD=4/5 are present in the raw data but not defined by the consulted codebooks, so they remain unmapped; EDAD=999 and unusual values receive no year-based age group. Source codes are retained and no person or admission IDs were exported."},
        {"domain": "Brazil ICU/death/bed-days", "period": "2008–2023 panel source; main table window 2008–2022", "definition": "Daily severity builder uses MORTE=1 for death, UTI_MES_TO>0 for ICU, and sums DIAS_PERM for bed-days by primary diagnosis prefix. A separate strict urgent N20 severe definition also treats UTI_INT_TO>0 as ICU and defines severe as ICU or death.", "source_or_field": "build_daily_severity.py; build_beddays.py; build_n20_severe.py", "limits": "Broad panel ICU and strict-N20 ICU definitions differ; keep their labels separate."},
        {"domain": "Negative controls: Brazil", "period": "2008–2023 source extraction", "definition": "Principal diagnosis prefixes H25/H26 for cataract; S02, S12, S22, S32, S42, S52, S62, S72, S82 and S92 for fracture.", "source_or_field": "build_negative_controls.py; DIAG_PRINC", "limits": "Three-character code groups; cataract is a negative-control outcome, not proof of no confounding."},
        {"domain": "Negative controls: Mexico", "period": "2018–2022", "definition": "AFECPRIN prefixes H25/H26 for cataract; S02, S12, S22, S32, S42, S52, S62, S72, S82 and S92 for fracture; valid states ENTIDAD 1–32.", "source_or_field": "extract_mexico_controls.py", "limits": "State-day aggregates; no individual rows exported."},
        {"domain": "Brazil heat/smoke categories", "period": "Threshold reference 2008–2020; application 2008–2022", "definition": "Per-health-region tmean P90 and fire PM2.5 P75 estimated in the 2008–2020 baseline; heat and smoke use strict greater-than; both means heat=1 and smoke=1.", "source_or_field": "build_composite_panel.py", "limits": "Four descriptive groups are neither, heat only, smoke only and both; composite is the both indicator."},
        {"domain": "Brazil humidity and humidex", "period": "2008–2022", "definition": "RH and humidex are merged by health_region_2023 and local_date from land-weighted contribution columns; the helper merges precomputed weighted contributions.", "source_or_field": "build_composite_humidity.py", "limits": "The current descriptive pass does not rebuild the upstream grid contribution weights."},
        {"domain": "Mexico state temperature", "period": "2018–2022", "definition": "ERA5-Land hourly 2 m temperature converted K to °C, averaged by day, then cos(latitude)-weighted over grid points within state polygons; 32 states × 1,826 days. NASA POWER at state-capital coordinates is a distinct secondary product.", "source_or_field": "mexico_strict_phenotype.py; build_mex_state_day.py", "limits": "This state-level method is distinct from the municipality exposure-complete branch."},
        {"domain": "Mexico heat/smoke categories", "period": "Threshold reference 2018–2020; application 2018–2022", "definition": "State-specific tmean P90 and fire PM2.5 P75, each estimated through 2020; strict greater-than thresholds; both is heat=1 and smoke=1.", "source_or_field": "build_mexico_composite_panel.py", "limits": "Municipality panel categories are summarized separately from the 32-state model panel."},
        {"domain": "Mexico municipality exposure completeness", "period": "2018–2022", "definition": "Admission is climate-complete when its municipality-date has nonmissing tmean, fire PM2.5, heat and smoke; RH-complete is separately tabulated.", "source_or_field": "mexico_n20_municipality_panel.parquet; mexico_n20_municipality_day.parquet", "limits": "1,098 exclusions are decomposed by invalid state versus valid-state missing exposure in sample_flow_summary.csv."},
        {"domain": "Capacity", "period": "CNES 2018–2022; tertile baseline 2018", "definition": "Monthly SUS beds from CNES LT, summed to 450 health regions; supply is SUS beds / panel population × 100,000. Fixed low/mid/high groups use the 2018 region-level mean tertiles.", "source_or_field": "build_capacity.py; capacity_report_zh.md", "limits": "Area-level supply proxy is not local hospital capacity or occupancy."},
        {"domain": "Sample-flow meaning", "period": "Brazil service and lag models", "definition": "Full-window area-day inputs and lag-0–7 complete inputs are separate. All-zero fixed-effect strata and one-observation strata are separate: the former are reported in acute_repair_model_summary.csv; the complete calendar has no one-observation base stratum (minimum 4 rows full-window, 3 after lag restriction).", "source_or_field": "service/sampleflow.csv; service/sample_audit.csv; acute_repair_model_summary.csv", "limits": "Rows dropped by a fit are never relabelled as singleton strata."},
    ]
    write_csv(pd.DataFrame(definitions), "methods_definitions.csv")

    # Aggregate JSON manifest; source identifiers and person-level values are omitted.
    row_total = int(len(panel))
    json_summary = {
        "status": "descriptive_cache_audit_complete",
        "scope": "aggregate-only; no association models rerun; no person IDs written",
        "brazil_primary": {
            "period": "2008-01-01 through 2022-12-31", "area_day_rows": row_total,
            "health_regions": int(panel.hrid.nunique()), "calendar_days": int(panel.date.nunique()),
            "N20_N23_admissions": main_events,
            "diagnosis_totals": {col: int(panel[col].sum()) for col in diagnosis_cols},
            "exposure_complete_area_days_tmean_fire_rh": int(panel[["tmean", "fire_pm25", "rh"]].notna().all(axis=1).sum()),
            "population_years": pop_year.to_dict(orient="records"),
            "age_group_cache": age_rows[0].to_dict(orient="records"),
            "sex_case_distribution": "See demographic_source_audit.json and demographic_summary_table.csv for source-code counts.",
            "age_group_source_distribution": "See demographic_source_audit.json and demographic_summary_table.csv for source-level age group counts.",
        },
        "brazil_clinical_2018_2022": {
            "broad_N20_N23_total": 431374, "broad_urgent": 315504, "broad_elective": 115870,
            "strict_N20_total": 321201, "strict_urgent": 226402, "strict_elective": 94799,
        },
        "mexico_baseline_2018_2022": {
            "raw_cache_records": mx_id_rows, "valid_state_admissions": valid_state_cases,
            "invalid_state_admissions": invalid_state_cases,
            "diagnosis_annual": mx_annual.to_dict(orient="records"),
            "municipality_exposure_complete_valid_state_admissions": valid_complete,
            "municipality_exposure_complete_with_rh_valid_state_admissions": valid_rh_complete,
            "valid_state_cases_without_climate_exposure": no_climate,
            "valid_state_municipalities_with_complete_exposure_support": municipality_n,
            "municipality_panel_rows": int(len(mxc)), "municipality_panel_cases_all_geographies": mx_panel_total,
            "climate_complete_cases_all_geographies": mx_climate_total_allgeo,
            "rh_complete_cases_all_geographies": mx_rh_total_allgeo,
            "event_day_join_mismatches": recon,
            "urgent_elective": "unavailable in current cache; not inferred from procedure",
            "age_sex": "See demographic_source_audit.json and demographic_summary_table.csv for codebook-labelled aggregate counts.",
            "source_2023_boundary_additions": {"broad_N20_N23": int(additions["broad_N20_N23"]), "strict_N20": int(additions["strict_N20"])},
        },
        "sample_flow": {
            "brazil_2018_2022_full_area_days": br_full_rows, "brazil_2018_2022_lag0_7_complete_area_days": br_lag_rows,
            "brazil_sp_drainage_full_window_events": 84764, "brazil_sp_drainage_lag0_7_events": 84572,
            "brazil_sp_removal_full_window_events": 71833, "brazil_sp_removal_lag0_7_events": 71684,
            "brazil_sp_removal_lag0_7_fitted_rows": 158873, "brazil_sp_removal_categorical_fitted_rows": 159591,
            "brazil_base_strata_full": full_strata, "brazil_min_rows_per_stratum_full": full_min,
            "brazil_single_observation_strata_full": full_singletons,
            "brazil_base_strata_lag": lag_strata, "brazil_min_rows_per_stratum_lag": lag_min,
            "brazil_single_observation_strata_lag": lag_singletons,
        },
        "source_files": {
            "brazil_primary_panel": PANEL.name, "brazil_clinical_aggregate": CLINICAL.name,
            "mexico_baseline_case_cache": MX_SOURCE.name, "mexico_id_cache_person_ids_not_exported": MX_ID_SOURCE.name,
            "mexico_municipality_exposure_panel": MX_MUNI_PANEL.name,
        },
    }
    (OUT / "descriptive_results.json").write_text(json.dumps(json_summary, ensure_ascii=False, indent=2), encoding="utf-8")

    main_line = f"2008–2022年巴西主面板包含{row_total:,}个卫生区日、{panel.hrid.nunique():,}个卫生区；宽口径N20–N23累计{main_events:,}次住院。"
    pop_text = f"巴西暴露面板提供450个卫生区的年度人口字段；墨西哥保留的州日和市镇日暴露面板没有人口字段。"
    age_text = "年龄/性别来源聚合与映射说明见demographic_source_audit.json及demographic_summary_table.csv；Brazil保留SEXO原始代码，Mexico仅对CVEEDAD=3按EDAD生成年龄组。"
    flow_text = f"2018–2022巴西面板全窗为{br_full_rows:,}个区日；构造完整0–7日滞后后为{br_lag_rows:,}个区日。SP引流全窗计数84,764、滞后完整计数84,572；SP取石/碎石相关记录全窗71,833、滞后完整模型71,684（连续模型15.9万拟合区日行；分类模型159,591行）。"
    mx_text = f"墨西哥基线缓存{mx_id_rows:,}条记录中，有效州码病例{valid_state_cases:,}，无效/越界州码{invalid_state_cases:,}；有效州病例中{valid_complete:,}条具备市镇日温度与火源暴露，覆盖{municipality_n:,}个市镇。"
    type_text = "急诊/择期严格分口径：巴西宽N20–N23为315,504/115,870（合计431,374）；严格N20为226,402/94,799（合计321,201）。"
    md = f"""# 描述性病例、暴露与样本流程审计

本目录只汇总既有聚合缓存、既有审计输出和已保存拟合样本信息；没有重新估计关联模型，没有写出个人标识或病例级记录。

## 主要病例与人口

{main_line}

`brazil_primary_annual.csv`列出年度诊断事件数；`brazil_population_by_year.csv`列出面板中可得的卫生区年度人口汇总；`brazil_age_group_distribution_cache.csv`为现有年龄组病例聚合。

{pop_text}

{age_text}

## 急诊/择期与墨西哥流程

{type_text}

{mx_text}

原始墨西哥基线记录与有效州样本是两个分母：40,907条缓存记录中492条州码无效或越界；有效州样本40,415条，其中市镇日温度/火源暴露完整39,809条，差额606条。两种原因互斥，故相对原始缓存总排除1,098条。附加来源的33条（含严格N20 29条）单列为边界敏感性记录，不并入市镇暴露分支。

## 暴露分布与类别支持

`brazil_exposure_distribution.csv`和`mexico_municipality_exposure_distribution.csv`按地区日暴露分布报告均值、标准差、分位数和缺失数。类别支持表逐列给出四组地区日、病例数和有病例日；墨西哥另列32州模型面板和市镇完整暴露分支，二者不可混作同一空间总体。

{flow_text}

`service_sampleflow_reference.csv`、`service_model_events_reference.csv`与`service_model_sample_flow.csv`保留已有拟合行与拟合事件流；`acute_repair_conditional_poisson_sample_flow.csv`单列全零固定效应层、相应删行数和拟合器额外删行数。完整日历的基线分层每层至少4行，lag-0–7输入每层至少3行；没有单观测行层。全零固定效应层与singleton不是同一概念，样本流中未将拟合删行统称为singleton。

## 方法定义

`methods_definitions.csv`列出病例口径、急择期、负对照ICD组、ICU/死亡/床日、RH、温度/火源阈值和容量分组的操作定义及解读边界。巴西主面板与墨西哥州级温度使用不同空间赋值方法；墨西哥主州级方法为ERA5-Land日均温度并按州内网格点`cos(latitude)`加权。巴西RH/Humidex来自已加权贡献字段，本次不重新构造上游网格权重。

以上输出不重新估计关联模型；年龄/性别来源计数仅为去标识聚合，具体年龄标签遵循来源代码本及现有字段定义。
"""
    (OUT / "descriptive_summary_zh.md").write_text(md, encoding="utf-8")

    # When a validated aggregate-only source audit is available, restore its detailed
    # age/sex tables and method notes after rebuilding the descriptive cache outputs.
    audit_path = OUT / "demographic_source_audit.json"
    joint_path = OUT / "mexico_demographics_joint.csv"
    if audit_path.exists() and joint_path.exists():
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        mexico_audit = audit.get("Mexico", {})
        if audit.get("status") == "complete" and mexico_audit.get("age_codebook_sources"):
            from finalize_descriptive_audit import finalize
            finalize()

    print(json.dumps({
        "status": "complete", "brazil_area_days": row_total,
        "brazil_N20_N23_2008_2022": main_events,
        "brazil_clinical_broad": 431374, "brazil_clinical_strict_N20": 321201,
        "mexico_raw": mx_id_rows, "mexico_valid_state": valid_state_cases,
        "mexico_complete_exposure": valid_complete, "mexico_exposure_missing_valid_state": no_climate,
        "mexico_exposure_municipalities": municipality_n, "mexico_join_mismatches": recon,
        "brazil_full_days_2018_2022": br_full_rows, "brazil_lag_days_2018_2022": br_lag_rows,
        "output_dir": str(OUT),
    }, ensure_ascii=False))


if __name__ == "__main__":
    sys.path.insert(0, str(BASE / "pylibs"))
    main()
