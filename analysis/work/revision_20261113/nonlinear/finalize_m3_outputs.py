"""Validate saved M3 fit checkpoints and finalize reports/status without refits."""
from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_m3_nonlinear as analysis


EXPECTED = [
    "brazil_strict_N20_linear",
    "brazil_strict_N20_spline3df",
    "mexico_strict_N20_linear",
    "mexico_strict_N20_spline3df",
    "brazil_primary_linear_temp_fire_two_way",
    "brazil_primary_spline3df_temp_fire_two_way",
    "brazil_primary_same_day_compound_two_way",
]


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> None:
    out = analysis.OUT
    previous_status = json.loads((out / "run_status.json").read_text(encoding="utf-8"))
    analysis.input_manifest()
    records: dict[str, dict] = {}
    for key in EXPECTED:
        path = out / f"fit_{key}.json"
        if not path.is_file():
            raise FileNotFoundError(f"Saved model checkpoint is missing: {path}")
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("source_hashes") != analysis.SOURCE_HASHES:
            raise ValueError(f"Source hash mismatch in checkpoint {key}")
        if record.get("spec_sha256") != analysis.SPEC_SHA256:
            raise ValueError(f"Fixed specification hash mismatch in checkpoint {key}")
        covariance = np.asarray(record["covariance"], dtype=float)
        if covariance.shape != (len(record["coefficient_names"]),) * 2:
            raise ValueError(f"Invalid full covariance dimensions in checkpoint {key}")
        records[key] = record
        analysis.FIT_RECORDS[key] = record

    linear = records["brazil_primary_linear_temp_fire_two_way"]
    compound = records["brazil_primary_same_day_compound_two_way"]
    if linear["model_rows"] != 1_481_520:
        raise ValueError(f"Brazil primary historical model N changed: {linear['model_rows']}")
    if compound["model_rows"] != 1_484_094:
        raise ValueError(f"Brazil same-day compound historical model N changed: {compound['model_rows']}")

    linear_beta = dict(zip(linear["coefficient_names"], linear["coefficients"]))
    temperature_rr = math.exp(sum(linear_beta[f"t{i}"] for i in range(8)))
    fire_rr = math.exp(sum(linear_beta[f"f{i}"] for i in range(8)))
    old_linear = pd.read_csv(analysis.OLD_LINEAR).query("outcome == 'N20_N23'").iloc[0]
    if not np.isclose(temperature_rr, float(old_linear.heat_RR_per_C), rtol=0, atol=1e-8):
        raise ValueError(f"Frozen temperature RR anchor failed: {temperature_rr}")
    if not np.isclose(fire_rr, float(old_linear.smoke_RR_per10), rtol=0, atol=1e-8):
        raise ValueError(f"Frozen fire RR anchor failed: {fire_rr}")

    compound_beta = dict(zip(compound["coefficient_names"], compound["coefficients"]))
    old_compound = pd.read_csv(analysis.OLD_COMPOUND).query("outcome == 'N20_N23'")
    compound_terms = {"heat_only": "heat_only", "smoke_only": "smoke_only", "composite": "both"}
    compound_anchors: dict[str, float] = {}
    for _, old in old_compound.iterrows():
        group = str(old.group)
        estimate = math.exp(compound_beta[compound_terms[group]])
        compound_anchors[group] = estimate
        if not np.isclose(estimate, float(old.RR), rtol=0, atol=1e-8):
            raise ValueError(f"Frozen compound point estimate anchor failed for {group}: {estimate}")

    source_audit = json.loads((out / "source_sample_audit.json").read_text(encoding="utf-8"))
    knot_info = json.loads((out / "knots_and_percentiles.json").read_text(encoding="utf-8"))
    primary_audit = source_audit["brazil_primary"]
    primary_knots = knot_info["Brazil_primary_2008_2022"]
    full_events = float(primary_audit["source_events_full_window"])
    lag_events = float(primary_audit["events_after_lags"])
    if full_events != 1_132_163 or lag_events != 1_130_926:
        raise ValueError(f"Brazil primary event flow did not reconcile: {full_events}, {lag_events}")
    primary_knots["events_full_window"] = full_events
    primary_knots["events_after_lags"] = lag_events
    analysis.save_json(out / "knots_and_percentiles.json", knot_info)

    fit_flow = pd.read_csv(out / "sample_flow.csv")
    if set(fit_flow.model) != set(EXPECTED) or len(fit_flow) != len(EXPECTED):
        raise ValueError("Exactly seven saved model sample-flow records were not found")
    raw_vs_psd_needed = [
        "psd_covariance_diagnostics.csv", "covariance_psd_sensitivity.json",
        "temperature_curves_psd.csv", "temperature_contrasts_psd.csv",
        "joint_nonlinearity_tests_psd.csv", "raw_vs_psd_curve_uncertainty.csv",
        "raw_vs_psd_percentile_uncertainty.csv",
    ]
    missing = [name for name in raw_vs_psd_needed if not (out / name).is_file()]
    if missing:
        raise FileNotFoundError("Explicit PSD sensitivity outputs missing: " + ", ".join(missing))

    analysis.write_method_report(
        source_audit["strict_country_frames"], primary_audit, knot_info,
        pd.read_csv(out / "joint_nonlinearity_tests.csv").to_dict(orient="records"),
        pd.read_csv(out / "historical_vs_updated_inference.csv"),
    )

    psd_tests = pd.read_csv(out / "joint_nonlinearity_tests_psd.csv")
    psd_diag = pd.read_csv(out / "psd_covariance_diagnostics.csv")
    p90 = pd.read_csv(out / "raw_vs_psd_percentile_uncertainty.csv")
    p90 = p90[p90.percentile.eq("P90")]
    psd_test_records = psd_tests[["family", "F", "df_num", "df_denom", "p_F", "raw_p_F", "raw_test_status"]].to_dict(orient="records")
    for row in psd_test_records:
        for key, value in list(row.items()):
            if pd.isna(value):
                row[key] = None
    check_summary = {
        "model_fit_count": len(records),
        "model_rows": {key: records[key]["model_rows"] for key in EXPECTED},
        "historical_primary_temperature_RR": temperature_rr,
        "historical_primary_fire_RR": fire_rr,
        "historical_same_day_compound_RRs": compound_anchors,
        "primary_events_full_window": full_events,
        "primary_events_after_lag_window": lag_events,
        "all_historical_sample_and_point_estimate_anchors_passed": True,
        "all_full_covariance_dimensions_valid": True,
        "explicit_psd_sensitivity_models": psd_diag.model.tolist(),
        "psd_joint_nonlinearity_F_and_p": psd_test_records,
        "P90_raw_vs_psd": p90.to_dict(orient="records"),
        "models_refit_during_finalization": False,
    }
    completed = now()
    analysis.save_json(out / "finalization_audit.json", {
        "status": "complete",
        "report_only_recovery": True,
        "matrix_only_psd_sensitivity": True,
        "initial_report_failure": previous_status.get("error"),
        "resolution": "Corrected report sample-flow field mapping; finalized report and PSD sensitivity from saved full coefficient/covariance checkpoints; no models refit.",
        "completed_utc": completed,
        "checks": check_summary,
    })
    metadata = {
        "status": "complete",
        "code_revision": analysis.CODE_REVISION,
        "specification_sha256": analysis.SPEC_SHA256,
        "completed_utc": completed,
        "threads": analysis.THREADS,
        "input_hashes": analysis.SOURCE_HASHES,
        "models_completed": EXPECTED,
        "historical_point_estimates_reproduced": True,
        "point_estimate_and_psd_checks": check_summary,
        "outputs": sorted(path.name for path in out.iterdir() if path.is_file()),
    }
    analysis.save_json(out / "analysis_metadata.json", metadata)
    status = {
        "status": "complete",
        "phase": "complete",
        "analysis_pid": previous_status.get("pid"),
        "analysis_process_running": False,
        "finalizer_pid": os.getpid(),
        "started_utc": previous_status.get("started_utc"),
        "completed_utc": completed,
        "elapsed_seconds_analysis": previous_status.get("elapsed_seconds"),
        "models_completed": EXPECTED,
        "log": str(out / "run.log"),
        "stderr": str(out / "run.err.log"),
        "manifest": str(out / "run_manifest.json"),
        "note": "Seven model fits and all historical anchors verified; report and explicit PSD sensitivity generated only from saved checkpoints.",
    }
    analysis.save_json(out / "run_status.json", status)
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    manifest.update({
        "status": "complete",
        "analysis_process_ended": True,
        "finalizer_pid": os.getpid(),
        "completed_utc": completed,
        "finalization_audit": str(out / "finalization_audit.json"),
    })
    analysis.save_json(out / "run_manifest.json", manifest)
    print(json.dumps({"status": "complete", "checks": check_summary}, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
