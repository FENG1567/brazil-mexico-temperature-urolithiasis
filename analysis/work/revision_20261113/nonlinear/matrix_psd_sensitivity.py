"""Explicit eigenvalue-zero PSD sensitivity from saved M3 spline fits only.

This script never refits a model. It preserves raw covariance artifacts and
writes a separate, fully labelled post-fit covariance sensitivity.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_m3_nonlinear as analysis


OUT = analysis.OUT
SPLINE_MODELS = {
    "brazil_strict_N20_spline3df": "Brazil",
    "mexico_strict_N20_spline3df": "Mexico",
    "brazil_primary_spline3df_temp_fire_two_way": "Brazil_primary_2008_2022",
}


def relative_negative_scale(eigenvalues: np.ndarray, *, numerical_tolerance: float = 0.0) -> float:
    if not len(eigenvalues):
        return 0.0
    scale = float(np.max(np.abs(eigenvalues)))
    minimum = float(np.min(eigenvalues))
    if scale == 0 or minimum >= -numerical_tolerance:
        return 0.0
    return abs(minimum) / scale


def project_eigenvalue_zero(covariance: np.ndarray) -> tuple[np.ndarray, dict[str, float | int]]:
    raw = (covariance + covariance.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(raw)
    clipped = np.maximum(eigenvalues, 0.0)
    projected = (eigenvectors * clipped) @ eigenvectors.T
    projected = (projected + projected.T) / 2.0
    projected_eigenvalues = np.linalg.eigvalsh(projected)
    tolerance = float(np.max(np.abs(clipped))) * 1e-12 if len(clipped) else 0.0
    difference = projected - raw
    raw_norm = float(np.linalg.norm(raw, ord="fro"))
    projected_min = float(projected_eigenvalues.min())
    projected_max = float(projected_eigenvalues.max())
    diagnostics = {
        "raw_minimum_eigenvalue": float(eigenvalues.min()),
        "raw_maximum_eigenvalue": float(eigenvalues.max()),
        "raw_negative_eigenvalue_count": int((eigenvalues < 0).sum()),
        "raw_relative_negative_eigenvalue_scale": relative_negative_scale(eigenvalues),
        "psd_minimum_eigenvalue": projected_min,
        "psd_maximum_eigenvalue": projected_max,
        "psd_negative_eigenvalue_count_beyond_tolerance": int((projected_eigenvalues < -tolerance).sum()),
        "psd_relative_negative_eigenvalue_scale": relative_negative_scale(projected_eigenvalues, numerical_tolerance=tolerance),
        "eigenvalue_zero_rule": "all negative eigenvalues set explicitly to zero; original covariance remains unchanged",
        "negative_eigenvalues_set_to_zero": int((eigenvalues < 0).sum()),
        "frobenius_correction_norm": float(np.linalg.norm(difference, ord="fro")),
        "relative_frobenius_correction": float(np.linalg.norm(difference, ord="fro") / raw_norm) if raw_norm else 0.0,
        "max_absolute_covariance_change": float(np.max(np.abs(difference))),
        "trace_change": float(np.trace(projected) - np.trace(raw)),
        "relative_trace_change": float((np.trace(projected) - np.trace(raw)) / abs(np.trace(raw))) if np.trace(raw) else 0.0,
        "max_absolute_diagonal_variance_change": float(np.max(np.abs(np.diag(projected) - np.diag(raw)))),
        "numerical_eigenvalue_tolerance": tolerance,
    }
    return projected, diagnostics


def main() -> None:
    analysis.input_manifest()
    expected_all = [
        "brazil_strict_N20_linear", "brazil_strict_N20_spline3df",
        "mexico_strict_N20_linear", "mexico_strict_N20_spline3df",
        "brazil_primary_linear_temp_fire_two_way",
        "brazil_primary_spline3df_temp_fire_two_way",
        "brazil_primary_same_day_compound_two_way",
    ]
    records: dict[str, dict] = {}
    for key in expected_all:
        checkpoint = OUT / f"fit_{key}.json"
        if not checkpoint.is_file():
            raise FileNotFoundError(f"Model checkpoint missing; no refit will be attempted: {checkpoint}")
        record = json.loads(checkpoint.read_text(encoding="utf-8"))
        if record.get("source_hashes") != analysis.SOURCE_HASHES or record.get("spec_sha256") != analysis.SPEC_SHA256:
            raise ValueError(f"Checkpoint source/spec hashes changed: {key}")
        covariance = np.asarray(record["covariance"], dtype=float)
        if covariance.shape != (len(record["coefficient_names"]),) * 2:
            raise ValueError(f"Saved coefficient covariance has invalid shape: {key}")
        analysis.FIT_RECORDS[key] = record
        records[key] = record

    # Refresh the raw diagnostic table using the corrected convention: positive
    # definite covariance has relative negative eigenvalue scale exactly zero.
    analysis.export_fit_tables()

    knot_info = json.loads((OUT / "knots_and_percentiles.json").read_text(encoding="utf-8"))
    sample_audit = json.loads((OUT / "source_sample_audit.json").read_text(encoding="utf-8"))
    primary_events_full = float(sample_audit["brazil_primary"]["source_events_full_window"])
    primary_events_after_lag = float(sample_audit["brazil_primary"]["events_after_lags"])
    knot_info["Brazil_primary_2008_2022"]["events_full_window"] = primary_events_full
    knot_info["Brazil_primary_2008_2022"]["events_after_lags"] = primary_events_after_lag
    (OUT / "knots_and_percentiles.json").write_text(
        json.dumps(knot_info, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )

    raw_test_rows = pd.read_csv(OUT / "joint_nonlinearity_tests.csv").to_dict(orient="records")
    raw_test_by_family = {str(row["family"]): row for row in raw_test_rows}
    raw_curves = pd.read_csv(OUT / "temperature_curves.csv")
    raw_contrasts = pd.read_csv(OUT / "temperature_and_category_contrasts.csv")

    psd_rows: list[dict] = []
    psd_covariances: dict[str, dict] = {}
    psd_test_rows: list[dict] = []
    psd_curve_rows: list[dict] = []
    psd_contrast_rows: list[dict] = []
    for model_key, family in SPLINE_MODELS.items():
        original = records[model_key]
        raw_covariance = np.asarray(original["covariance"], dtype=float)
        psd_covariance, diagnostic = project_eigenvalue_zero(raw_covariance)
        psd_record = copy.deepcopy(original)
        psd_record["covariance"] = psd_covariance.tolist()
        psd_record["covariance_version"] = "explicit eigenvalue-zero PSD sensitivity"
        diagnostic.update({"model": model_key, "family": family, "raw_covariance_overwritten": False,
                           "method": "symmetrize, eigendecompose full coefficient covariance, set every negative eigenvalue to zero, reconstruct"})
        psd_rows.append(diagnostic)
        psd_covariances[model_key] = {
            "coefficient_names": original["coefficient_names"],
            "coefficients": original["coefficients"],
            "covariance_psd_eigenzero": psd_covariance.tolist(),
            "diagnostics": diagnostic,
        }
        psd_test = analysis.nonlinear_wald_test(psd_record, family)
        raw_test = raw_test_by_family[family]
        psd_test.update({
            "covariance_version": "explicit full-matrix eigenvalue-zero PSD sensitivity",
            "raw_test_status": raw_test.get("status"),
            "raw_F": raw_test.get("F"),
            "raw_p_F": raw_test.get("p_F"),
            "raw_restriction_covariance_min_eigenvalue": raw_test.get("restriction_covariance_min_eigenvalue"),
            "raw_covariance_preserved": True,
        })
        psd_test_rows.append(psd_test)
        contrasts, curves = analysis.spline_curve(psd_record, knot_info[family], family, "temperature")
        for row in contrasts:
            row["covariance_version"] = "explicit full-matrix eigenvalue-zero PSD sensitivity"
        for row in curves:
            row["covariance_version"] = "explicit full-matrix eigenvalue-zero PSD sensitivity"
        psd_contrast_rows.extend(contrasts)
        psd_curve_rows.extend(curves)

    pd.DataFrame(psd_rows).to_csv(OUT / "psd_covariance_diagnostics.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(psd_test_rows).to_csv(OUT / "joint_nonlinearity_tests_psd.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(psd_curve_rows).to_csv(OUT / "temperature_curves_psd.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(psd_contrast_rows).to_csv(OUT / "temperature_contrasts_psd.csv", index=False, encoding="utf-8-sig")
    (OUT / "covariance_psd_sensitivity.json").write_text(
        json.dumps(psd_covariances, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )

    raw_curve_spline = raw_curves[raw_curves.model.isin(SPLINE_MODELS)].copy()
    psd_curve_frame = pd.DataFrame(psd_curve_rows)
    curve_compare = raw_curve_spline.merge(psd_curve_frame, on=["model", "curve_index"],
                                           suffixes=("_raw", "_psd"), validate="one_to_one")
    curve_compare_out = pd.DataFrame({
        "model": curve_compare.model,
        "curve_index": curve_compare.curve_index,
        "temperature_C": curve_compare.exposure_value_C_raw,
        "RR_raw": curve_compare.RR_raw,
        "RR_psd": curve_compare.RR_psd,
        "SE_raw": curve_compare.se_raw,
        "SE_psd": curve_compare.se_psd,
        "CI_low_raw": curve_compare.CI_low_raw,
        "CI_high_raw": curve_compare.CI_high_raw,
        "CI_low_psd": curve_compare.CI_low_psd,
        "CI_high_psd": curve_compare.CI_high_psd,
        "CI_width_raw": curve_compare.CI_high_raw - curve_compare.CI_low_raw,
        "CI_width_psd": curve_compare.CI_high_psd - curve_compare.CI_low_psd,
        "CI_width_change_psd_minus_raw": (curve_compare.CI_high_psd - curve_compare.CI_low_psd) - (curve_compare.CI_high_raw - curve_compare.CI_low_raw),
        "p_raw": curve_compare.p_t_raw,
        "p_psd": curve_compare.p_t_psd,
        "t_df": curve_compare.t_df_raw,
    })
    curve_compare_out.to_csv(OUT / "raw_vs_psd_curve_uncertainty.csv", index=False, encoding="utf-8-sig")

    raw_percentiles = raw_contrasts[
        raw_contrasts.model.isin(SPLINE_MODELS) & raw_contrasts.percentile.notna()
    ].copy()
    psd_percentiles = pd.DataFrame(psd_contrast_rows)
    percentile_compare = raw_percentiles.merge(psd_percentiles, on=["model", "percentile"],
                                                suffixes=("_raw", "_psd"), validate="one_to_one")
    percentile_compare_out = pd.DataFrame({
        "model": percentile_compare.model,
        "percentile": percentile_compare.percentile,
        "exposure_value_C": percentile_compare.exposure_value_C_raw,
        "RR_raw": percentile_compare.RR_raw,
        "RR_psd": percentile_compare.RR_psd,
        "SE_raw": percentile_compare.se_raw,
        "SE_psd": percentile_compare.se_psd,
        "CI_low_raw": percentile_compare.CI_low_raw,
        "CI_high_raw": percentile_compare.CI_high_raw,
        "CI_low_psd": percentile_compare.CI_low_psd,
        "CI_high_psd": percentile_compare.CI_high_psd,
        "p_raw": percentile_compare.p_t_raw,
        "p_psd": percentile_compare.p_t_psd,
        "t_df": percentile_compare.t_df_raw,
    })
    percentile_compare_out.to_csv(OUT / "raw_vs_psd_percentile_uncertainty.csv", index=False, encoding="utf-8-sig")

    print(json.dumps({
        "status": "complete",
        "models_used": list(SPLINE_MODELS),
        "model_fits_recomputed": False,
        "primary_events_full_window_corrected": primary_events_full,
        "primary_events_after_lag": primary_events_after_lag,
        "raw_joint_tests": raw_test_rows,
        "psd_joint_tests": psd_test_rows,
        "psd_covariance_diagnostics": psd_rows,
        "p90_raw_vs_psd": percentile_compare_out[percentile_compare_out.percentile.eq("P90")].to_dict(orient="records"),
        "outputs": [
            "psd_covariance_diagnostics.csv", "covariance_psd_sensitivity.json",
            "temperature_curves_psd.csv", "temperature_contrasts_psd.csv",
            "joint_nonlinearity_tests_psd.csv", "raw_vs_psd_curve_uncertainty.csv",
            "raw_vs_psd_percentile_uncertainty.csv",
        ],
    }, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
