"""Finite M3 natural-spline sensitivity and Brazil primary inference update.

Writes all new outputs only to outputs/revision_20261113/nonlinear. The
specification is fixed in prespecified_specification.json before estimation.
"""
from __future__ import annotations

import gc
import hashlib
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

for _name in (
    "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_name] = "4"

ROOT = Path(__file__).resolve().parents[3]
WORK = ROOT / "work" / "revision_20261113" / "nonlinear"
OUT = ROOT / "outputs" / "revision_20261113" / "nonlinear"
SPEC_PATH = WORK / "prespecified_specification.json"
BRAZIL_PANEL = Path(r"F:\analysis\climate_biliary_20261102\work\panel_composite_humidity.parquet")
MEXICO_PANEL = Path(r"F:\analysis\climate_biliary_20261102\work\mexico_composite_state_day_panel.parquet")
MEXICO_ERA5 = Path(r"D:\数据库集合\新增数据\ERA5_crossnational\mex_state_day_2015_2023.csv")
MEXICO_DAILY = ROOT / "outputs" / "upgrade_20261112" / "mexico_boundary_daily_counts.csv"
HOLIDAY_CSV = ROOT / "outputs" / "upgrade_20261112" / "holiday_calendar_dates.csv"
OLD_LINEAR = ROOT / "work" / "frozen_snapshot" / "tables" / "n20_composite_continuous_lag.csv"
OLD_COMPOUND = ROOT / "work" / "frozen_snapshot" / "tables" / "n20_composite_groups_corrected.csv"
CODE_REVISION = "nonlinear-m3-20261113-v1"
THREADS = 4
STATUS_PATH = OUT / "run_status.json"
LOG_PATH = OUT / "run.log"

OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, r"F:/analysis/climate_biliary_20261102/work/pylibs")
sys.path.insert(0, str(ROOT / "work" / "upgrade_20261112"))

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import f as f_dist
from scipy.stats import t as student_t


T0 = time.time()
FIT_RECORDS: dict[str, dict[str, Any]] = {}
SOURCE_HASHES: dict[str, dict[str, Any]] = {}
SPEC_SHA256 = ""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (pd.Timestamp,)):
        return value.isoformat()
    if pd.isna(value):
        return None
    raise TypeError(f"Cannot JSON encode {type(value).__name__}")


def save_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_default, allow_nan=False),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def log(message: str, **values: Any) -> None:
    record = {"time": utc_now(), "elapsed_seconds": round(time.time() - T0, 2), **values}
    print(message, json.dumps(record, ensure_ascii=False, default=json_default), flush=True)


def write_status(status: str, phase: str, **values: Any) -> None:
    payload = {
        "status": status,
        "phase": phase,
        "pid": os.getpid(),
        "started_utc": datetime.fromtimestamp(T0, timezone.utc).isoformat(),
        "updated_utc": utc_now(),
        "elapsed_seconds": round(time.time() - T0, 2),
        "log": str(LOG_PATH),
        "protocol": str(SPEC_PATH),
        **values,
    }
    save_json(STATUS_PATH, payload)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_manifest() -> None:
    paths = [SPEC_PATH, BRAZIL_PANEL, MEXICO_PANEL, MEXICO_ERA5,
             MEXICO_DAILY, HOLIDAY_CSV, OLD_LINEAR, OLD_COMPOUND,
             ROOT / "work" / "upgrade_20261112" / "clinical_code_day_aggregate.parquet"]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required M3 inputs missing: " + "; ".join(missing))
    for path in paths:
        SOURCE_HASHES[str(path)] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    global SPEC_SHA256
    SPEC_SHA256 = SOURCE_HASHES[str(SPEC_PATH)]["sha256"]
    save_json(OUT / "input_manifest.json", SOURCE_HASHES)


def checkpoint_path(model_key: str) -> Path:
    return OUT / f"fit_{model_key}.json"


def relative_negative_eigenvalue_scale(eigenvalues: np.ndarray) -> float:
    scale = float(np.max(np.abs(eigenvalues))) if len(eigenvalues) else 0.0
    if scale == 0 or float(np.min(eigenvalues)) >= 0:
        return 0.0
    return float(abs(np.min(eigenvalues)) / scale)


def export_fit_tables() -> None:
    coefficient_rows: list[dict[str, Any]] = []
    diagnostic_rows: list[dict[str, Any]] = []
    flow_rows: list[dict[str, Any]] = []
    full_covariances: dict[str, Any] = {}
    for key, record in FIT_RECORDS.items():
        names = record["coefficient_names"]
        beta = np.asarray(record["coefficients"], dtype=float)
        covariance = np.asarray(record["covariance"], dtype=float)
        df = int(record["t_df"])
        critical = float(student_t.ppf(0.975, df))
        for index, name in enumerate(names):
            variance = float(covariance[index, index])
            if variance >= 0 and np.isfinite(variance):
                se = float(np.sqrt(variance))
                p = float(2 * student_t.sf(abs(beta[index] / se), df)) if se > 0 else 1.0
                lower = float(beta[index] - critical * se)
                upper = float(beta[index] + critical * se)
                rr = float(np.exp(beta[index]))
                rr_low = float(np.exp(lower))
                rr_high = float(np.exp(upper))
            else:
                se = p = lower = upper = rr = rr_low = rr_high = None
            coefficient_rows.append({
                "model": key, "term": name, "beta": float(beta[index]), "se": se,
                "RR": rr, "CI_low": rr_low, "CI_high": rr_high, "p_t": p,
                "t_df": df, "model_rows": record["model_rows"],
                "clusters_1": record["clusters_1"], "clusters_2": record["clusters_2"],
            })
        eig = np.asarray(record["eigenvalues"], dtype=float)
        diagnostic_rows.append({
            "model": key,
            "dimension": len(names),
            "max_abs_asymmetry": record["max_abs_asymmetry"],
            "minimum_eigenvalue": float(eig.min()),
            "maximum_eigenvalue": float(eig.max()),
            "negative_eigenvalue_count": int((eig < 0).sum()),
            "relative_negative_eigenvalue_scale": relative_negative_eigenvalue_scale(eig),
            "psd_projection_applied": False,
        })
        flow_rows.append(record["sample_flow"])
        full_covariances[key] = {
            "coefficient_names": names,
            "coefficients": record["coefficients"],
            "covariance": record["covariance"],
            "covariance_diagnostic": {
                "minimum_eigenvalue": float(eig.min()),
                "maximum_eigenvalue": float(eig.max()),
                "negative_eigenvalue_count": int((eig < 0).sum()),
                "relative_negative_eigenvalue_scale": relative_negative_eigenvalue_scale(eig),
                "max_abs_asymmetry": record["max_abs_asymmetry"],
                "psd_projection_applied": False,
            },
            "formula": record["formula"],
            "vcov": record["vcov"],
        }
    pd.DataFrame(coefficient_rows).to_csv(OUT / "coefficients.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(diagnostic_rows).to_csv(OUT / "covariance_diagnostics.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(flow_rows).to_csv(OUT / "sample_flow.csv", index=False, encoding="utf-8-sig")
    save_json(OUT / "covariance_matrices.json", full_covariances)


def load_or_fit(model_key: str, outcome: str, terms: list[str], frame: pd.DataFrame,
                fixed_effect: str, cluster_1: str, cluster_2: str,
                source_rows: int, full_window_events: float,
                model_family: str, outcome_label: str) -> dict[str, Any]:
    formula = f"{outcome} ~ " + " + ".join(terms) + f" | {fixed_effect}"
    cache = checkpoint_path(model_key)
    if cache.exists():
        try:
            old = json.loads(cache.read_text(encoding="utf-8"))
            if (old.get("code_revision") == CODE_REVISION
                    and old.get("spec_sha256") == SPEC_SHA256
                    and old.get("source_hashes") == SOURCE_HASHES
                    and old.get("formula") == formula
                    and old.get("analysis_rows") == len(frame)):
                FIT_RECORDS[model_key] = old
                log("CHECKPOINT_REUSED", model=model_key, model_rows=old["model_rows"])
                export_fit_tables()
                return old
        except Exception as exc:
            log("CHECKPOINT_IGNORED", model=model_key, error=repr(exc))

    log("FIT_START", model=model_key, input_rows=len(frame), formula=formula)
    write_status("running", f"fitting:{model_key}", current_model=model_key,
                 completed_models=list(FIT_RECORDS))
    started = time.time()
    required = list(dict.fromkeys([outcome, *terms, fixed_effect, cluster_1, cluster_2]))
    model_data = frame.loc[:, required].copy()
    if model_data.isna().any().any():
        raise ValueError(f"Unexpected missing model data in {model_key}")

    y_sum = model_data.groupby(fixed_effect, observed=True)[outcome].sum()
    fe_size = model_data.groupby(fixed_effect, observed=True).size()
    zero_levels = y_sum.index[y_sum.eq(0)]
    zero_rows = int(fe_size.reindex(zero_levels).sum())
    nonzero_levels = y_sum.index[y_sum.gt(0)]
    singleton_levels = fe_size.index[(fe_size.eq(1)) & fe_size.index.isin(nonzero_levels)]
    singleton_rows = int(fe_size.reindex(singleton_levels).sum())
    input_events = float(model_data[outcome].sum())
    input_rows = int(len(model_data))

    model = pf.fepois(
        formula, data=model_data,
        vcov={"CRV1": f"{cluster_1} + {cluster_2}"},
        copy_data=False, store_data=True,
    )
    names = [str(name) for name in model.coef().index]
    beta = model.coef().to_numpy(dtype=float)
    raw_covariance = np.asarray(model._vcov, dtype=float)
    if raw_covariance.shape != (len(names), len(names)) or not np.isfinite(raw_covariance).all():
        raise ValueError(f"Invalid covariance dimensions or values for {model_key}: {raw_covariance.shape}")
    max_asymmetry = float(np.max(np.abs(raw_covariance - raw_covariance.T)))
    covariance = (raw_covariance + raw_covariance.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(covariance)
    retained = model._data
    if retained is None or retained.empty:
        raise ValueError(f"PyFixest did not retain analyzed rows for {model_key}")
    g1 = int(retained[cluster_1].nunique())
    g2 = int(retained[cluster_2].nunique())
    g_intersection = int(retained[[cluster_1, cluster_2]].drop_duplicates().shape[0])
    t_df = min(g1 - 1, g2 - 1)
    if t_df < 1:
        raise ValueError(f"Insufficient retained clusters for {model_key}: {g1}, {g2}")
    model_rows = int(model._N)
    if len(retained) != model_rows:
        raise ValueError(f"Stored row count differs from estimated N for {model_key}")
    events_model = float(retained[outcome].sum())
    nonzero_rows = input_rows - zero_rows
    other_drops = input_rows - model_rows - zero_rows
    sample_flow = {
        "model": model_key,
        "family": model_family,
        "outcome": outcome_label,
        "source_rows_before_lag": int(source_rows),
        "analysis_rows_after_complete_lag_frame": input_rows,
        "source_events_full_window": float(full_window_events),
        "events_after_lag_before_fit": input_events,
        "zero_event_strata": int(len(zero_levels)),
        "zero_event_strata_rows_expected_removed": zero_rows,
        "rows_remaining_after_zero_strata": nonzero_rows,
        "positive_singleton_strata": int(len(singleton_levels)),
        "positive_singleton_rows_expected_removed": singleton_rows,
        "other_solver_dropped_rows_after_zero_strata": int(other_drops),
        "model_rows_retained": model_rows,
        "events_model_retained": events_model,
        "G_cluster_1_retained": g1,
        "G_cluster_2_retained": g2,
        "G_intersection_retained": g_intersection,
        "cluster_1": cluster_1,
        "cluster_2": cluster_2,
        "t_df": t_df,
    }
    record = {
        "model_key": model_key,
        "code_revision": CODE_REVISION,
        "spec_sha256": SPEC_SHA256,
        "source_hashes": SOURCE_HASHES,
        "formula": formula,
        "outcome": outcome,
        "outcome_label": outcome_label,
        "family": model_family,
        "fixed_effect": fixed_effect,
        "cluster_1": cluster_1,
        "cluster_2": cluster_2,
        "vcov": f"CRV1:{cluster_1}+{cluster_2}; inclusion-exclusion intersection retained",
        "coefficient_names": names,
        "coefficients": beta.tolist(),
        "covariance": covariance.tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "max_abs_asymmetry": max_asymmetry,
        "analysis_rows": input_rows,
        "model_rows": model_rows,
        "events_model": events_model,
        "clusters_1": g1,
        "clusters_2": g2,
        "clusters_intersection": g_intersection,
        "t_df": t_df,
        "fit_seconds": round(time.time() - started, 3),
        "sample_flow": sample_flow,
    }
    FIT_RECORDS[model_key] = record
    save_json(cache, record)
    del model, model_data
    gc.collect()
    export_fit_tables()
    log("FIT_COMPLETE", model=model_key, model_rows=model_rows, events=events_model,
        G1=g1, G2=g2, G_intersection=g_intersection,
        min_eigenvalue=float(eigenvalues.min()), negative_eigenvalues=int((eigenvalues < 0).sum()),
        seconds=record["fit_seconds"])
    return record


def estimate(record: dict[str, Any], weights: dict[str, float], estimand: str,
             reference: str | None = None, exposure_value: float | None = None) -> dict[str, Any]:
    names = record["coefficient_names"]
    beta = np.asarray(record["coefficients"], dtype=float)
    covariance = np.asarray(record["covariance"], dtype=float)
    vector = np.asarray([weights.get(name, 0.0) for name in names], dtype=float)
    if not any(name in weights for name in names):
        raise ValueError(f"No requested coefficient in {record['model_key']}: {estimand}")
    b = float(vector @ beta)
    variance = float(vector @ covariance @ vector)
    if variance < 0 or not np.isfinite(variance):
        return {
            "model": record["model_key"], "estimand": estimand,
            "reference": reference, "exposure_value_C": exposure_value,
            "beta": b, "variance": variance, "se": None, "RR": None,
            "CI_low": None, "CI_high": None, "p_t": None,
            "t_df": record["t_df"], "status": "negative_or_nonfinite_contrast_variance",
        }
    se = float(np.sqrt(variance))
    df = int(record["t_df"])
    critical = float(student_t.ppf(0.975, df))
    p = float(2 * student_t.sf(abs(b / se), df)) if se > 0 else 1.0
    return {
        "model": record["model_key"], "estimand": estimand,
        "reference": reference, "exposure_value_C": exposure_value,
        "beta": b, "variance": variance, "se": se,
        "RR": float(np.exp(b)), "CI_low": float(np.exp(b - critical * se)),
        "CI_high": float(np.exp(b + critical * se)), "p_t": p,
        "t_df": df, "status": "ok",
    }


def rcs_nonlinear(x: np.ndarray, knots: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Return the two restricted-cubic nonlinear columns for four fixed knots.

    The raw x column is supplied separately as the explicit linear component.
    The restricted-power basis has linear tails outside the boundary knots.
    """
    k = np.asarray(knots, dtype=float)
    if k.shape != (4,) or not np.all(np.diff(k) > 0):
        raise ValueError(f"Natural spline knots must be four strictly increasing values: {knots}")
    denom = float((k[3] - k[0]) ** 2)
    out: list[np.ndarray] = []
    for j in (0, 1):
        term = (
            np.maximum(x - k[j], 0.0) ** 3
            - np.maximum(x - k[2], 0.0) ** 3 * ((k[3] - k[j]) / (k[3] - k[2]))
            + np.maximum(x - k[3], 0.0) ** 3 * ((k[2] - k[j]) / (k[3] - k[2]))
        ) / denom
        out.append(term)
    return out[0], out[1]


def pooled_temperature_quantiles(frame: pd.DataFrame, lag_names: list[str]) -> dict[str, float]:
    values = frame[lag_names].to_numpy(dtype=float, copy=False).reshape(-1)
    if not np.isfinite(values).all():
        raise ValueError("The matched lag analysis frame contains missing/nonfinite temperature")
    probs = np.array([0.05, 1/3, 2/3, 0.95, 0.10, 0.50, 0.90])
    q = np.quantile(values, probs)
    knots = [float(q[0]), float(q[1]), float(q[2]), float(q[3])]
    if not np.all(np.diff(knots) > 0):
        raise ValueError(f"Temperature knot quantiles are not distinct: {knots}")
    return {
        "knots_P5_P33_333_P66_667_P95_C": knots,
        "P5_C": float(q[0]), "P10_C": float(q[4]), "P50_C": float(q[5]),
        "P90_C": float(q[6]), "P95_C": float(q[3]),
        "quantile_probabilities": probs.tolist(),
        "pooled_lag_temperature_values": int(values.size),
    }


def add_spline_terms(frame: pd.DataFrame, lag_names: list[str], knots: list[float]) -> list[str]:
    terms: list[str] = []
    for lag, column in enumerate(lag_names):
        values = frame[column].to_numpy(dtype=float, copy=False)
        nl1, nl2 = rcs_nonlinear(values, knots)
        col1, col2 = f"nl1_{lag}", f"nl2_{lag}"
        frame[col1] = nl1
        frame[col2] = nl2
        terms.extend([column, col1, col2])
        del values, nl1, nl2
    return terms


def linear_lag_estimate(record: dict[str, Any], prefix: str, label: str) -> dict[str, Any]:
    names = record["coefficient_names"]
    selected = [f"{prefix}{lag}" for lag in range(8)]
    if not set(selected).issubset(names):
        raise ValueError(f"Cumulative lag terms missing from {record['model_key']}")
    return estimate(record, {term: 1.0 for term in selected}, label)


def spline_curve(record: dict[str, Any], knot_info: dict[str, Any],
                 family: str, temperature_label: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    p50 = float(knot_info["P50_C"])
    qpoints = [("P10", float(knot_info["P10_C"])),
               ("P50", p50), ("P90", float(knot_info["P90_C"])),
               ("P95", float(knot_info["P95_C"]))]
    contrasts: list[dict[str, Any]] = []
    for percentile, value in qpoints:
        basis_value = [value, *rcs_nonlinear(np.asarray([value]), knot_info["knots_P5_P33_333_P66_667_P95_C"])]
        basis_ref = [p50, *rcs_nonlinear(np.asarray([p50]), knot_info["knots_P5_P33_333_P66_667_P95_C"])]
        weights: dict[str, float] = {}
        for lag in range(8):
            for index, component in enumerate(("linear", "nl1", "nl2")):
                weights[f"t{lag}" if component == "linear" else f"{component}_{lag}"] = float(basis_value[index][0] - basis_ref[index][0]) if index else float(basis_value[index] - basis_ref[index])
        row = estimate(record, weights, f"{temperature_label}_{percentile}_vs_P50", reference="P50", exposure_value=value)
        row.update({"percentile": percentile, "country_or_family": family})
        contrasts.append(row)

    curve_rows: list[dict[str, Any]] = []
    curve_x = np.linspace(float(knot_info["P5_C"]), float(knot_info["P95_C"]), 100)
    for i, value in enumerate(curve_x, start=1):
        basis_value = [float(value), *[float(z[0]) for z in rcs_nonlinear(np.asarray([value]), knot_info["knots_P5_P33_333_P66_667_P95_C"])] ]
        basis_ref = [p50, *[float(z[0]) for z in rcs_nonlinear(np.asarray([p50]), knot_info["knots_P5_P33_333_P66_667_P95_C"])] ]
        weights = {}
        for lag in range(8):
            weights[f"t{lag}"] = basis_value[0] - basis_ref[0]
            weights[f"nl1_{lag}"] = basis_value[1] - basis_ref[1]
            weights[f"nl2_{lag}"] = basis_value[2] - basis_ref[2]
        row = estimate(record, weights, f"{temperature_label}_curve", reference="P50", exposure_value=float(value))
        row.update({"curve_index": i, "country_or_family": family})
        curve_rows.append(row)
    return contrasts, curve_rows


def nonlinear_wald_test(record: dict[str, Any], family: str) -> dict[str, Any]:
    names = record["coefficient_names"]
    selected = [name for name in names if name.startswith("nl1_") or name.startswith("nl2_")]
    if len(selected) != 16:
        raise ValueError(f"Expected 16 nonlinear terms in {record['model_key']}, found {len(selected)}")
    indices = [names.index(name) for name in selected]
    beta = np.asarray(record["coefficients"], dtype=float)[indices]
    covariance = np.asarray(record["covariance"], dtype=float)[np.ix_(indices, indices)]
    eig = np.linalg.eigvalsh((covariance + covariance.T) / 2.0)
    rank = int(np.linalg.matrix_rank(covariance))
    if rank < 1:
        return {"model": record["model_key"], "family": family, "restrictions": 16,
                "rank": rank, "F": None, "df_num": None, "df_denom": record["t_df"],
                "p_F": None, "status": "singular_restriction_covariance"}
    statistic = float(beta @ np.linalg.pinv(covariance) @ beta / rank)
    if statistic < 0 or not np.isfinite(statistic) or eig.min() < -1e-10 * max(1.0, float(np.max(np.abs(eig)))):
        return {"model": record["model_key"], "family": family, "restrictions": 16,
                "rank": rank, "F": statistic if np.isfinite(statistic) else None,
                "df_num": rank, "df_denom": record["t_df"], "p_F": None,
                "restriction_covariance_min_eigenvalue": float(eig.min()),
                "status": "non_psd_restriction_covariance_no_projection"}
    return {
        "model": record["model_key"], "family": family,
        "null": "all 16 added nonlinear coefficients are zero, conditional on temperature linear terms and unchanged covariates",
        "restrictions": 16, "rank": rank, "F": statistic,
        "df_num": rank, "df_denom": int(record["t_df"]),
        "p_F": float(f_dist.sf(statistic, rank, int(record["t_df"]))),
        "restriction_covariance_min_eigenvalue": float(eig.min()),
        "status": "ok",
    }


def load_holidays() -> dict[str, set[pd.Timestamp]]:
    calendar = pd.read_csv(HOLIDAY_CSV, dtype={"country": "string"})
    if not {"country", "date"}.issubset(calendar.columns):
        raise ValueError("Holiday calendar lacks country/date columns")
    output: dict[str, set[pd.Timestamp]] = {}
    for country in ("BR", "MX"):
        dates = pd.to_datetime(calendar.loc[calendar.country.eq(country), "date"], errors="raise").dt.normalize()
        if dates.empty or dates.duplicated().any():
            raise ValueError(f"Invalid existing {country} public-date calendar")
        output[country] = set(dates)
    return output


def prepare_strict_frames(holidays: dict[str, set[pd.Timestamp]]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    import mexico_strict_phenotype as strict

    # The established loader is reused, with output redirected before the load call.
    strict.OUT = OUT
    strict.db.execute(f"SET threads={THREADS}")
    for variable in (
        "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
        "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = str(THREADS)
    brazil = strict.load_brazil()
    brazil["temperature"] = pd.to_numeric(brazil.tmean, errors="raise").astype(float)
    brazil["strict_N20"] = pd.to_numeric(brazil.strict_N20, errors="raise").astype("int64")
    brazil["day_id"] = brazil.date.dt.strftime("%Y%m%d")
    brazil["holiday"] = brazil.date.isin(holidays["BR"]).astype("int8")
    if brazil.duplicated(["state", "date"]).any():
        raise ValueError("Brazil strict state-day loader produced duplicates")

    connection = duckdb.connect()
    connection.execute(f"SET threads={THREADS}")
    mex_panel = connection.execute(
        "SELECT state, CAST(date AS DATE) AS date, stratum "
        "FROM read_parquet(?) WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'",
        [str(MEXICO_PANEL)],
    ).df()
    connection.close()
    mex_panel["state"] = mex_panel.state.astype(str).str.zfill(2)
    mex_panel["date"] = pd.to_datetime(mex_panel.date).dt.normalize()
    mex_panel["stratum"] = mex_panel.stratum.astype(str)
    strict.assert_complete_state_calendar(mex_panel, "Mexico", 32)

    mex_daily = pd.read_csv(MEXICO_DAILY, dtype={"state": "string"})
    required = {"state", "date", "admissions"}
    if not required.issubset(mex_daily.columns):
        raise ValueError(f"Mexico expanded daily input missing {sorted(required - set(mex_daily.columns))}")
    mex_daily["state"] = mex_daily.state.astype(str).str.zfill(2)
    mex_daily["date"] = pd.to_datetime(mex_daily.date, errors="raise").dt.normalize()
    mex_daily["strict_N20"] = pd.to_numeric(mex_daily.admissions, errors="raise").astype("int64")
    mex_daily = mex_daily[["state", "date", "strict_N20"]]
    if mex_daily.duplicated(["state", "date"]).any():
        raise ValueError("Mexico expanded daily strict N20 input has duplicate state-days")
    mexico = mex_panel.merge(mex_daily, on=["state", "date"], how="left", validate="one_to_one")
    if mexico.strict_N20.isna().any():
        raise ValueError("Mexico expanded admissions did not cover every state-day")
    mexico["strict_N20"] = mexico.strict_N20.astype("int64")

    era5 = pd.read_csv(MEXICO_ERA5, dtype={"state": "string"})
    needed_era5 = {"state", "year", "month", "day", "tmean"}
    if not needed_era5.issubset(era5.columns):
        raise ValueError(f"Mexico ERA5 input missing {sorted(needed_era5 - set(era5.columns))}")
    era5["state"] = era5.state.astype(str).str.zfill(2)
    era5["date"] = pd.to_datetime({
        "year": pd.to_numeric(era5.year, errors="raise").astype(int),
        "month": pd.to_numeric(era5.month, errors="raise").astype(int),
        "day": pd.to_numeric(era5.day, errors="raise").astype(int),
    }, errors="raise").dt.normalize()
    era5["temperature"] = pd.to_numeric(era5.tmean, errors="raise").astype(float)
    era5 = era5[["state", "date", "temperature"]]
    if era5.duplicated(["state", "date"]).any():
        raise ValueError("Mexico ERA5 input has duplicate state-days")
    mexico = mexico.merge(era5, on=["state", "date"], how="left", validate="one_to_one")
    if mexico.temperature.isna().any() or len(mexico) != 32 * 1826:
        raise ValueError("Mexico ERA5 does not fully cover expanded strict daily frame")
    mexico["day_id"] = mexico.date.dt.strftime("%Y%m%d")
    mexico["holiday"] = mexico.date.isin(holidays["MX"]).astype("int8")
    if mexico.duplicated(["state", "date"]).any():
        raise ValueError("Mexico strict state-day analysis frame has duplicates")

    source_audit = {
        "Brazil": {
            "outcome_definition": "strict N20; state-day aggregation from established population-weighted health-region loader",
            "source_rows": len(brazil), "state_clusters_pre_lag": int(brazil.state.nunique()),
            "dates_pre_lag": int(brazil.date.nunique()), "events_full_window": int(brazil.strict_N20.sum()),
            "temperature": "existing Brazil ERA5-Land state temperature assignment",
            "holiday_dates": len(holidays["BR"]),
        },
        "Mexico": {
            "outcome_definition": "strict N20 = admissions from audited expanded daily counts, source year includes 2023 rows dated through 2022",
            "source_rows": len(mexico), "state_clusters_pre_lag": int(mexico.state.nunique()),
            "dates_pre_lag": int(mexico.date.nunique()), "events_full_window": int(mexico.strict_N20.sum()),
            "temperature": "ERA5 cross-national state/day tmean Celsius",
            "holiday_dates": len(holidays["MX"]),
        },
    }
    return brazil, mexico, source_audit


def prepare_strict_lags(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame[["state", "date", "day_id", "stratum", "strict_N20", "temperature", "holiday"]].copy()
    out = out.sort_values(["state", "date"]).reset_index(drop=True)
    for lag in range(8):
        out[f"t{lag}"] = out.groupby("state", sort=False).temperature.shift(lag)
    lag_terms = [f"t{lag}" for lag in range(8)]
    out = out.dropna(subset=lag_terms).copy()
    out["state"] = out.state.astype("category")
    out["day_id"] = out.day_id.astype("category")
    out["stratum"] = out.stratum.astype("category")
    return out


def fit_strict_country(country: str, source_frame: pd.DataFrame,
                       source_audit: dict[str, Any],
                       curve_rows: list[dict[str, Any]],
                       contrast_rows: list[dict[str, Any]],
                       test_rows: list[dict[str, Any]],
                       knot_rows: dict[str, Any]) -> None:
    lagged = prepare_strict_lags(source_frame)
    t_names = [f"t{lag}" for lag in range(8)]
    knot_info = pooled_temperature_quantiles(lagged, t_names)
    knot_info.update({"country": country, "analysis_rows_after_lags": len(lagged),
                      "events_after_lags": int(lagged.strict_N20.sum()),
                      "events_full_window": int(source_audit["events_full_window"])})
    knot_rows[country] = knot_info
    fixed = "stratum"
    c1, c2 = "state", "day_id"
    linear = load_or_fit(
        f"{country.lower()}_strict_N20_linear", "strict_N20", t_names + ["holiday"],
        lagged, fixed, c1, c2, int(source_audit["source_rows"]),
        float(source_audit["events_full_window"]), "strict_N20", f"{country} strict N20",
    )
    contrast_rows.append(linear_lag_estimate(linear, "t", "cumulative_temperature_lag0_7_per_C"))
    nl_terms = add_spline_terms(lagged, t_names, knot_info["knots_P5_P33_333_P66_667_P95_C"])
    spline = load_or_fit(
        f"{country.lower()}_strict_N20_spline3df", "strict_N20", nl_terms + ["holiday"],
        lagged, fixed, c1, c2, int(source_audit["source_rows"]),
        float(source_audit["events_full_window"]), "strict_N20", f"{country} strict N20",
    )
    test_rows.append(nonlinear_wald_test(spline, country))
    crows, curves = spline_curve(spline, knot_info, country, "temperature")
    contrast_rows.extend(crows)
    curve_rows.extend(curves)
    del linear, spline, lagged
    gc.collect()


def prepare_brazil_primary() -> tuple[pd.DataFrame, dict[str, Any]]:
    con = duckdb.connect()
    con.execute(f"SET threads={THREADS}")
    frame = con.execute(
        "SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date, "
        "CAST(N20_N23 AS DOUBLE) AS N20_N23, CAST(tmean AS DOUBLE) AS temperature, "
        "CAST(fire_pm25 AS DOUBLE)/10.0 AS fire10, CAST(heat AS INTEGER) AS heat, "
        "CAST(smoke AS INTEGER) AS smoke "
        "FROM read_parquet(?) WHERE date >= DATE '2008-01-01' AND date < DATE '2023-01-01' "
        "ORDER BY hrid,date", [str(BRAZIL_PANEL)],
    ).df()
    con.close()
    frame["hrid"] = frame.hrid.astype(str)
    frame["date"] = pd.to_datetime(frame.date, errors="raise").dt.normalize()
    frame["ym"] = frame.date.dt.strftime("%Y%m")
    frame["dow"] = frame.date.dt.day_of_week.astype("int8")
    frame["stratum"] = frame.hrid + "_" + frame.ym + "_" + frame.dow.astype(str)
    frame["day_id"] = frame.date.dt.strftime("%Y%m%d")
    if frame.duplicated(["hrid", "date"]).any():
        raise ValueError("Brazil primary regional panel has duplicate region-days")
    if len(frame) != 450 * 5479 or frame.hrid.nunique() != 450:
        raise ValueError(f"Brazil primary 2008-2022 frame changed: {len(frame)} rows, {frame.hrid.nunique()} regions")
    if not frame.groupby("hrid", sort=False).date.diff().dropna().dt.days.eq(1).all():
        raise ValueError("Brazil primary 2008-2022 calendar has gaps")
    missing = frame[["N20_N23", "temperature", "fire10", "heat", "smoke"]].isna().sum().to_dict()
    # No row is removed for temperature percentiles; fire inputs are complete through 2022.
    if any(int(v) for v in missing.values()):
        raise ValueError(f"Unexpected missing Brazil primary input fields: {missing}")
    source_events_full_window = float(frame.N20_N23.sum())
    for lag in range(8):
        frame[f"t{lag}"] = frame.groupby("hrid", sort=False).temperature.shift(lag)
        frame[f"f{lag}"] = frame.groupby("hrid", sort=False).fire10.shift(lag)
    lag_terms = [f"t{lag}" for lag in range(8)] + [f"f{lag}" for lag in range(8)]
    before_lag = len(frame)
    frame = frame.dropna(subset=lag_terms).copy()
    if frame.empty:
        raise ValueError("Brazil primary matched lag frame is empty")
    frame["hrid"] = frame.hrid.astype("category")
    frame["day_id"] = frame.day_id.astype("category")
    frame["stratum"] = frame.stratum.astype("category")
    audit = {
        "source_rows_before_lag": before_lag,
        "source_events_full_window": source_events_full_window,
        "regions": int(frame.hrid.nunique()), "calendar_dates_after_lag": int(frame.day_id.nunique()),
        "analysis_rows_after_complete_lags": len(frame),
        "events_after_lags": float(frame.N20_N23.sum()),
        "lag_window": "filtered to 2008-01-01 through 2022-12-31 before lag creation; lags 0-7 shifted within hrid; first 7 days per hrid excluded",
    }
    return frame, audit


def fit_brazil_primary(curve_rows: list[dict[str, Any]], contrast_rows: list[dict[str, Any]],
                       test_rows: list[dict[str, Any]], knot_rows: dict[str, Any],
                       primary_audit: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    frame, lag_audit = prepare_brazil_primary()
    primary_audit.update(lag_audit)
    t_names, f_names = [f"t{i}" for i in range(8)], [f"f{i}" for i in range(8)]
    knot_info = pooled_temperature_quantiles(frame, t_names)
    knot_info.update({"country": "Brazil_primary_2008_2022", "analysis_rows_after_lags": len(frame),
                      "events_after_lags": float(frame.N20_N23.sum()),
                      "events_full_window": float(lag_audit["source_events_full_window"])})
    knot_rows["Brazil_primary_2008_2022"] = knot_info
    linear = load_or_fit(
        "brazil_primary_linear_temp_fire_two_way", "N20_N23", t_names + f_names,
        frame, "stratum", "hrid", "day_id", int(lag_audit["source_rows_before_lag"]),
        float(lag_audit["source_events_full_window"]), "Brazil_primary", "Brazil N20-N23 primary",
    )
    temp = linear_lag_estimate(linear, "t", "cumulative_temperature_lag0_7_per_C")
    fire = linear_lag_estimate(linear, "f", "cumulative_fire_lag0_7_per_10ug_m3")
    contrast_rows.extend([temp, fire])
    if linear["model_rows"] != 1_481_520:
        raise ValueError(f"Brazil primary linear model rows should reproduce 1,481,520; got {linear['model_rows']}")
    old = pd.read_csv(OLD_LINEAR).query("outcome == 'N20_N23'").iloc[0]
    if not np.isclose(temp["RR"], float(old.heat_RR_per_C), rtol=0, atol=1e-8):
        raise ValueError(f"Brazil primary temperature RR did not reproduce frozen estimate: {temp['RR']} vs {old.heat_RR_per_C}")
    if not np.isclose(fire["RR"], float(old.smoke_RR_per10), rtol=0, atol=1e-8):
        raise ValueError(f"Brazil primary fire RR did not reproduce frozen estimate: {fire['RR']} vs {old.smoke_RR_per10}")

    nl_terms = add_spline_terms(frame, t_names, knot_info["knots_P5_P33_333_P66_667_P95_C"])
    spline_terms = nl_terms + f_names
    spline = load_or_fit(
        "brazil_primary_spline3df_temp_fire_two_way", "N20_N23", spline_terms,
        frame, "stratum", "hrid", "day_id", int(lag_audit["source_rows_before_lag"]),
        float(lag_audit["source_events_full_window"]), "Brazil_primary", "Brazil N20-N23 primary",
    )
    test_rows.append(nonlinear_wald_test(spline, "Brazil_primary_2008_2022"))
    crows, curves = spline_curve(spline, knot_info, "Brazil_primary_2008_2022", "temperature")
    contrast_rows.extend(crows)
    curve_rows.extend(curves)
    del linear, spline, frame
    gc.collect()
    return temp, fire


def fit_brazil_compound(contrast_rows: list[dict[str, Any]], primary_audit: dict[str, Any]) -> dict[str, Any]:
    con = duckdb.connect()
    con.execute(f"SET threads={THREADS}")
    frame = con.execute(
        "SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date, "
        "CAST(N20_N23 AS DOUBLE) AS N20_N23, CAST(heat AS INTEGER) AS heat, "
        "CAST(smoke AS INTEGER) AS smoke "
        "FROM read_parquet(?) WHERE date >= DATE '2008-01-01' AND date < DATE '2023-01-01'",
        [str(BRAZIL_PANEL)],
    ).df()
    con.close()
    frame["date"] = pd.to_datetime(frame.date, errors="raise").dt.normalize()
    frame["hrid"] = frame.hrid.astype(str)
    frame["day_id"] = frame.date.dt.strftime("%Y%m%d")
    frame["ym"] = frame.date.dt.strftime("%Y%m")
    frame["dow"] = frame.date.dt.dayofweek.astype("int8")
    frame["stratum"] = frame.hrid + "_" + frame.ym + "_" + frame.dow.astype(str)
    frame["heat_only"] = ((frame.heat == 1) & (frame.smoke == 0)).astype("int8")
    frame["smoke_only"] = ((frame.heat == 0) & (frame.smoke == 1)).astype("int8")
    frame["both"] = ((frame.heat == 1) & (frame.smoke == 1)).astype("int8")
    if frame[["N20_N23", "heat", "smoke"]].isna().any().any():
        raise ValueError("Missing Brazil same-day compound inputs")
    frame["hrid"] = frame.hrid.astype("category")
    frame["day_id"] = frame.day_id.astype("category")
    frame["stratum"] = frame.stratum.astype("category")
    record = load_or_fit(
        "brazil_primary_same_day_compound_two_way", "N20_N23",
        ["heat_only", "smoke_only", "both"], frame, "stratum", "hrid", "day_id",
        len(frame), float(frame.N20_N23.sum()), "Brazil_compound", "Brazil N20-N23 same-day categories",
    )
    for term, label in [("heat_only", "heat_only_vs_neither"),
                        ("smoke_only", "smoke_only_vs_neither"),
                        ("both", "compound_vs_neither")]:
        contrast_rows.append(estimate(record, {term: 1.0}, label))
    primary_audit["same_day_compound_source_rows"] = len(frame)
    primary_audit["same_day_compound_events_full_window"] = float(frame.N20_N23.sum())
    primary_audit["same_day_compound_model_rows"] = record["model_rows"]
    del frame
    gc.collect()
    return record


def save_derived_tables(contrast_rows: list[dict[str, Any]], curve_rows: list[dict[str, Any]],
                        test_rows: list[dict[str, Any]], knot_rows: dict[str, Any],
                        source_audit: dict[str, Any], primary_audit: dict[str, Any]) -> None:
    pd.DataFrame(contrast_rows).to_csv(OUT / "temperature_and_category_contrasts.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(OUT / "temperature_curves.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(test_rows).to_csv(OUT / "joint_nonlinearity_tests.csv", index=False, encoding="utf-8-sig")
    save_json(OUT / "knots_and_percentiles.json", knot_rows)
    save_json(OUT / "source_sample_audit.json", {"strict_country_frames": source_audit,
                                                  "brazil_primary": primary_audit})
    export_fit_tables()


def historical_inference(primary_temp: dict[str, Any], primary_fire: dict[str, Any],
                         compound_record: dict[str, Any]) -> pd.DataFrame:
    old_linear = pd.read_csv(OLD_LINEAR).query("outcome == 'N20_N23'").iloc[0]
    updated = {
        "temperature_lag0_7": primary_temp,
        "fire_lag0_7": primary_fire,
    }
    rows: list[dict[str, Any]] = []
    for key, term, old_est, old_low, old_high, old_p in [
        ("temperature_lag0_7", "heat", old_linear.heat_RR_per_C, old_linear.heat_low, old_linear.heat_high, old_linear.heat_p),
        ("fire_lag0_7", "smoke", old_linear.smoke_RR_per10, old_linear.smoke_low, old_linear.smoke_high, old_linear.smoke_p),
    ]:
        current = updated[key]
        rows.append({
            "analysis": key, "contrast": term, "historical_RR": float(old_est),
            "historical_CI_low": float(old_low), "historical_CI_high": float(old_high),
            "historical_p": float(old_p), "historical_inference": "archived one-way health-region cluster; unchanged and not overwritten",
            "updated_RR": current["RR"], "updated_CI_low": current["CI_low"],
            "updated_CI_high": current["CI_high"], "updated_p_t": current["p_t"],
            "updated_se": current["se"], "updated_t_df": current["t_df"],
            "updated_inference": "two-way CRV1 hrid + calendar date with intersection term; retained-cluster t reference",
            "updated_model": current["model"],
        })
    old_compound = pd.read_csv(OLD_COMPOUND).query("outcome == 'N20_N23'")
    old_map = {"heat_only": "heat_only", "smoke_only": "smoke_only", "composite": "both"}
    for _, old in old_compound.iterrows():
        term = old_map[str(old.group)]
        current = estimate(compound_record, {term: 1.0}, f"{old.group}_vs_neither")
        rows.append({
            "analysis": "same_day_compound_categories", "contrast": str(old.group),
            "historical_RR": float(old.RR), "historical_CI_low": float(old.CI_low),
            "historical_CI_high": float(old.CI_high), "historical_p": float(old.p),
            "historical_inference": "archived one-way health-region cluster; unchanged and not overwritten",
            "updated_RR": current["RR"], "updated_CI_low": current["CI_low"],
            "updated_CI_high": current["CI_high"], "updated_p_t": current["p_t"],
            "updated_se": current["se"], "updated_t_df": current["t_df"],
            "updated_inference": "two-way CRV1 hrid + calendar date with intersection term; retained-cluster t reference",
            "updated_model": compound_record["model_key"],
        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "historical_vs_updated_inference.csv", index=False, encoding="utf-8-sig")
    return result


def write_method_report(source_audit: dict[str, Any], primary_audit: dict[str, Any],
                        knot_rows: dict[str, Any], tests: list[dict[str, Any]],
                        historical: pd.DataFrame) -> None:
    def fmt(value: Any, digits: int = 4) -> str:
        return "NA" if value is None or pd.isna(value) else f"{float(value):.{digits}f}"

    def p_fmt(value: Any) -> str:
        if value is None or pd.isna(value):
            return "NA"
        return "<0.0001" if float(value) < 0.0001 else f"{float(value):.4f}"

    tests_by = {row["family"]: row for row in tests}
    contrasts = pd.read_csv(OUT / "temperature_and_category_contrasts.csv")
    samples = pd.read_csv(OUT / "sample_flow.csv")
    diagnostics = pd.read_csv(OUT / "covariance_diagnostics.csv")
    lines = [
        "# M3 温度非线性敏感性与巴西主分析推断更新",
        "",
        f"运行状态：complete。规格文件在模型结果计算前固定，SHA-256：`{SPEC_SHA256}`。所有新增产物限于本目录；原始数据和历史输出未覆盖。",
        "",
        "## 数据与分析样本",
        "",
        f"- Brazil strict N20：2018–2022年，{source_audit['Brazil']['state_clusters_pre_lag']}州、{source_audit['Brazil']['source_rows']:,}州日、{source_audit['Brazil']['events_full_window']:,}例；温度沿用既有巴西健康区人口加权州级赋值。",
        f"- Mexico strict N20：2018–2022年，{source_audit['Mexico']['state_clusters_pre_lag']}州、{source_audit['Mexico']['source_rows']:,}州日、{source_audit['Mexico']['events_full_window']:,}例；病例使用扩展日计数文件中的`admissions`，其中包含来源年2023、入院日期仍在2018–2022的记录；温度使用ERA5州日` tmean`。",
        f"- Brazil primary：健康区日面板筛选为2008-01-01至2022-12-31后才生成lag，{primary_audit['regions']:,}个健康区；完整lag 0–7前有{primary_audit['source_rows_before_lag']:,}行，分析frame有{primary_audit['analysis_rows_after_complete_lags']:,}行。温度每滞后日用°C，火源PM₂.₅按每10 µg/m³缩放。该裁窗/造lag顺序保持原主分析口径。",
        "- Brazil strict及Mexico strict均先限制2018–2022窗口，再在州内生成lag 0–7；线性与样条模型在相同完整lag分析frame上拟合。分层固定效应为州×年月×星期几，原有国家公共假日指示项保留。巴西primary延用健康区×年月×星期几分层。",
        "",
        "## 预设模型与推断",
        "",
        "每个温度暴露使用三自由度自然立方（restricted cubic）样条：原始°C线性列加两个非线性列。每个国家在同一匹配分析frame上，将lag 0–7温度合并后，不看结局计算P5、P33.333、P66.667、P95作为边界/内部结点；边界外观测全部保留，按自然样条线性外推。每个lag复用相同暴露基函数。",
        "- strict N20模型为温度lag 0–7；唯一变化是把8个线性温度项换成每lag的线性+2个自然样条列。样条模型联合检验16个非线性系数为零，条件于线性温度项及不变的节假日项。",
        "- Brazil 2008–2022 primary保持N20–N23结局、温度及火源PM₂.₅ lag 0–7和原健康区×年月×星期几分层；样条敏感性只替换温度滞后基函数，火源滞后保持不变。",
        "- 所有新模型的CRV1协方差按州/健康区与日历日期双向聚类，含交集项。对称化后报告协方差特征值诊断；未将矩阵投影到PSD。对比和系数采用双侧95%区间及t分布，df=min(G1−1,G2−1)。曲线为相对P50的点态95%区间，不是同时置信带。",
        "",
        "## 匹配样本流与聚类支持",
        "",
        "详细逐模型行数、事件数、零事件层、拟合保留行数和保留簇数见`sample_flow.csv`。全部系数见`coefficients.csv`；各模型完整协方差矩阵见`covariance_matrices.json`。",
        "",
        "| 模型 | 拟合行数 | 拟合事件 | G1 | G2 | 交集簇 | t df |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in samples.iterrows():
        lines.append(f"| {row['model']} | {int(row['model_rows_retained']):,} | {int(row['events_model_retained']):,} | {int(row['G_cluster_1_retained']):,} | {int(row['G_cluster_2_retained']):,} | {int(row['G_intersection_retained']):,} | {int(row['t_df'])} |")
    lines.extend(["", "## 温度函数形式结果", "", "| 分析 | 联合非线性F检验 | 分子df | 分母df | p值 |", "|---|---:|---:|---:|---:|"])
    for family in ["Brazil", "Mexico", "Brazil_primary_2008_2022"]:
        row = tests_by.get(family, {})
        lines.append(f"| {family} | {fmt(row.get('F'))} | {fmt(row.get('df_num'), 0)} | {fmt(row.get('df_denom'), 0)} | {p_fmt(row.get('p_F'))} |")
    invalid_tests = [row for row in tests if row.get("status") != "ok"]
    if invalid_tests:
        lines.append("")
        for row in invalid_tests:
            lines.append(
                f"**推断限制：** {row.get('family')}的16维非线性限制协方差含负特征值"
                f"（最小值 {fmt(row.get('restriction_covariance_min_eigenvalue'), 7)}）；"
                "不对协方差作PSD投影，因此保留Wald F统计量但不报告其p值。"
            )
    lines.extend(["", "温度百分位对比与100点曲线逐点结果见`temperature_and_category_contrasts.csv`及`temperature_curves.csv`。曲线覆盖各自分析frame pooled lag-0–7温度P5–P95；百分位对比为P10/P50/P90/P95相对P50，均对lag 0–7基函数系数作联合对比。", "", "## 巴西primary历史估计与双向推断"])
    for _, row in historical.iterrows():
        lines.append(f"- {row['analysis']} / {row['contrast']}：旧估计RR {fmt(row['historical_RR'])}（{fmt(row['historical_CI_low'])}–{fmt(row['historical_CI_high'])}，p={p_fmt(row['historical_p'])}）；双向推断RR {fmt(row['updated_RR'])}（{fmt(row['updated_CI_low'])}–{fmt(row['updated_CI_high'])}，p={p_fmt(row['updated_p_t'])}，df={int(row['updated_t_df'])}）。旧列只作历史参照，独立保存于`historical_vs_updated_inference.csv`。")
    lines.extend(["", "## 协方差诊断与解释边界", "", "| 模型 | 最小特征值 | 最大特征值 | 负特征值数 | 相对负值尺度 |", "|---|---:|---:|---:|---:|"])
    for _, row in diagnostics.iterrows():
        lines.append(f"| {row['model']} | {float(row['minimum_eigenvalue']):.6g} | {float(row['maximum_eigenvalue']):.6g} | {int(row['negative_eigenvalue_count'])} | {float(row['relative_negative_eigenvalue_scale']):.6g} |")
    psd_diagnostic_path = OUT / "psd_covariance_diagnostics.csv"
    psd_test_path = OUT / "joint_nonlinearity_tests_psd.csv"
    if psd_diagnostic_path.exists() and psd_test_path.exists():
        psd_diagnostics = pd.read_csv(psd_diagnostic_path)
        psd_tests = pd.read_csv(psd_test_path)
        psd_contrasts = pd.read_csv(OUT / "temperature_contrasts_psd.csv")
        p90 = psd_contrasts.loc[psd_contrasts.percentile.eq("P90")]
        lines.extend([
            "",
            "## 显式全矩阵PSD敏感性",
            "",
            "额外进行一项明确标记的协方差敏感性：对三个样条模型的完整对称系数协方差矩阵作特征分解，将负特征值设为零后重构。原始矩阵、原始区间和原始检验均保留且未覆盖；PSD结果是补充推断敏感性，不替代主输出。该处理也未按结果方向或显著性选择。",
            "",
            "| 样条模型 | 原始相对负特征值尺度 | PSD Frobenius改变量/原矩阵 | PSD最大对角方差改变量 |",
            "|---|---:|---:|---:|",
        ])
        for _, row in psd_diagnostics.iterrows():
            lines.append(f"| {row['model']} | {float(row['raw_relative_negative_eigenvalue_scale']):.6g} | {float(row['relative_frobenius_correction']):.6g} | {float(row['max_absolute_diagonal_variance_change']):.6g} |")
        lines.extend([
            "",
            "| 样条模型 | PSD敏感性Wald F | 分子df | 分母df | PSD p值 | 原始协方差p值 |",
            "|---|---:|---:|---:|---:|---:|",
        ])
        for _, row in psd_tests.iterrows():
            lines.append(f"| {row['family']} | {fmt(row.get('F'))} | {fmt(row.get('df_num'), 0)} | {fmt(row.get('df_denom'), 0)} | {p_fmt(row.get('p_F'))} | {p_fmt(row.get('raw_p_F'))} |")
        lines.extend(["", "| 分析 | P90 vs P50 RR（95% CI） | PSD p值 |", "|---|---:|---:|"])
        for _, row in p90.iterrows():
            lines.append(f"| {row['country_or_family']} | {fmt(row['RR'])}（{fmt(row['CI_low'])}–{fmt(row['CI_high'])}） | {p_fmt(row['p_t'])} |")
        lines.extend([
            "",
            "三模型完整PSD协方差见`covariance_psd_sensitivity.json`；PSD曲线、百分位对比、联合非线性检验分别见`temperature_curves_psd.csv`、`temperature_contrasts_psd.csv`、`joint_nonlinearity_tests_psd.csv`。raw-vs-PSD区间变动见`raw_vs_psd_curve_uncertainty.csv`与`raw_vs_psd_percentile_uncertainty.csv`。Mexico原始16维限制协方差未给出可用的Wald p值；上表另列的是对完整协方差作显式特征值归零后的敏感性检验，不把原始p值缺失解读为无非线性。",
        ])
    lines.extend([
        "",
        "本分析只检验预设的有限函数形式并更新聚类推断，不按显著性搜索规格。严格N20与Brazil primary的结局范围和时间窗不同；结果描述的是行政住院利用关联，不能推断结石发生时间或因果效应。样条曲线的逐点区间不控制整条曲线的同时覆盖率。",
        "",
        f"## 复现与环境\n\n脚本：`work/revision_20261113/nonlinear/run_m3_nonlinear.py`。运行Python：3.12；分析线程上限：{THREADS}。输入路径、大小及SHA-256见`input_manifest.json`。状态、进程号、当前阶段和日志位置见`run_status.json`。",
    ])
    (OUT / "方法报告_中文.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    write_status("running", "input_hashes")
    input_manifest()
    log("INPUTS_HASHED", files=len(SOURCE_HASHES), spec_sha256=SPEC_SHA256)
    holidays = load_holidays()
    write_status("running", "prepare_strict_frames")
    brazil_strict, mexico_strict, strict_audit = prepare_strict_frames(holidays)
    log("STRICT_FRAMES_READY", brazil=strict_audit["Brazil"], mexico=strict_audit["Mexico"])

    contrast_rows: list[dict[str, Any]] = []
    curve_rows: list[dict[str, Any]] = []
    test_rows: list[dict[str, Any]] = []
    knot_rows: dict[str, Any] = {}
    primary_audit: dict[str, Any] = {}

    fit_strict_country("Brazil", brazil_strict, strict_audit["Brazil"],
                       curve_rows, contrast_rows, test_rows, knot_rows)
    save_derived_tables(contrast_rows, curve_rows, test_rows, knot_rows, strict_audit, primary_audit)
    fit_strict_country("Mexico", mexico_strict, strict_audit["Mexico"],
                       curve_rows, contrast_rows, test_rows, knot_rows)
    save_derived_tables(contrast_rows, curve_rows, test_rows, knot_rows, strict_audit, primary_audit)
    del brazil_strict, mexico_strict
    gc.collect()

    write_status("running", "Brazil_primary_temperature_fire")
    primary_temp, primary_fire = fit_brazil_primary(
        curve_rows, contrast_rows, test_rows, knot_rows, primary_audit,
    )
    save_derived_tables(contrast_rows, curve_rows, test_rows, knot_rows, strict_audit, primary_audit)
    write_status("running", "Brazil_primary_same_day_compound")
    compound = fit_brazil_compound(contrast_rows, primary_audit)
    historical = historical_inference(primary_temp, primary_fire, compound)
    save_derived_tables(contrast_rows, curve_rows, test_rows, knot_rows, strict_audit, primary_audit)
    write_method_report(strict_audit, primary_audit, knot_rows, test_rows, historical)

    # Fixed historical point-estimate anchors ensure that the original lag window was recreated.
    old_compound = pd.read_csv(OLD_COMPOUND).query("outcome == 'N20_N23'")
    new_compound_map = {"heat_only": "heat_only", "smoke_only": "smoke_only", "composite": "both"}
    for _, old in old_compound.iterrows():
        term = new_compound_map[str(old.group)]
        current = estimate(compound, {term: 1.0}, str(old.group))
        if not np.isclose(current["RR"], float(old.RR), rtol=0, atol=1e-8):
            raise ValueError(f"Historical compound point estimate failed reproduction for {old.group}: {current['RR']} vs {old.RR}")

    metadata = {
        "status": "complete",
        "code_revision": CODE_REVISION,
        "specification_sha256": SPEC_SHA256,
        "started_utc": datetime.fromtimestamp(T0, timezone.utc).isoformat(),
        "completed_utc": utc_now(),
        "elapsed_seconds": round(time.time() - T0, 2),
        "threads": THREADS,
        "input_hashes": SOURCE_HASHES,
        "models_completed": list(FIT_RECORDS),
        "outputs": sorted(path.name for path in OUT.iterdir() if path.is_file()),
        "historical_point_estimates_reproduced": True,
    }
    save_json(OUT / "analysis_metadata.json", metadata)
    write_status("complete", "complete", completed_models=list(FIT_RECORDS),
                 elapsed_seconds=metadata["elapsed_seconds"], outputs=metadata["outputs"])
    log("ANALYSIS_COMPLETE", models=len(FIT_RECORDS), elapsed_seconds=metadata["elapsed_seconds"],
        outputs=metadata["outputs"])


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        failure = {
            "status": "failed", "phase": "exception", "pid": os.getpid(),
            "updated_utc": utc_now(), "elapsed_seconds": round(time.time() - T0, 2),
            "error": repr(exc), "traceback": traceback.format_exc(),
            "log": str(LOG_PATH),
        }
        save_json(STATUS_PATH, failure)
        log("ANALYSIS_FAILED", error=repr(exc), traceback=traceback.format_exc())
        raise
