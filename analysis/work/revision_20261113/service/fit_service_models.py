"""Pre-specified revision Brazil service-model reanalysis.

Reads established aggregate panels and writes only to the revision output
directory. Each new model is fitted once; contrasts reuse its full covariance.
"""
from __future__ import annotations

import os

for _name in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMBA_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_name] = "4"

import hashlib
import json
import math
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

sys.path.insert(0, "F:/analysis/climate_biliary_20261102/work/pylibs")

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[3]
WORK_DIR = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "revision_20261113" / "service"
DATA = Path("F:/analysis/climate_biliary_20261102/work")
EXPOSURE_PATH = DATA / "panel_composite_humidity.parquet"
CLINICAL_PATH = ROOT / "work" / "upgrade_20261112" / "clinical_code_day_aggregate.parquet"
PRICE_PATH = ROOT / "outputs" / "upgrade_20261112" / "ipca_monthly.csv"
SP_LINK_PATH = ROOT / "outputs" / "upgrade_20261112" / "sp_linked_daily_outcomes.csv"
HOLIDAY_PATH = ROOT / "outputs" / "upgrade_20261112" / "holiday_calendar_dates.csv"
LEGACY_RESULTS_PATH = ROOT / "outputs" / "upgrade_20261112" / "acute_repair_model_results.csv"
STATUS_PATH = OUT / "attempt_02_status.json"
CLUSTERS = ["hrid", "day_id"]
GROUPS = ["neither", "heat_only", "smoke_only", "both"]
CELL_KEYS = [f"{group}_h{holiday}" for group in GROUPS for holiday in (0, 1)]
SERVICE_OUTCOMES = [
    "urgent_N20",
    "admission_attributed_beddays",
    "reimbursement_2022BRL",
]
EXPECTED_COUNTS = {
    "broad_urgent": 315_504,
    "broad_elective": 115_870,
    "strict_urgent": 226_402,
    "strict_elective": 94_799,
}
ABSOLUTE_DRAWS = 4096
ABSOLUTE_SEED = 20261113
RESULT_ROWS: list[dict[str, Any]] = []
SUMMARY_ROWS: list[dict[str, Any]] = []
SAMPLE_AUDIT_ROWS: list[dict[str, Any]] = []
EVENT_ROWS: list[dict[str, Any]] = []
SAMPLEFLOW_ROWS: list[dict[str, Any]] = []
FIT_RECORDS: dict[str, dict[str, Any]] = {}
ABSOLUTE_ROWS: list[dict[str, Any]] = []
STARTED = datetime.now(timezone.utc).isoformat()
OUT.mkdir(parents=True, exist_ok=True)


def log(*parts: Any) -> None:
    print(*parts, flush=True)


def json_scalar(value: Any) -> Any:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if pd.isna(value):
        return None
    raise TypeError(f"Cannot JSON encode {type(value).__name__}")


def save_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_scalar),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def persist_tables() -> None:
    pd.DataFrame(RESULT_ROWS).to_csv(OUT / "results.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(SUMMARY_ROWS).to_csv(OUT / "model_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(SAMPLE_AUDIT_ROWS).to_csv(OUT / "sample_audit.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(EVENT_ROWS).to_csv(OUT / "model_events.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(SAMPLEFLOW_ROWS).to_csv(OUT / "sampleflow.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(ABSOLUTE_ROWS).to_csv(
        OUT / "absolute_service_contrasts.csv", index=False, encoding="utf-8-sig"
    )
    save_json(OUT / "model_fits.json", FIT_RECORDS)


def add_flow(
    stage: str,
    frame: pd.DataFrame,
    note: str,
    event_columns: list[str] | None = None,
) -> None:
    row: dict[str, Any] = {
        "stage": stage,
        "rows": int(len(frame)),
        "health_regions": int(frame.hrid.nunique()) if "hrid" in frame else None,
        "dates": int(frame.date.nunique()) if "date" in frame else None,
        "min_date": frame.date.min().date().isoformat() if "date" in frame and len(frame) else None,
        "max_date": frame.date.max().date().isoformat() if "date" in frame and len(frame) else None,
        "note": note,
    }
    for column in event_columns or []:
        row[f"events_{column}"] = float(frame[column].sum())
    SAMPLEFLOW_ROWS.append(row)
    pd.DataFrame(SAMPLEFLOW_ROWS).to_csv(OUT / "sampleflow.csv", index=False, encoding="utf-8-sig")


def read_parquet(query: str, path: Path) -> pd.DataFrame:
    with duckdb.connect() as connection:
        connection.execute("SET threads=4")
        return connection.execute(query, [str(path)]).df()


def load_holiday_dates() -> tuple[set[pd.Timestamp], pd.DataFrame, str]:
    if not HOLIDAY_PATH.exists():
        raise FileNotFoundError(f"Existing national holiday calendar missing: {HOLIDAY_PATH}")
    calendar = pd.read_csv(HOLIDAY_PATH, dtype={"country": "string"})
    required = {"country", "date"}
    if not required.issubset(calendar.columns):
        raise ValueError(f"Holiday calendar columns missing: {sorted(required - set(calendar.columns))}")
    brazil = calendar.loc[calendar.country.eq("BR")].copy()
    brazil["date"] = pd.to_datetime(brazil.date, errors="raise").dt.normalize()
    if brazil.empty or brazil.date.duplicated().any():
        raise ValueError("Existing Brazil national holiday calendar is empty or has duplicate dates")
    digest = hashlib.sha256(HOLIDAY_PATH.read_bytes()).hexdigest()
    keep = [
        column
        for column in ["country", "date", "holiday_name", "date_role", "calendar_source"]
        if column in brazil
    ]
    brazil[keep].to_csv(OUT / "holiday_calendar_used.csv", index=False, encoding="utf-8-sig")
    return set(brazil.date), brazil, digest


def prepare_base_panel() -> tuple[pd.DataFrame, dict[str, Any]]:
    query = """
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               tmean, fire_pm25, heat, smoke
        FROM read_parquet(?)
        WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
    """
    panel = read_parquet(query, EXPOSURE_PATH)
    panel["hrid"] = panel.hrid.astype(str)
    panel["date"] = pd.to_datetime(panel.date, errors="raise").dt.normalize()
    panel = panel.sort_values(["hrid", "date"]).reset_index(drop=True)
    if panel.duplicated(["hrid", "date"]).any():
        raise ValueError("Exposure panel has duplicate health-region dates")
    if panel[["tmean", "fire_pm25", "heat", "smoke"]].isna().any().any():
        raise ValueError("Exposure or category flags are missing in the 2018-2022 panel")
    if len(panel) != 450 * 1826 or panel.hrid.nunique() != 450 or panel.date.nunique() != 1826:
        raise ValueError(
            f"Expected complete 450-region x 1,826-day panel, got rows={len(panel)}, "
            f"regions={panel.hrid.nunique()}, dates={panel.date.nunique()}"
        )
    if not panel.groupby("hrid", sort=False).date.diff().dropna().dt.days.eq(1).all():
        raise ValueError("Exposure panel has within-region calendar gaps")
    if not set(panel.heat.unique()).issubset({0, 1}) or not set(panel.smoke.unique()).issubset({0, 1}):
        raise ValueError("Heat and smoke indicators are not binary")
    panel["ym"] = panel.date.dt.strftime("%Y%m")
    panel["dow"] = panel.date.dt.dayofweek.astype(int)
    panel["day_id"] = panel.date.dt.strftime("%Y%m%d")
    panel["stratum"] = panel.hrid + "_" + panel.ym + "_" + panel.dow.astype(str)
    panel["h_only"] = ((panel.heat.eq(1)) & (panel.smoke.eq(0))).astype("int8")
    panel["s_only"] = ((panel.heat.eq(0)) & (panel.smoke.eq(1))).astype("int8")
    panel["both"] = ((panel.heat.eq(1)) & (panel.smoke.eq(1))).astype("int8")
    panel["group"] = np.select(
        [panel.h_only.eq(1), panel.s_only.eq(1), panel.both.eq(1)],
        ["heat_only", "smoke_only", "both"],
        default="neither",
    )
    panel["group"] = pd.Categorical(panel.group, categories=GROUPS)
    holiday_dates, holiday_frame, holiday_sha = load_holiday_dates()
    panel["holiday"] = panel.date.isin(holiday_dates).astype("int8")
    add_flow(
        "exposure_panel_2018_2022",
        panel,
        "Complete established Brazil health-region by calendar-day exposure panel.",
    )
    audit = {
        "exposure_panel_path": str(EXPOSURE_PATH),
        "exposure_panel_sha256": hashlib.sha256(EXPOSURE_PATH.read_bytes()).hexdigest(),
        "rows": int(len(panel)),
        "regions": int(panel.hrid.nunique()),
        "dates": int(panel.date.nunique()),
        "holiday_calendar_path": str(HOLIDAY_PATH),
        "holiday_calendar_sha256": holiday_sha,
        "brazil_holiday_dates_in_panel": int(panel.holiday.sum()),
        "brazil_calendar_dates": int(len(holiday_frame)),
        "fixed_effect": "health region x year-month x day-of-week",
    }
    return panel, audit


def prepare_service_panel(panel: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    query = """
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               CAST(competence AS VARCHAR) AS competence,
               admissions, beddays, nominal_cost
        FROM read_parquet(?)
        WHERE IDENT = '1' AND diagnosis = 'N20' AND CAR_INT = '02'
          AND date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
    """
    clinical = read_parquet(query, CLINICAL_PATH)
    clinical["date"] = pd.to_datetime(clinical.date, errors="raise").dt.normalize()
    prices = pd.read_csv(PRICE_PATH, dtype={"month": "string"})
    required = {"month", "factor_to_2022_12"}
    if not required.issubset(prices.columns):
        raise ValueError(f"IPCA table columns missing: {sorted(required - set(prices.columns))}")
    factors = dict(zip(prices.month.astype(str), pd.to_numeric(prices.factor_to_2022_12, errors="raise")))
    clinical["competence"] = clinical.competence.astype(str).str.strip()
    clinical["price_factor"] = clinical.competence.map(factors)
    if clinical.price_factor.isna().any():
        missing = clinical.loc[clinical.price_factor.isna(), "competence"].drop_duplicates().tolist()
        raise ValueError(f"Price index is missing competence months: {missing[:12]}")
    clinical["reimbursement_2022BRL"] = (
        pd.to_numeric(clinical.nominal_cost, errors="raise")
        * pd.to_numeric(clinical.price_factor, errors="raise")
    )
    clinical["admissions"] = pd.to_numeric(clinical.admissions, errors="raise")
    clinical["beddays"] = pd.to_numeric(clinical.beddays, errors="raise")
    daily = (
        clinical.groupby(["hrid", "date"], as_index=False, observed=True)
        .agg(
            urgent_N20=("admissions", "sum"),
            admission_attributed_beddays=("beddays", "sum"),
            reimbursement_2022BRL=("reimbursement_2022BRL", "sum"),
        )
    )
    if daily.duplicated(["hrid", "date"]).any():
        raise ValueError("Strict N20 daily outcome aggregate has duplicate region dates")
    if int(daily.urgent_N20.sum()) != 226_402:
        raise ValueError(f"Strict N20 urgent admissions differ from established total: {daily.urgent_N20.sum()}")
    merged = panel.merge(daily, on=["hrid", "date"], how="left", validate="one_to_one")
    for outcome in SERVICE_OUTCOMES:
        merged[outcome] = pd.to_numeric(merged[outcome], errors="coerce").fillna(0.0)
    if len(merged) != len(panel) or merged[SERVICE_OUTCOMES].isna().any().any():
        raise ValueError("Strict N20 outcomes did not restore to the complete region-day panel")
    add_flow(
        "strict_N20_urgent_outcomes_joined",
        merged,
        "IDENT=1, diagnosis=N20, CAR_INT=02; sparse outcome dates zero-filled after validated one-to-one join.",
        SERVICE_OUTCOMES,
    )
    audit = {
        "clinical_aggregate_path": str(CLINICAL_PATH),
        "clinical_aggregate_sha256": hashlib.sha256(CLINICAL_PATH.read_bytes()).hexdigest(),
        "clinical_source_rows": int(len(clinical)),
        "clinical_region_days": int(len(daily)),
        "clinical_urgent_admissions": int(daily.urgent_N20.sum()),
        "clinical_admission_attributed_beddays": float(daily.admission_attributed_beddays.sum()),
        "clinical_reimbursement_2022BRL": float(daily.reimbursement_2022BRL.sum()),
        "price_index_path": str(PRICE_PATH),
        "price_index_sha256": hashlib.sha256(PRICE_PATH.read_bytes()).hexdigest(),
        "price_reference_month": "2022-12",
        "input_panel_rows_after_zero_fill": int(len(merged)),
        "zero_filled_region_days": int(merged.urgent_N20.eq(0).sum()),
    }
    return merged, audit


def prepare_drainage_panel(panel: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    query = """
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS VARCHAR) AS date,
               any_drainage
        FROM read_csv_auto(?, header=true, all_varchar=true)
    """
    linked = read_parquet(query, SP_LINK_PATH)
    linked["hrid"] = linked.hrid.astype(str).str.strip()
    linked["date"] = pd.to_datetime(
        linked.date.astype(str).str.strip().str.replace(r"\.0$", "", regex=True),
        format="%Y%m%d",
        errors="coerce",
    )
    if linked[["hrid", "date", "any_drainage"]].isna().any().any():
        raise ValueError("SP link has missing or invalid keys/outcomes")
    if linked.duplicated(["hrid", "date"]).any():
        raise ValueError("SP link has duplicate region dates")
    linked["any_drainage"] = pd.to_numeric(linked.any_drainage, errors="raise")
    if (linked.any_drainage < 0).any() or not np.allclose(linked.any_drainage, np.round(linked.any_drainage)):
        raise ValueError("SP any_drainage must be a nonnegative integer count")
    if int(linked.any_drainage.sum()) != 84_764:
        raise ValueError(f"SP any_drainage total changed from established link: {linked.any_drainage.sum()}")
    linked["date"] = linked.date.dt.normalize()
    outside = linked.merge(panel[["hrid", "date"]], on=["hrid", "date"], how="left", indicator=True)
    if outside._merge.ne("both").any():
        raise ValueError("SP linked outcomes include region dates outside the exposure panel")
    merged = panel.merge(
        linked[["hrid", "date", "any_drainage"]],
        on=["hrid", "date"],
        how="left",
        validate="one_to_one",
    )
    zero_filled = int(merged.any_drainage.isna().sum())
    merged["any_drainage"] = pd.to_numeric(merged.any_drainage, errors="coerce").fillna(0.0)
    if len(merged) != 821_700 or int(merged.any_drainage.sum()) != 84_764:
        raise ValueError("SP daily link did not restore to the validated full panel")
    for lag in range(8):
        merged[f"t{lag}"] = merged.groupby("hrid", sort=False).tmean.shift(lag)
        merged[f"f{lag}"] = merged.groupby("hrid", sort=False).fire_pm25.shift(lag) / 10.0
    lag_terms = [f"t{lag}" for lag in range(8)] + [f"f{lag}" for lag in range(8)]
    lagged = merged.dropna(subset=lag_terms).copy()
    if len(lagged) != 818_550:
        raise ValueError(f"Expected 818,550 drainage lag-complete rows; got {len(lagged)}")
    add_flow(
        "SP_any_drainage_linked_zero_filled",
        merged,
        f"Existing SP linked outcome; restored {zero_filled} unlinked region-days as zero.",
        ["any_drainage"],
    )
    add_flow(
        "SP_any_drainage_lag0_7_complete",
        lagged,
        "First seven days within each region removed only to form complete chronological lag 0-7 predictors.",
        ["any_drainage"],
    )
    audit = {
        "SP_link_path": str(SP_LINK_PATH),
        "SP_link_sha256": hashlib.sha256(SP_LINK_PATH.read_bytes()).hexdigest(),
        "linked_positive_region_days": int(len(linked)),
        "zero_filled_region_days": zero_filled,
        "full_panel_rows": int(len(merged)),
        "full_panel_regions": int(merged.hrid.nunique()),
        "full_panel_days": int(merged.date.nunique()),
        "any_drainage_events": int(merged.any_drainage.sum()),
        "lag_complete_rows": int(len(lagged)),
        "lag_complete_events": int(lagged.any_drainage.sum()),
    }
    return lagged, audit


def prepare_urgent_elective_panel(panel: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    query = """
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
          SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') AND CAR_INT='02'
                   THEN admissions ELSE 0 END) AS broad_urgent,
          SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') AND CAR_INT='01'
                   THEN admissions ELSE 0 END) AS broad_elective,
          SUM(CASE WHEN diagnosis='N20' AND CAR_INT='02'
                   THEN admissions ELSE 0 END) AS strict_urgent,
          SUM(CASE WHEN diagnosis='N20' AND CAR_INT='01'
                   THEN admissions ELSE 0 END) AS strict_elective
        FROM read_parquet(?)
        WHERE IDENT='1' AND date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
        GROUP BY 1,2
    """
    daily = read_parquet(query, CLINICAL_PATH)
    daily["hrid"] = daily.hrid.astype(str)
    daily["date"] = pd.to_datetime(daily.date, errors="raise").dt.normalize()
    if daily.duplicated(["hrid", "date"]).any():
        raise ValueError("Urgent/elective aggregate has duplicate region dates")
    frame = panel.merge(daily, on=["hrid", "date"], how="left", validate="one_to_one")
    for outcome, expected in EXPECTED_COUNTS.items():
        frame[outcome] = pd.to_numeric(frame[outcome], errors="coerce").fillna(0).astype("int64")
        total = int(frame[outcome].sum())
        if total != expected:
            raise ValueError(f"{outcome} total differs from established source count: {total} != {expected}")
    if len(frame) != len(panel):
        raise ValueError("Urgent/elective daily merge changed the region-day sample")
    add_flow(
        "urgent_elective_categories_joined",
        frame,
        "Broad N20-N23 and strict N20 urgent/elective totals restored to full region-day panel.",
        list(EXPECTED_COUNTS),
    )
    audit = {
        "daily_aggregate_region_days": int(len(daily)),
        "panel_rows_after_join": int(len(frame)),
        "input_events": {name: int(frame[name].sum()) for name in EXPECTED_COUNTS},
    }
    return frame, audit


def cluster_counts(model_data: pd.DataFrame) -> tuple[dict[str, int], int]:
    absent = [name for name in CLUSTERS if name not in model_data]
    if absent:
        raise ValueError(f"Fitted sample does not retain cluster columns: {absent}")
    counts = {name: int(model_data[name].nunique()) for name in CLUSTERS}
    df = min(counts.values()) - 1
    if df <= 0:
        raise ValueError(f"Invalid small-cluster degrees of freedom: counts={counts}, df={df}")
    return counts, df


def covariance_diagnostics(covariance: np.ndarray) -> dict[str, Any]:
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError(f"Covariance must be square, got {covariance.shape}")
    if not np.isfinite(covariance).all():
        raise ValueError("Covariance has non-finite entries")
    asymmetry = float(np.max(np.abs(covariance - covariance.T))) if covariance.size else 0.0
    if asymmetry > 1e-8:
        raise ValueError(f"Fitted covariance is materially asymmetric: {asymmetry}")
    symmetric = (covariance + covariance.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(symmetric)
    maximum = float(np.max(np.abs(eigenvalues))) if eigenvalues.size else 0.0
    minimum = float(np.min(eigenvalues)) if eigenvalues.size else 0.0
    tolerance = max(1e-12, maximum * 1e-10)
    return {
        "dimension": int(covariance.shape[0]),
        "max_absolute_asymmetry": asymmetry,
        "min_eigenvalue": minimum,
        "max_eigenvalue": float(np.max(eigenvalues)) if eigenvalues.size else 0.0,
        "max_absolute_eigenvalue": maximum,
        "negative_eigenvalue_count": int(np.sum(eigenvalues < 0)),
        "negative_eigenvalue_count_below_tolerance": int(np.sum(eigenvalues < -tolerance)),
        "relative_min_eigenvalue": minimum / maximum if maximum else None,
        "psd_with_1e-10_relative_tolerance": bool(minimum >= -tolerance),
        "eigenvalue_tolerance": tolerance,
        "eigenvalues": eigenvalues.tolist(),
    }


def base_result_row(
    model_id: str,
    model_family: str,
    estimand: str,
    names: list[str],
    beta: np.ndarray,
    covariance: np.ndarray,
    weights: dict[str, float],
    df: int,
    summary: dict[str, Any],
) -> dict[str, Any]:
    vector = np.zeros(len(names), dtype=float)
    for name, value in weights.items():
        if name not in names:
            raise ValueError(f"Required contrast term {name!r} absent in {model_id}: {names}")
        vector[names.index(name)] = value
    estimate = float(vector @ beta)
    variance_raw = float(vector @ covariance @ vector)
    tolerance = max(1e-12, float(np.max(np.abs(covariance))) * 1e-10)
    if variance_raw < -tolerance:
        raise ValueError(f"Negative variance for {model_id}/{estimand}: {variance_raw}")
    variance = max(variance_raw, 0.0)
    se = math.sqrt(variance)
    if se == 0:
        p_value = 0.0 if estimate != 0 else 1.0
        lower_beta = upper_beta = estimate
        statistic = math.copysign(math.inf, estimate) if estimate else 0.0
    else:
        critical = float(student_t.ppf(0.975, df))
        statistic = estimate / se
        p_value = float(2 * student_t.sf(abs(statistic), df))
        lower_beta = estimate - critical * se
        upper_beta = estimate + critical * se
    return {
        "model_id": model_id,
        "model_family": model_family,
        "estimand": estimand,
        "beta": estimate,
        "se": se,
        "RR_or_ratio": float(np.exp(estimate)),
        "CI_low": float(np.exp(lower_beta)),
        "CI_high": float(np.exp(upper_beta)),
        "statistic_t": statistic,
        "p": p_value,
        "inference": "two-way CRV1; Student t",
        "df": int(df),
        "contrast_terms": json.dumps(list(weights), ensure_ascii=False),
        "input_rows": summary["input_rows"],
        "model_rows": summary["model_rows"],
        "input_events": summary["input_events"],
        "model_events": summary["model_events"],
        "retained_hrid_clusters": summary["retained_cluster_counts"]["hrid"],
        "retained_date_clusters": summary["retained_cluster_counts"]["day_id"],
        "formula": summary["formula"],
        "cluster": summary["cluster"],
        "source_status": "new_two_way_fit",
        "note": "",
    }


def add_estimate_row(
    model_id: str,
    model_family: str,
    estimand: str,
    names: list[str],
    beta: np.ndarray,
    covariance: np.ndarray,
    weights: dict[str, float],
    df: int,
    summary: dict[str, Any],
) -> dict[str, Any]:
    row = base_result_row(
        model_id, model_family, estimand, names, beta, covariance, weights, df, summary
    )
    RESULT_ROWS.append(row)
    return row


def summarize_model(
    model_id: str,
    family: str,
    formula: str,
    frame: pd.DataFrame,
    outcome: str,
    model: Any,
    model_data: pd.DataFrame,
    kept_clusters: dict[str, int],
    df: int,
    diagnostics: dict[str, Any],
    elapsed: float,
) -> dict[str, Any]:
    return {
        "model_id": model_id,
        "model_family": family,
        "formula": formula,
        "fixed_effect": "health region x year-month x day-of-week",
        "cluster": "CRV1, two-way by hrid + day_id",
        "inference": "Student t with df=min(retained hrid clusters, retained date clusters)-1",
        "df": int(df),
        "input_rows": int(len(frame)),
        "model_rows": int(model._N),
        "rows_dropped_by_fit": int(len(frame) - int(model._N)),
        "input_events": float(frame[outcome].sum()),
        "model_events": float(model_data[outcome].sum()),
        "input_cluster_counts": {name: int(frame[name].nunique()) for name in CLUSTERS},
        "retained_cluster_counts": kept_clusters,
        "strata_input": int(frame.stratum.nunique()) if "stratum" in frame else None,
        "strata_model": int(model_data.stratum.nunique()) if "stratum" in model_data else None,
        "covariance_eigen_diagnostics": diagnostics,
        "elapsed_seconds": round(elapsed, 3),
    }


def fit_model(
    model_id: str,
    family: str,
    formula: str,
    frame: pd.DataFrame,
    outcome: str,
    expected_terms: list[str],
) -> tuple[Any, list[str], np.ndarray, np.ndarray, pd.DataFrame, dict[str, Any]]:
    started = time.time()
    log("START", model_id, "rows=", len(frame), "events=", float(frame[outcome].sum()), "formula=", formula)
    missing_cluster_columns = [name for name in CLUSTERS if name not in frame.columns]
    if missing_cluster_columns:
        raise ValueError(
            f"{model_id}: two-way cluster columns missing from model input: {missing_cluster_columns}"
        )
    model = pf.fepois(formula, data=frame, vcov={"CRV1": "hrid + day_id"})
    coefficient_series = model.coef()
    names = [str(name) for name in coefficient_series.index]
    beta = coefficient_series.to_numpy(dtype=float)
    covariance_raw = np.asarray(model._vcov, dtype=float)
    if covariance_raw.shape != (len(names), len(names)):
        raise ValueError(
            f"{model_id}: covariance shape {covariance_raw.shape} does not match {len(names)} coefficients"
        )
    missing = sorted(set(expected_terms) - set(names))
    if missing:
        raise ValueError(f"{model_id}: expected coefficients omitted by estimator: {missing}")
    model_data = getattr(model, "_data", None)
    if model_data is None or len(model_data) != int(model._N):
        raise ValueError(f"{model_id}: fitted sample rows are not available for complete audit")
    model_data = model_data.copy()
    kept_clusters, df = cluster_counts(model_data)
    diagnostics = covariance_diagnostics(covariance_raw)
    covariance = (covariance_raw + covariance_raw.T) / 2.0
    summary = summarize_model(
        model_id, family, formula, frame, outcome, model, model_data,
        kept_clusters, df, diagnostics, time.time() - started,
    )
    record = {
        **summary,
        "coefficient_names": names,
        "coefficients": beta.tolist(),
        "covariance_raw": covariance_raw.tolist(),
        "covariance_used_for_linear_contrasts": covariance.tolist(),
        "covariance_symmetrized_for_rounding_only": True,
        "eigen_diagnostics": diagnostics,
    }
    FIT_RECORDS[model_id] = record
    SUMMARY_ROWS.append(summary)
    SAMPLE_AUDIT_ROWS.append(
        {
            **{key: value for key, value in summary.items() if key != "covariance_eigen_diagnostics"},
            "covariance_min_eigenvalue": diagnostics["min_eigenvalue"],
            "covariance_relative_min_eigenvalue": diagnostics["relative_min_eigenvalue"],
            "covariance_negative_eigenvalue_count": diagnostics["negative_eigenvalue_count"],
            "covariance_negative_eigenvalue_count_below_tolerance": diagnostics[
                "negative_eigenvalue_count_below_tolerance"
            ],
        }
    )
    if outcome == "count" and "outcome" in model_data:
        input_counts = frame.groupby("outcome", observed=True)[outcome].sum()
        model_counts = model_data.groupby("outcome", observed=True)[outcome].sum()
        input_rows = frame.groupby("outcome", observed=True).size()
        model_rows = model_data.groupby("outcome", observed=True).size()
        for endpoint in sorted(set(input_counts.index.astype(str)) | set(model_counts.index.astype(str))):
            EVENT_ROWS.append(
                {
                    "model_id": model_id,
                    "endpoint": endpoint,
                    "input_events": float(input_counts.get(endpoint, 0.0)),
                    "model_events": float(model_counts.get(endpoint, 0.0)),
                    "input_rows": int(input_rows.get(endpoint, 0)),
                    "model_rows": int(model_rows.get(endpoint, 0)),
                }
            )
    else:
        EVENT_ROWS.append(
            {
                "model_id": model_id,
                "endpoint": outcome,
                "input_events": summary["input_events"],
                "model_events": summary["model_events"],
                "input_rows": summary["input_rows"],
                "model_rows": summary["model_rows"],
            }
        )
    SAMPLEFLOW_ROWS.extend(
        [
            {
                "stage": f"{model_id}_input_to_fit",
                "rows": summary["input_rows"],
                "health_regions": summary["input_cluster_counts"]["hrid"],
                "dates": summary["input_cluster_counts"]["day_id"],
                "events": summary["input_events"],
                "note": "Rows passed to the specified fixed-effect Poisson model.",
            },
            {
                "stage": f"{model_id}_retained_by_fit",
                "rows": summary["model_rows"],
                "health_regions": kept_clusters["hrid"],
                "dates": kept_clusters["day_id"],
                "events": summary["model_events"],
                "note": "Estimator retained sample; model_events summed from model._data.",
            },
        ]
    )
    for name in names:
        if family == "service_category":
            display_name = {
                "h_only": "heat_only_vs_neither",
                "s_only": "smoke_only_vs_neither",
                "both": "both_vs_neither",
                "holiday": "holiday_vs_nonholiday",
            }.get(name, name)
        else:
            display_name = name
        RESULT_ROWS.append(
            base_result_row(
                model_id, family, display_name, names, beta, covariance, {name: 1.0}, df, summary
            )
        )
    persist_tables()
    log(
        "FIT_OK", model_id, "rows=", summary["model_rows"], "events=", summary["model_events"],
        "clusters=", kept_clusters, "df=", df, "cov_min_eigen=", diagnostics["min_eigenvalue"],
        "seconds=", summary["elapsed_seconds"],
    )
    return model, names, beta, covariance, model_data, summary


def absolute_contrast(
    model_id: str,
    outcome: str,
    model_data: pd.DataFrame,
    names: list[str],
    beta: np.ndarray,
    covariance: np.ndarray,
    summary: dict[str, Any],
) -> None:
    required = ["h_only", "s_only", "both"]
    if any(term not in names for term in required):
        raise ValueError(f"{model_id}: heat-only contrast requires all category terms")
    holiday_in_model = "holiday" in names
    draw_terms = required + (["holiday"] if holiday_in_model else [])
    positions = [names.index(term) for term in draw_terms]
    point_beta = beta[positions]
    raw_cov = covariance[np.ix_(positions, positions)]
    raw_cov = (raw_cov + raw_cov.T) / 2.0
    raw_eigenvalues, raw_eigenvectors = np.linalg.eigh(raw_cov)
    max_eigen = float(np.max(np.abs(raw_eigenvalues))) if raw_eigenvalues.size else 0.0
    tolerance = max(1e-12, max_eigen * 1e-10)
    psd_eigenvalues = np.maximum(raw_eigenvalues, 0.0)
    psd_cov = (raw_eigenvectors * psd_eigenvalues) @ raw_eigenvectors.T
    psd_cov = (psd_cov + psd_cov.T) / 2.0
    projection_frobenius = float(np.linalg.norm(psd_cov - raw_cov, ord="fro"))

    data = model_data[["stratum", "h_only", "s_only", "both", "holiday", outcome]].copy()
    data["category"] = np.select(
        [data.h_only.eq(1), data.s_only.eq(1), data.both.eq(1)],
        ["heat_only", "smoke_only", "both"],
        default="neither",
    )
    data["cell_key"] = data.category.astype(str) + "_h" + data.holiday.astype(str)
    counts_by_stratum = (
        data.groupby(["stratum", "cell_key"], observed=True)
        .size()
        .unstack("cell_key", fill_value=0)
        .reindex(columns=CELL_KEYS, fill_value=0)
    )
    observed_by_stratum = data.groupby("stratum", observed=True)[outcome].sum()
    counts_by_stratum["observed_events"] = observed_by_stratum.reindex(counts_by_stratum.index)
    profiles = counts_by_stratum.reset_index().groupby(
        CELL_KEYS, as_index=False, observed=True
    )["observed_events"].sum()
    exposure_counts = profiles.loc[:, CELL_KEYS].to_numpy(dtype=float)
    profile_events = profiles["observed_events"].to_numpy(dtype=float)
    if exposure_counts.shape[0] == 0 or np.any(exposure_counts.sum(axis=1) <= 0):
        raise ValueError(f"{model_id}: no retained exposure profiles for absolute contrast")
    if not np.isclose(profile_events.sum(), float(data[outcome].sum()), rtol=0, atol=1e-6):
        raise ValueError(f"{model_id}: exposure-profile aggregation changed observed model events")
    heat0 = CELL_KEYS.index("heat_only_h0")
    heat1 = CELL_KEYS.index("heat_only_h1")

    def difference_batch(draw_batch: np.ndarray) -> np.ndarray:
        bh = draw_batch[:, 0]
        bs = draw_batch[:, 1]
        bb = draw_batch[:, 2]
        bholiday = draw_batch[:, 3] if holiday_in_model else np.zeros(len(draw_batch))
        cat_weights = np.column_stack(
            [np.ones(len(draw_batch)), np.exp(bh), np.exp(bs), np.exp(bb)]
        )
        holiday_weights = np.column_stack([np.ones(len(draw_batch)), np.exp(bholiday)])
        all_weights = (cat_weights[:, :, None] * holiday_weights[:, None, :]).reshape(
            len(draw_batch), -1
        )
        denominator = all_weights @ exposure_counts.T
        if np.any(denominator <= 0) or not np.isfinite(denominator).all():
            raise ValueError(f"{model_id}: invalid factual stratum calibration denominator")
        changed_exposure = (
            exposure_counts[:, heat0][None, :]
            + exposure_counts[:, heat1][None, :] * np.exp(bholiday)[:, None]
        ) * np.expm1(bh)[:, None]
        per_profile_ratio = changed_exposure / denominator
        return per_profile_ratio @ profile_events

    point_difference = float(difference_batch(point_beta[None, :])[0])
    rng = np.random.default_rng(ABSOLUTE_SEED)
    standard_normal = rng.standard_normal((ABSOLUTE_DRAWS, len(draw_terms)))
    # Explicit eigenvalue-zero PSD projection; no implicit RNG covariance repair.
    transform = raw_eigenvectors @ np.diag(np.sqrt(psd_eigenvalues))
    sampled_beta = point_beta + standard_normal @ transform.T
    estimates = np.empty(ABSOLUTE_DRAWS, dtype=float)
    batch_size = 32
    for start in range(0, ABSOLUTE_DRAWS, batch_size):
        stop = min(start + batch_size, ABSOLUTE_DRAWS)
        estimates[start:stop] = difference_batch(sampled_beta[start:stop])
        if start == 0 or stop == ABSOLUTE_DRAWS or stop % 512 == 0:
            log(
                "ABSOLUTE_PROGRESS", model_id, outcome,
                f"{stop}/{ABSOLUTE_DRAWS}",
                "unique_exposure_profiles=", len(profiles),
            )
    psd_low, psd_high = np.quantile(estimates, [0.025, 0.975])

    def difference_at(point: np.ndarray) -> float:
        return float(difference_batch(point[None, :])[0])

    gradient = np.empty(len(point_beta), dtype=float)
    for index in range(len(point_beta)):
        step = max(1e-6, abs(float(point_beta[index])) * 1e-5)
        plus = point_beta.copy()
        minus = point_beta.copy()
        plus[index] += step
        minus[index] -= step
        gradient[index] = (difference_at(plus) - difference_at(minus)) / (2.0 * step)
    raw_delta_variance = float(gradient @ raw_cov @ gradient)
    raw_scale = max(1e-12, float(np.max(np.abs(raw_cov))) * float(gradient @ gradient) * 1e-10)
    if raw_delta_variance < -raw_scale:
        raw_delta_se = None
        raw_delta_low = None
        raw_delta_high = None
        raw_delta_status = "undefined: original covariance gives negative delta-method variance"
    else:
        raw_delta_se = math.sqrt(max(raw_delta_variance, 0.0))
        raw_delta_low = point_difference - 1.96 * raw_delta_se
        raw_delta_high = point_difference + 1.96 * raw_delta_se
        raw_delta_status = "normal delta-method interval using the unmodified covariance"
    observed_model_total = float(data[outcome].sum())
    observed_input_total = float(summary["input_events"])
    if observed_model_total <= 0 or observed_input_total <= 0:
        raise ValueError(f"{model_id}: cannot scale contrast to zero observed outcome total")
    row = {
        "model_id": model_id,
        "outcome": outcome,
        "contrast": "heat-only days reassigned to neither; other categories and holiday retained",
        "observed_input_total": observed_input_total,
        "observed_model_total": observed_model_total,
        "model_rows": int(summary["model_rows"]),
        "retained_hrid_clusters": summary["retained_cluster_counts"]["hrid"],
        "retained_date_clusters": summary["retained_cluster_counts"]["day_id"],
        "point_model_associated_difference": point_difference,
        "PSD_draw_interval_low": float(psd_low),
        "PSD_draw_interval_high": float(psd_high),
        "draws": ABSOLUTE_DRAWS,
        "draw_distribution": "joint multivariate normal using explicitly eigenvalue-zero PSD covariance",
        "random_seed": ABSOLUTE_SEED,
        "raw_delta_method_interval_low": raw_delta_low,
        "raw_delta_method_interval_high": raw_delta_high,
        "raw_delta_method_se": raw_delta_se,
        "raw_delta_method_status": raw_delta_status,
        "raw_covariance_submatrix_min_eigenvalue": float(np.min(raw_eigenvalues)),
        "raw_covariance_submatrix_negative_eigenvalues": int(np.sum(raw_eigenvalues < 0)),
        "PSD_projection_negative_eigenvalues_set_to_zero": int(np.sum(raw_eigenvalues < 0)),
        "PSD_projection_frobenius_change": projection_frobenius,
        "PSD_projection_tolerance": tolerance,
        "PSD_covariance_submatrix": json.dumps(psd_cov.tolist()),
        "raw_covariance_submatrix": json.dumps(raw_cov.tolist()),
        "difference_over_observed_model_total": point_difference / observed_model_total,
        "PSD_fraction_interval_low": float(psd_low / observed_model_total),
        "PSD_fraction_interval_high": float(psd_high / observed_model_total),
        "difference_over_observed_input_total": point_difference / observed_input_total,
        "input_fraction_interval_low": float(psd_low / observed_input_total),
        "input_fraction_interval_high": float(psd_high / observed_input_total),
        "strata_with_observed_model_events": int(observed_by_stratum.gt(0).sum()),
        "unique_stratum_exposure_profiles": int(len(profiles)),
        "stratum_intercept_method": (
            "For each coefficient draw, calibrate each retained stratum intercept under factual exposure "
            "to its observed stratum outcome total; hold it fixed after heat-only is reassigned to neither."
        ),
        "holiday_handling": (
            "Holiday is held at each observed row value and may vary within stratum."
            if holiday_in_model else "No holiday term in the original category specification."
        ),
        "conditioning": "Condition on observed retained stratum totals; parameter uncertainty only.",
        "interpretation": "Model-associated noncausal contrast, not an attributable fraction.",
        "formula": summary["formula"],
    }
    FIT_RECORDS[model_id]["absolute_sampling_diagnostic"] = {
        "coefficient_order": draw_terms,
        "raw_covariance_submatrix": raw_cov.tolist(),
        "raw_submatrix_eigenvalues": raw_eigenvalues.tolist(),
        "PSD_covariance_submatrix": psd_cov.tolist(),
        "PSD_eigenvalues": psd_eigenvalues.tolist(),
        "negative_eigenvalues_set_to_zero_explicitly": int(np.sum(raw_eigenvalues < 0)),
        "projection_frobenius_norm": projection_frobenius,
        "point_contrast_gradient": gradient.tolist(),
        "raw_delta_method_variance": raw_delta_variance,
        "raw_delta_method_interval": [raw_delta_low, raw_delta_high],
        "PSD_draw_interval": [float(psd_low), float(psd_high)],
    }
    ABSOLUTE_ROWS.append(row)
    persist_tables()
    log(
        "ABSOLUTE_OK", model_id, outcome,
        "difference=", point_difference,
        "PSD_interval=", [float(psd_low), float(psd_high)],
        "raw_interval=", [raw_delta_low, raw_delta_high],
        "fraction_retained_total=", point_difference / observed_model_total,
    )


def run_service_models(service: pd.DataFrame) -> None:
    specs = [
        ("service_original_two_way", False, False, True),
        ("service_holiday_two_way", True, False, True),
        ("service_holiday_excluding_2020_2021", True, True, False),
    ]
    for spec_id, include_holiday, exclude_covid, calculate_absolute in specs:
        frame = service.copy()
        if exclude_covid:
            frame = frame.loc[frame.date.dt.year.isin([2018, 2019, 2022])].copy()
        regressors = ["h_only", "s_only", "both"] + (["holiday"] if include_holiday else [])
        for outcome in SERVICE_OUTCOMES:
            model_id = f"{spec_id}_{outcome}"
            formula = f"{outcome} ~ {' + '.join(regressors)} | stratum"
            model, names, beta, covariance, model_data, summary = fit_model(
                model_id, "service_category", formula, frame, outcome, regressors
            )
            persist_tables()
            if calculate_absolute:
                absolute_contrast(model_id, outcome, model_data, names, beta, covariance, summary)
            del model, model_data


def run_drainage_models(lagged: pd.DataFrame) -> None:
    lag_terms = [f"t{lag}" for lag in range(8)] + [f"f{lag}" for lag in range(8)]
    specs = [
        ("drainage_continuous_two_way", False, False),
        ("drainage_continuous_holiday_two_way", True, False),
        ("drainage_continuous_holiday_excluding_2020_2021", True, True),
    ]
    for model_id, include_holiday, exclude_covid in specs:
        frame = lagged
        if exclude_covid:
            frame = lagged.loc[lagged.date.dt.year.isin([2018, 2019, 2022])].copy()
        terms = lag_terms + (["holiday"] if include_holiday else [])
        formula = f"any_drainage ~ {' + '.join(terms)} | stratum"
        model, names, beta, covariance, model_data, summary = fit_model(
            model_id, "drainage_continuous_lag0_7", formula, frame, "any_drainage", terms
        )
        add_estimate_row(
            model_id, "drainage_continuous_lag0_7", "cumulative_temperature_lag0_7_per_1C",
            names, beta, covariance, {f"t{lag}": 1.0 for lag in range(8)}, summary["df"], summary,
        )
        add_estimate_row(
            model_id, "drainage_continuous_lag0_7", "cumulative_fire_lag0_7_per_10ug_m3",
            names, beta, covariance, {f"f{lag}": 1.0 for lag in range(8)}, summary["df"], summary,
        )
        persist_tables()
        del model, model_data


def make_stacked(frame: pd.DataFrame, phenotype: str, exclude_covid: bool) -> pd.DataFrame:
    endpoint_columns = {
        "broad": ["broad_elective", "broad_urgent"],
        "strict": ["strict_elective", "strict_urgent"],
    }[phenotype]
    selected = frame
    if exclude_covid:
        selected = frame.loc[frame.date.dt.year.isin([2018, 2019, 2022])]
    long = selected[
        ["hrid", "date", "day_id", "stratum", "group", "holiday", *endpoint_columns]
    ].melt(
        id_vars=["hrid", "date", "day_id", "stratum", "group", "holiday"],
        value_vars=endpoint_columns,
        var_name="endpoint_source",
        value_name="count",
    )
    source_to_outcome = {
        f"{phenotype}_elective": "elective",
        f"{phenotype}_urgent": "urgent",
    }
    long["outcome"] = pd.Categorical(
        long.endpoint_source.map(source_to_outcome), categories=["elective", "urgent"]
    )
    long["endpoint_stratum"] = long.stratum + "_" + long.outcome.astype(str)
    if long.outcome.isna().any() or long["count"].isna().any():
        raise ValueError(f"{phenotype} urgent/elective stack has invalid records")
    return long


def run_joint_models(joint: pd.DataFrame) -> None:
    for phenotype in ["broad", "strict"]:
        for exclude_covid in [False, True]:
            suffix = "_excluding_2020_2021" if exclude_covid else ""
            model_id = f"{phenotype}_urgent_elective_holiday_two_way{suffix}"
            frame = make_stacked(joint, phenotype, exclude_covid)
            formula = "count ~ C(outcome)*C(group) + C(outcome):holiday | stratum^outcome"
            expected = [
                "C(group)[T.heat_only]",
                "C(group)[T.smoke_only]",
                "C(group)[T.both]",
                "C(outcome)[T.urgent]:C(group)[T.heat_only]",
                "C(outcome)[T.urgent]:C(group)[T.smoke_only]",
                "C(outcome)[T.urgent]:C(group)[T.both]",
                "C(outcome)[elective]:holiday",
                "C(outcome)[urgent]:holiday",
            ]
            model, names, beta, covariance, model_data, summary = fit_model(
                model_id, f"{phenotype}_urgent_elective_joint", formula, frame, "count", expected
            )
            add_estimate_row(
                model_id,
                f"{phenotype}_urgent_elective_joint",
                "heat_only_urgent_RR_over_elective_RR",
                names,
                beta,
                covariance,
                {"C(outcome)[T.urgent]:C(group)[T.heat_only]": 1.0},
                summary["df"],
                summary,
            )
            persist_tables()
            del model, model_data, frame


def preserve_legacy_joint_results() -> None:
    legacy = pd.read_csv(LEGACY_RESULTS_PATH)
    target = legacy.loc[
        legacy.contrast.eq("heat_only_urgent_RR_over_elective_RR")
        & legacy.model.isin(
            [
                "broad_N20_N23_urgent_elective_joint_categorical",
                "strict_N20_urgent_elective_joint_categorical",
            ]
        )
    ]
    if len(target) != 2:
        raise ValueError(f"Expected two existing legacy joint contrasts, found {len(target)}")
    for row in target.to_dict("records"):
        family = "broad_urgent_elective_joint" if row["model"].startswith("broad_") else "strict_urgent_elective_joint"
        RESULT_ROWS.append(
            {
                "model_id": row["model"] + "_legacy",
                "model_family": family,
                "estimand": row["contrast"],
                "beta": float(np.log(float(row["RR"]))),
                "se": None,
                "RR_or_ratio": float(row["RR"]),
                "CI_low": float(row["CI_low"]),
                "CI_high": float(row["CI_high"]),
                "statistic_t": None,
                "p": float(row["p"]),
                "inference": "legacy existing normal inference",
                "df": None,
                "contrast_terms": row.get("contrast_terms"),
                "input_rows": None,
                "model_rows": None,
                "input_events": float(row["input_events"]),
                "model_events": float(row["model_events"]),
                "retained_hrid_clusters": None,
                "retained_date_clusters": None,
                "formula": row["formula"],
                "cluster": row["cluster"],
                "source_status": "legacy_existing_result_not_refitted",
                "note": "Prior no-holiday health-region-only fit preserved as legacy; not revised inference.",
            }
        )
    persist_tables()


def write_methods_and_manifest(input_audits: dict[str, Any], status: str) -> None:
    methods = """# revision 巴西服务模型重分析方法

## 预定模型与样本

仅分析 2018–2022 年巴西 450 个卫生区的每日数据。模型输入沿用既有的完整区域日暴露面板、已完成的急诊与手术日汇总、2022 年 12 月 IPCA 固定价指数，以及既有国家假日历。新模型不扩展国家或病种，也不触碰 2008 年主模型。模型规格在查看本轮结果前写入 work/revision_20261113/service/analysis_specification.json，不按 p 值选择规格。

严格 N20 急诊定义为 IDENT=1、主诊断 N20、CAR_INT=02。床日按同一批入院记录的行政字段汇总，作为入院归属床日；金额按申报月份乘既有 IPCA 因子折为 2022 年 BRL 后按卫生区和入院日求和。零病例区域日显式补零。

固定效应为卫生区 × 年月 × 星期几。热、烟暴露四类别以 neither 为参考，分别纳入 heat-only、smoke-only、both 指示项。急诊数、床日、2022BRL 报销总额均估计三组模型：原类别规格；加入既有国家 holiday 的修订主规格；加入 holiday 并仅保留 2018–2019、2022 入院日期的疫情敏感性。前两组均有绝对关联差异计算，第三组不计算该差异。

引流结局使用已链接的 any_drainage 区域日计数，连续纳入日均温度和火源 PM2.5 各自 lag 0–7 项（火源 PM2.5 按每 10 µg/m³ 缩放）。共同样本先去除各区最初 7 天形成 lag，再估计无假日、加 holiday、及加 holiday 且排除 2020–2021 入院日的模型。温度与火源累计对比分别使用相应 8 个 lag 系数之和及完整协方差。疫情敏感性按结局日期排除 2020–2021，2022 年初的滞后暴露仍沿原始日历计算。

急/择期 joint 模型将 N20–N23 广义或严格 N20 每日急诊和择期计数堆叠；固定效应按卫生区 × 年月 × 星期几 × 结局分别设置。主对比为热-only 日急诊 RR 与择期 RR 之比。仅对其重估 holiday 双向聚类和疫情敏感性；既有无 holiday、卫生区单向聚类数值原样标为 legacy。

## 推断

所有本轮模型一次拟合，协方差采用卫生区与日历日期双向 CRV1。每个结果使用该拟合样本内实际保留的卫生区簇数和日期簇数，df = min(G_region, G_date) − 1，并按 Student t 分布计算区间和 p 值。模型公式、输入/保留行数、输入/保留事件数、保留簇数、完整系数和原始协方差矩阵、协方差特征值诊断均保存。相对效应是条件于卫生区 × 年月 × 星期几固定效应的模型关联比率，不是总体边际效应或因果效应。

## 绝对关联对比

仅对原类别与加 holiday 的核心服务模型，各结局从拟合的暴露/假日系数协方差子矩阵抽取 4096 次联合参数样本。每次抽样按事实暴露计算每个保留层的截距，使事实情景期望和等于观测层总量；随后把 heat-only 记录改为 neither，保留其他暴露和每行 holiday 指示，并固定该次已校准截距。假日可在同一固定效应层内变化，且在反事实设计矩阵中逐行保留。

若用于抽样的原始协方差子矩阵不是半正定，保留原始矩阵并显式特征分解，将负特征值置零后构造 PSD 抽样矩阵；绝不依赖 NumPy 隐式修复。绝对对比同时报告原始协方差的数值 delta-method 区间（若标量方差非负）和 PSD 抽样区间。PSD 区间是指定投影下的参数不确定性近似；原始 delta 区间是对未改矩阵的局部诊断，二者差别与投影幅度一并保留。所有差异另报占保留模型事件数和输入观察总量比例。上述结果条件于观测层总量，不含基线计数不确定性；是非因果模型关联对比，不是 PAF。

## 软件与资源

使用 Python 3.12、DuckDB、Pandas、NumPy、SciPy 与 pyfixest；DuckDB 和每个进程的 OpenMP/BLAS/Numba 线程均限为 4。使用 bio-clinical-biostatistics-effect-measures 与 bio-clinical-biostatistics-logistic-regression 技能对条件效应、总体边际效应和相对效应作区分。资源盘点显示有效 CPU 容量为 16，磁盘空闲约 313 GiB；首次检测 RAM 可用量为 unknown，后续主控资源检查约为 48 GiB 总量、24 GiB 空闲。
"""
    (OUT / "methods_zh.md").write_text(methods, encoding="utf-8")
    save_json(
        OUT / "analysis_manifest.json",
        {
            "status": status,
            "started_at_utc": STARTED,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": str(WORK_DIR / "fit_service_models.py"),
            "specification": str(WORK_DIR / "analysis_specification.json"),
            "python_version": sys.version,
            "package_versions": {
                "duckdb": duckdb.__version__,
                "pandas": pd.__version__,
                "numpy": np.__version__,
                "pyfixest": pf.__version__,
            },
            "thread_limits": {
                "duckdb": 4,
                "OMP": os.environ.get("OMP_NUM_THREADS"),
                "MKL": os.environ.get("MKL_NUM_THREADS"),
                "NUMBA": os.environ.get("NUMBA_NUM_THREADS"),
                "OPENBLAS": os.environ.get("OPENBLAS_NUM_THREADS"),
            },
            "inputs": input_audits,
            "outputs": {
                "results": str(OUT / "results.csv"),
                "absolute_contrasts": str(OUT / "absolute_service_contrasts.csv"),
                "sample_audit": str(OUT / "sample_audit.csv"),
                "sampleflow": str(OUT / "sampleflow.csv"),
                "model_events": str(OUT / "model_events.csv"),
                "model_fits_and_covariances": str(OUT / "model_fits.json"),
                "model_summary": str(OUT / "model_summary.csv"),
                "methods": str(OUT / "methods_zh.md"),
            },
        },
    )


def write_status(status: str, **fields: Any) -> None:
    save_json(
        STATUS_PATH,
        {
            "status": status,
            "pid": os.getpid(),
            "started_at_utc": STARTED,
            "updated_at_utc": datetime.now(timezone.utc).isoformat(),
            "completed_models": list(FIT_RECORDS),
            "absolute_contrasts_completed": len(ABSOLUTE_ROWS),
            **fields,
        },
    )


def main() -> None:
    write_status("running", current_step="load_inputs")
    panel, panel_audit = prepare_base_panel()
    service, service_audit = prepare_service_panel(panel)
    drainage, drainage_audit = prepare_drainage_panel(panel)
    joint, joint_audit = prepare_urgent_elective_panel(panel)
    input_audits = {
        "base_panel_and_holiday": panel_audit,
        "strict_N20_service": service_audit,
        "SP_any_drainage": drainage_audit,
        "urgent_elective_joint": joint_audit,
    }
    write_status("running", current_step="service_models", input_audits=input_audits)
    run_service_models(service)
    del service
    write_status("running", current_step="drainage_models", input_audits=input_audits)
    run_drainage_models(drainage)
    del drainage
    write_status("running", current_step="urgent_elective_models", input_audits=input_audits)
    run_joint_models(joint)
    del joint
    write_status("running", current_step="preserve_legacy_results", input_audits=input_audits)
    preserve_legacy_joint_results()
    persist_tables()
    write_methods_and_manifest(input_audits, "complete")
    write_status("complete", input_audits=input_audits)
    log("COMPLETE", OUT)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        failure = {
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "failed_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        try:
            write_status("failed", **failure)
        finally:
            persist_tables()
        log("FAILED", repr(exc))
        raise
