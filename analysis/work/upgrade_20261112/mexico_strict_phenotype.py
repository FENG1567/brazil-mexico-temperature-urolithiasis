"""Compare broad N20-N23 and strict N20 lag effects in Mexico and Brazil.

Uses only pre-existing de-identified state/day sources and writes the
mexico_strict_* outputs. No row-level identifiers are loaded or printed.
"""
import os
os.environ["OMP_NUM_THREADS"] = "8"
os.environ["NUMBA_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, "F:/analysis/climate_biliary_20261102/work/pylibs")

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import norm, t as student_t


BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "work/upgrade_20261112"
OUT = BASE / "outputs/upgrade_20261112"
SOURCE_WORK = Path("F:/analysis/climate_biliary_20261102/work")
BRAZIL_PANEL_PATH = SOURCE_WORK / "panel_composite_humidity.parquet"
BRAZIL_CLINICAL_PATH = WORK / "clinical_code_day_aggregate.parquet"
MEXICO_PANEL_PATH = SOURCE_WORK / "mexico_composite_state_day_panel.parquet"
MEXICO_CLINICAL_PATH = SOURCE_WORK / "mexico_n20_proced.parquet"
MEXICO_NASA_PATH = SOURCE_WORK / "mexico_state_day_temperature.parquet"
MEXICO_ERA5_PATH = Path(
    "D:/数据库集合/新增数据/ERA5_crossnational/mex_state_day_2015_2023.csv"
)
WINDOW_START = "2018-01-01"
WINDOW_END_EXCLUSIVE = "2023-01-01"
WINDOW_DAYS = 1826
EXPECTED = {
    "brazil_broad": 431_374,
    "brazil_strict_all": 321_201,
    "brazil_strict_urgent": 226_402,
    "brazil_strict_elective": 94_799,
    "mexico_broad": 40_415,
    "mexico_strict": 32_417,
}

OUT.mkdir(parents=True, exist_ok=True)
db = duckdb.connect()
db.execute("SET threads=8")

RESULTS = []
SUMMARIES = []
COUNTS = []
FULL_MODELS = {}
SOURCE_COVERAGE = {}


def sql_path(path):
    return Path(path).as_posix().replace("'", "''")


def parquet_columns(path):
    return [row[0] for row in db.execute(
        f"DESCRIBE SELECT * FROM read_parquet('{sql_path(path)}')"
    ).fetchall()]


def persist():
    pd.DataFrame(RESULTS).to_csv(
        OUT / "mexico_strict_model_results.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(SUMMARIES).to_csv(
        OUT / "mexico_strict_model_summary.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(COUNTS).to_csv(
        OUT / "mexico_strict_endpoint_counts_by_state.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (OUT / "mexico_strict_full_covariance.json").write_text(
        json.dumps(
            {
                "window": "2018-01-01 through 2022-12-31; filtered before lags",
                "contrast": "sum of temperature coefficients for lags 0-7",
                "models": FULL_MODELS,
            },
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )


def assert_complete_state_calendar(frame, country, expected_states):
    if frame.duplicated(["state", "date"]).any():
        raise ValueError(f"Duplicate {country} state/date rows")
    states = sorted(frame.state.astype(str).unique().tolist())
    if len(states) != expected_states:
        raise ValueError(f"Expected {expected_states} {country} states, found {len(states)}")
    dates_per_state = frame.groupby("state", observed=True).date.agg(["min", "max", "nunique"])
    if not dates_per_state["nunique"].eq(WINDOW_DAYS).all():
        raise ValueError(f"Incomplete {country} calendar in 2018-2022")
    if not dates_per_state["min"].eq(pd.Timestamp(WINDOW_START)).all():
        raise ValueError(f"Unexpected {country} start date")
    if not dates_per_state["max"].eq(pd.Timestamp("2022-12-31")).all():
        raise ValueError(f"Unexpected {country} end date")
    diffs = frame.sort_values(["state", "date"]).groupby("state", observed=True).date.diff().dropna()
    if not diffs.dt.days.eq(1).all():
        raise ValueError(f"Nonconsecutive {country} calendar; refusing to create lags")


def add_state_counts(country, frame, endpoints):
    for endpoint, label in endpoints:
        grouped = frame.groupby("state", observed=True)
        agg = grouped[endpoint].agg(admissions="sum")
        days = grouped[endpoint].apply(lambda values: int((values > 0).sum()))
        total_days = grouped.size()
        for state in agg.index:
            COUNTS.append(
                {
                    "country": country,
                    "endpoint": label,
                    "state": str(state),
                    "state_days": int(total_days[state]),
                    "days_with_admissions": int(days[state]),
                    "admissions": int(agg.loc[state, "admissions"]),
                }
            )


def fit_lag_model(country, product, phenotype, outcome_column, outcome_label,
                  state_daily, full_window_events):
    key = f"{country.lower()}_{product}_{phenotype}_lag07"
    lag_terms = [f"t{lag}" for lag in range(8)]
    model_frame = state_daily[["state", "date", "stratum", *lag_terms, outcome_column]].copy()
    model_frame = model_frame.dropna(subset=lag_terms).copy()
    if len(model_frame) == 0:
        raise ValueError(f"Empty lag model sample: {key}")
    formula = outcome_column + " ~ " + " + ".join(lag_terms) + " | stratum"
    fe_sum = model_frame.groupby("stratum", observed=True)[outcome_column].sum()
    fe_size = model_frame.groupby("stratum", observed=True).size()
    zero_levels = fe_sum.index[fe_sum.eq(0)]
    zero_fe_rows = int(fe_size.reindex(zero_levels).sum())
    zero_fe_cells = int(len(zero_levels))
    input_n = len(model_frame)
    clusters = int(model_frame.state.nunique())
    input_events = float(model_frame[outcome_column].sum())
    print(
        f"START {key} input_n={input_n} input_events={int(input_events)} "
        f"zero_fe={zero_fe_cells} formula={formula}",
        flush=True,
    )
    started = time.time()
    model = pf.fepois(formula, data=model_frame, vcov={"CRV1": "state"})
    coef = model.coef()
    names = [str(name) for name in coef.index]
    values = coef.to_numpy(dtype=float)
    covariance = np.asarray(model._vcov, dtype=float)
    if covariance.shape != (len(names), len(names)):
        raise ValueError(f"Unexpected covariance shape in {key}: {covariance.shape}")
    weights = np.array([1.0 if name in lag_terms else 0.0 for name in names])
    if not weights.any():
        raise ValueError(f"No lag coefficients found in {key}: {names}")
    beta = float(weights @ values)
    variance = float(weights @ covariance @ weights)
    if variance < 0 and abs(variance) < 1e-12:
        variance = 0.0
    if variance < 0:
        raise ValueError(f"Negative cumulative variance in {key}: {variance}")
    se = float(np.sqrt(variance))
    z = beta / se if se > 0 else float("inf")
    p_normal = float(2 * norm.sf(abs(z))) if np.isfinite(z) else 0.0
    df = clusters - 1
    if df <= 0:
        raise ValueError(f"Insufficient state clusters in {key}: {clusters}")
    tcrit = float(student_t.ppf(0.975, df))
    p_t = float(2 * student_t.sf(abs(z), df)) if np.isfinite(z) else 0.0
    model_n = int(model._N)
    dropped_n = input_n - model_n
    summary = {
        "model": key,
        "country": country,
        "exposure_product": product,
        "phenotype": phenotype,
        "diagnosis_definition": outcome_label,
        "formula": formula,
        "fixed_effects": "state x calendar year-month x day-of-week stratum",
        "cluster": "CRV1 by state",
        "state_clusters": clusters,
        "t_df": df,
        "input_n": input_n,
        "model_n": model_n,
        "rows_dropped": dropped_n,
        "input_strata": int(model_frame.stratum.nunique()),
        "zero_fe_strata": zero_fe_cells,
        "zero_fe_rows_expected_dropped": zero_fe_rows,
        "other_rows_dropped_by_fit": dropped_n - zero_fe_rows,
        "events_full_window": int(full_window_events),
        "events_after_lag_window": int(input_events),
        "model_events": int(input_events) if dropped_n == zero_fe_rows else None,
        "seconds": round(time.time() - started, 2),
    }
    normal_low = float(np.exp(beta - 1.96 * se))
    normal_high = float(np.exp(beta + 1.96 * se))
    t_low = float(np.exp(beta - tcrit * se))
    t_high = float(np.exp(beta + tcrit * se))
    result = {
        **summary,
        "contrast": "cumulative_temperature_lag_0_7_per_C",
        "beta_sum": beta,
        "se_full_covariance": se,
        "RR_per_C": float(np.exp(beta)),
        "normal_CI_low": normal_low,
        "normal_CI_high": normal_high,
        "normal_p": p_normal,
        "t_G_minus_1_CI_low": t_low,
        "t_G_minus_1_CI_high": t_high,
        "t_G_minus_1_p": p_t,
        "contrast_terms": json.dumps(lag_terms),
    }
    RESULTS.append(result)
    SUMMARIES.append(summary)
    FULL_MODELS[key] = {
        **summary,
        "coefficient_names": names,
        "coefficients": values.tolist(),
        "covariance": covariance.tolist(),
        "cumulative_temperature_contrast": {
            "weights": weights.tolist(),
            "beta_sum": beta,
            "variance": variance,
            "se": se,
            "RR_per_C": result["RR_per_C"],
            "normal_CI": [normal_low, normal_high],
            "normal_p": p_normal,
            "t_G_minus_1_CI": [t_low, t_high],
            "t_G_minus_1_p": p_t,
        },
    }
    print(
        f"RESULT {key} RR={result['RR_per_C']:.9f} "
        f"normalCI=[{normal_low:.9f},{normal_high:.9f}] "
        f"tCI=[{t_low:.9f},{t_high:.9f}] p_t={p_t:.8g} "
        f"n={model_n} clusters={clusters}",
        flush=True,
    )
    persist()


def attach_lags_and_fit(country, product, exposure_frame, state_daily):
    exposure_frame = exposure_frame[["state", "date", "temperature"]].copy()
    exposure_frame["state"] = exposure_frame.state.astype(str).str.zfill(2)
    exposure_frame["date"] = pd.to_datetime(exposure_frame.date).dt.normalize()
    if exposure_frame.duplicated(["state", "date"]).any():
        raise ValueError(f"Duplicate {product} exposure state-days")
    if len(exposure_frame) != len(state_daily):
        raise ValueError(f"Exposure sample row count differs for {country}/{product}")
    joined = state_daily.merge(
        exposure_frame, on=["state", "date"], how="left", validate="one_to_one"
    )
    if joined.temperature.isna().any():
        raise ValueError(f"Missing {product} exposure on matched {country} state-days")
    joined = joined.sort_values(["state", "date"]).reset_index(drop=True)
    if joined.groupby("state", observed=True).date.diff().dropna().dt.days.ne(1).any():
        raise ValueError(f"Nonconsecutive dates for {country}/{product}; refusing lags")
    for lag in range(8):
        joined[f"t{lag}"] = joined.groupby("state", sort=False).temperature.shift(lag)
    product_frame = joined.drop(columns=["temperature"])
    full_events = {
        "broad_N20_N23": int(product_frame.broad_N20_N23.sum()),
        "strict_N20": int(product_frame.strict_N20.sum()),
    }
    for phenotype, outcome, label in [
        ("broad_N20_N23", "broad_N20_N23", "all admissions with principal diagnosis prefix N20-N23"),
        ("strict_N20", "strict_N20", "all admissions with principal diagnosis prefix N20"),
    ]:
        fit_lag_model(
            country, product, phenotype, outcome, label,
            product_frame, full_events[phenotype],
        )


def load_mexico():
    panel_path = sql_path(MEXICO_PANEL_PATH)
    panel = db.sql(
        f"""
        SELECT state, date, N20 AS broad_N20_N23, stratum
        FROM read_parquet('{panel_path}')
        WHERE date >= TIMESTAMP '{WINDOW_START}' AND date < TIMESTAMP '{WINDOW_END_EXCLUSIVE}'
        """
    ).df()
    panel["state"] = panel.state.astype(str).str.zfill(2)
    panel["date"] = pd.to_datetime(panel.date).dt.normalize()
    panel["broad_N20_N23"] = panel.broad_N20_N23.astype("int64")
    assert_complete_state_calendar(panel, "Mexico", 32)
    expected_stratum = (
        panel.state + "_" + panel.date.dt.strftime("%Y%m") + "_" + panel.date.dt.dayofweek.astype(str)
    )
    if not expected_stratum.equals(panel.stratum.astype(str)):
        raise ValueError("Mexico panel strata differ from state x year-month x weekday definition")

    clinical_path = sql_path(MEXICO_CLINICAL_PATH)
    clinical = db.sql(
        f"""
        SELECT LPAD(CAST(TRY_CAST(state AS INTEGER) AS VARCHAR), 2, '0') AS state,
               CAST(date AS DATE) AS date,
               COUNT(*) AS broad_source,
               SUM(CASE WHEN SUBSTR(UPPER(TRIM(COALESCE(AFECPRIN,''))),1,3)='N20'
                        THEN 1 ELSE 0 END) AS strict_N20
        FROM read_parquet('{clinical_path}')
        WHERE date >= DATE '{WINDOW_START}' AND date < DATE '{WINDOW_END_EXCLUSIVE}'
          AND TRY_CAST(state AS INTEGER) BETWEEN 1 AND 32
          AND SUBSTR(UPPER(TRIM(COALESCE(AFECPRIN,''))),1,3) IN ('N20','N21','N22','N23')
        GROUP BY 1,2
        """
    ).df()
    clinical["date"] = pd.to_datetime(clinical.date).dt.normalize()
    if clinical.duplicated(["state", "date"]).any():
        raise ValueError("Duplicate aggregated Mexico clinical state-days")
    joined = panel.merge(clinical, on=["state", "date"], how="left", validate="one_to_one")
    joined[["broad_source", "strict_N20"]] = joined[["broad_source", "strict_N20"]].fillna(0).astype("int64")
    mismatch = joined.broad_source.ne(joined.broad_N20_N23)
    if mismatch.any():
        raise ValueError(
            f"Mexico no-ID clinical source differs from broad panel on {int(mismatch.sum())} state-days"
        )
    if int(joined.broad_source.sum()) != EXPECTED["mexico_broad"]:
        raise ValueError(f"Unexpected Mexico broad count: {int(joined.broad_source.sum())}")
    if int(joined.strict_N20.sum()) != EXPECTED["mexico_strict"]:
        raise ValueError(f"Unexpected Mexico strict N20 count: {int(joined.strict_N20.sum())}")
    joined["strict_N20"] = joined.strict_N20.astype("int64")
    add_state_counts(
        "Mexico", joined,
        [("broad_N20_N23", "broad N20-N23"), ("strict_N20", "strict N20 all admissions")],
    )

    era5_path = sql_path(MEXICO_ERA5_PATH)
    era5 = db.sql(
        f"""
        SELECT LPAD(CAST(state AS VARCHAR), 2, '0') AS state,
               MAKE_DATE(CAST(year AS INTEGER), CAST(month AS INTEGER), CAST(day AS INTEGER)) AS date,
               CAST(tmean AS DOUBLE) AS temperature
        FROM read_csv_auto('{era5_path}', HEADER=TRUE, SAMPLE_SIZE=-1)
        WHERE year BETWEEN 2018 AND 2022
        """
    ).df()
    era5["date"] = pd.to_datetime(era5.date).dt.normalize()
    if len(era5) != len(joined) or era5.temperature.isna().any():
        raise ValueError("Mexico ERA5-Land exposure does not fully cover the state-day sample")
    if era5.duplicated(["state", "date"]).any():
        raise ValueError("Duplicate Mexico ERA5-Land state-day exposures")

    nasa_path = sql_path(MEXICO_NASA_PATH)
    nasa = db.sql(
        f"""
        SELECT LPAD(CAST(state AS VARCHAR), 2, '0') AS state,
               CAST(date AS DATE) AS date,
               CAST(tmean AS DOUBLE) AS temperature
        FROM read_parquet('{nasa_path}')
        WHERE date >= DATE '{WINDOW_START}' AND date < DATE '{WINDOW_END_EXCLUSIVE}'
        """
    ).df()
    nasa["date"] = pd.to_datetime(nasa.date).dt.normalize()
    if len(nasa) != len(joined) or nasa.temperature.isna().any():
        raise ValueError("Mexico NASA POWER exposure does not fully cover the state-day sample")
    if nasa.duplicated(["state", "date"]).any():
        raise ValueError("Duplicate Mexico NASA POWER state-day exposures")
    nasa_panel = db.sql(
        f"""SELECT state, CAST(date AS DATE) AS date, CAST(tmean AS DOUBLE) AS panel_tmean
            FROM read_parquet('{panel_path}')
            WHERE date >= DATE '{WINDOW_START}' AND date < DATE '{WINDOW_END_EXCLUSIVE}'
            ORDER BY state,date"""
    ).df()
    nasa_panel["date"] = pd.to_datetime(nasa_panel.date).dt.normalize()
    nasa_cmp = nasa.merge(nasa_panel, on=["state", "date"], how="left", validate="one_to_one")
    if nasa_cmp.panel_tmean.isna().any() or not np.allclose(
        nasa_cmp.temperature.to_numpy(), nasa_cmp.panel_tmean.to_numpy(), rtol=0, atol=1e-8
    ):
        raise ValueError("NASA POWER fallback differs from the historical Mexico panel tmean")

    SOURCE_COVERAGE["mexico"] = {
        "historical_panel": str(MEXICO_PANEL_PATH),
        "clinical_source": str(MEXICO_CLINICAL_PATH),
        "clinical_source_schema": parquet_columns(MEXICO_CLINICAL_PATH),
        "clinical_source_contains_personal_id": False,
        "clinical_source_script": "extract_mexico_proced.py; diagnosis-prefix filter N20-N23; no ID column",
        "window": "2018-01-01 through 2022-12-31",
        "panel_rows": len(joined),
        "state_count": int(joined.state.nunique()),
        "days_per_state": WINDOW_DAYS,
        "broad_N20_N23_total": int(joined.broad_N20_N23.sum()),
        "strict_N20_total": int(joined.strict_N20.sum()),
        "broad_daily_match_days": int((~mismatch).sum()),
        "broad_daily_mismatch_days": int(mismatch.sum()),
        "ERA5_Land": {
            "source": str(MEXICO_ERA5_PATH),
            "coverage_rows": len(era5),
            "matched_rows": len(joined),
            "missing_temperature": int(era5.temperature.isna().sum()),
            "source_definition": "CDS reanalysis-era5-land; hourly 2m temperature converted K to C, daily mean, cos(latitude)-area-weighted grid-cell mean within state polygons",
            "source_script": "build_mex_state_day.py",
            "actual_file_coverage": "2015-01-01 through 2022-12-31; filename says 2015_2023 but 2023 files are absent",
        },
        "NASA_POWER_supplement": {
            "source": str(MEXICO_NASA_PATH),
            "coverage_rows": len(nasa),
            "matched_rows": len(joined),
            "missing_temperature": int(nasa.temperature.isna().sum()),
            "source_definition": "NASA POWER daily T2M at state-capital coordinates",
            "source_script": "build_mexico_temperature.py",
            "matches_historical_panel_tmean": True,
        },
    }
    return joined, era5, nasa


def load_brazil():
    panel_path = sql_path(BRAZIL_PANEL_PATH)
    panel = db.sql(
        f"""
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               CAST(N20_N23 AS DOUBLE) AS broad_N20_N23,
               CAST(tmean AS DOUBLE) AS tmean,
               COALESCE(CAST(pop AS DOUBLE),0) AS pop
        FROM read_parquet('{panel_path}')
        WHERE date >= DATE '{WINDOW_START}' AND date < DATE '{WINDOW_END_EXCLUSIVE}'
        """
    ).df()
    panel["hrid"] = panel.hrid.astype(str)
    panel["state"] = panel.hrid.str[:2]
    panel["date"] = pd.to_datetime(panel.date).dt.normalize()
    if panel.duplicated(["hrid", "date"]).any():
        raise ValueError("Duplicate Brazil health-region-day input rows")

    clinical_path = sql_path(BRAZIL_CLINICAL_PATH)
    clinical = db.sql(
        f"""
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23')
                        THEN admissions ELSE 0 END) AS broad_source,
               SUM(CASE WHEN diagnosis='N20' AND CAR_INT IN ('01','02')
                        THEN admissions ELSE 0 END) AS strict_N20,
               SUM(CASE WHEN diagnosis='N20' AND CAR_INT='02'
                        THEN admissions ELSE 0 END) AS strict_N20_urgent,
               SUM(CASE WHEN diagnosis='N20' AND CAR_INT='01'
                        THEN admissions ELSE 0 END) AS strict_N20_elective
        FROM read_parquet('{clinical_path}')
        WHERE date >= DATE '{WINDOW_START}' AND date < DATE '{WINDOW_END_EXCLUSIVE}'
        GROUP BY 1,2
        """
    ).df()
    clinical["hrid"] = clinical.hrid.astype(str)
    clinical["date"] = pd.to_datetime(clinical.date).dt.normalize()
    if clinical.duplicated(["hrid", "date"]).any():
        raise ValueError("Duplicate Brazil clinical health-region-days")
    clinical_type_summary = db.sql(
        f"""
        SELECT SUM(CASE WHEN diagnosis='N20' THEN admissions ELSE 0 END) AS strict_all_types,
               SUM(CASE WHEN diagnosis='N20' AND CAR_INT IN ('01','02') THEN admissions ELSE 0 END) AS strict_01_02,
               COUNT(DISTINCT CASE WHEN diagnosis='N20' THEN CAR_INT END) AS n20_admission_type_levels
        FROM read_parquet('{clinical_path}')
        WHERE date >= DATE '{WINDOW_START}' AND date < DATE '{WINDOW_END_EXCLUSIVE}'
        """
    ).df().iloc[0]
    if int(clinical_type_summary.strict_all_types) != int(clinical_type_summary.strict_01_02):
        raise ValueError("N20 records exist outside CAR_INT 01/02; cannot label 01+02 as all admissions")
    clinical_join = panel[["hrid", "date", "broad_N20_N23"]].merge(
        clinical, on=["hrid", "date"], how="left", validate="one_to_one"
    )
    clinical_cols = ["broad_source", "strict_N20", "strict_N20_urgent", "strict_N20_elective"]
    clinical_join[clinical_cols] = clinical_join[clinical_cols].fillna(0).astype("int64")
    mismatch = clinical_join.broad_source.ne(clinical_join.broad_N20_N23)
    if mismatch.any():
        raise ValueError(
            f"Brazil clinical aggregate differs from broad source on {int(mismatch.sum())} HR-days"
        )
    if int(clinical_join.broad_source.sum()) != EXPECTED["brazil_broad"]:
        raise ValueError(f"Unexpected Brazil broad count: {int(clinical_join.broad_source.sum())}")
    if int(clinical_join.strict_N20.sum()) != EXPECTED["brazil_strict_all"]:
        raise ValueError(f"Unexpected Brazil strict N20 count: {int(clinical_join.strict_N20.sum())}")
    if int(clinical_join.strict_N20_urgent.sum()) != EXPECTED["brazil_strict_urgent"]:
        raise ValueError("Unexpected urgent strict N20 count in Brazil aggregate")
    if int(clinical_join.strict_N20_elective.sum()) != EXPECTED["brazil_strict_elective"]:
        raise ValueError("Unexpected elective strict N20 count in Brazil aggregate")
    if not clinical_join.strict_N20.eq(
        clinical_join.strict_N20_urgent + clinical_join.strict_N20_elective
    ).all():
        raise ValueError("Strict Brazil N20 total does not equal urgent plus elective counts")

    panel = panel.merge(
        clinical_join[["hrid", "date", "strict_N20"]],
        on=["hrid", "date"], how="left", validate="one_to_one",
    )
    panel["strict_N20"] = panel.strict_N20.fillna(0).astype("int64")
    panel["w"] = panel["pop"].fillna(0)

    def aggregate_state(g):
        weight_sum = float(g.w.sum())
        mean_temp = (
            float(np.average(g.tmean.to_numpy(), weights=g.w.to_numpy()))
            if weight_sum > 0 else float(g.tmean.mean())
        )
        return pd.Series(
            {
                "broad_N20_N23": float(g.broad_N20_N23.sum()),
                "strict_N20": float(g.strict_N20.sum()),
                "tmean": mean_temp,
                "population_weight_sum": weight_sum,
            }
        )

    state_daily = panel.groupby(["state", "date"], as_index=False).apply(
        aggregate_state, include_groups=False
    )
    state_daily["state"] = state_daily.state.astype(str).str.zfill(2)
    state_daily["date"] = pd.to_datetime(state_daily.date).dt.normalize()
    for col in ["broad_N20_N23", "strict_N20"]:
        state_daily[col] = state_daily[col].round().astype("int64")
    if int(state_daily.broad_N20_N23.sum()) != EXPECTED["brazil_broad"]:
        raise ValueError("Brazil state broad counts changed during state aggregation")
    if int(state_daily.strict_N20.sum()) != EXPECTED["brazil_strict_all"]:
        raise ValueError("Brazil state strict N20 counts changed during state aggregation")
    assert_complete_state_calendar(state_daily, "Brazil", 27)
    state_daily["ym"] = state_daily.date.dt.strftime("%Y%m")
    state_daily["dow"] = state_daily.date.dt.dayofweek
    state_daily["stratum"] = state_daily.state + "_" + state_daily.ym + "_" + state_daily.dow.astype(str)
    add_state_counts(
        "Brazil", state_daily,
        [("broad_N20_N23", "broad N20-N23"), ("strict_N20", "strict N20 all admissions")],
    )
    SOURCE_COVERAGE["brazil"] = {
        "historical_panel": str(BRAZIL_PANEL_PATH),
        "clinical_aggregate": str(BRAZIL_CLINICAL_PATH),
        "clinical_aggregate_schema": parquet_columns(BRAZIL_CLINICAL_PATH),
        "clinical_source_contains_personal_id": False,
        "clinical_definition": "diagnosis='N20'; CAR_INT 01 and 02 included as all admission types",
        "window": "2018-01-01 through 2022-12-31",
        "health_region_input_rows": len(panel),
        "health_regions": int(panel.hrid.nunique()),
        "state_count": int(state_daily.state.nunique()),
        "state_day_rows": len(state_daily),
        "days_per_state": WINDOW_DAYS,
        "broad_N20_N23_total": int(state_daily.broad_N20_N23.sum()),
        "strict_N20_all_admissions": int(state_daily.strict_N20.sum()),
        "strict_N20_urgent": int(clinical_join.strict_N20_urgent.sum()),
        "strict_N20_elective": int(clinical_join.strict_N20_elective.sum()),
        "strict_N20_CAR_INT_levels": int(clinical_type_summary.n20_admission_type_levels),
        "broad_daily_match_health_region_days": int((~mismatch).sum()),
        "broad_daily_mismatch_health_region_days": int(mismatch.sum()),
        "temperature_definition": "same as fit_same_caliber.py: within-state population-weighted mean of health-region tmean; unweighted mean fallback only if total population weight is zero",
        "population_weight_fallback_state_days": int((state_daily.population_weight_sum <= 0).sum()),
    }
    return state_daily


def main():
    mexico, mexico_era5, mexico_nasa = load_mexico()
    brazil = load_brazil()

    # Broad reproduction precedes the strict endpoint for each exposure set.
    brazil_exposure = brazil[["state", "date", "tmean"]].rename(
        columns={"tmean": "temperature"}
    )
    attach_lags_and_fit(
        "Brazil", "population_weighted_panel_temperature", brazil_exposure, brazil
    )
    attach_lags_and_fit(
        "Mexico", "ERA5_Land_area_weighted_state_daily",
        mexico_era5, mexico,
    )
    attach_lags_and_fit(
        "Mexico", "NASA_POWER_state_capital_daily_supplement",
        mexico_nasa, mexico,
    )

    (OUT / "mexico_strict_source_coverage.json").write_text(
        json.dumps(SOURCE_COVERAGE, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    note = """# 墨西哥与巴西的严格 N20 表型核查

原墨西哥面板中的 `N20` 按 `extract_mexico_n20.py` 的主诊断前三位筛选，包含 N20、N21、N22、N23。宽口径40,415例和原面板全部58,432个州日逐日一致。无个人ID的既有临床程序抽取文件提供严格诊断 `AFECPRIN` 前缀 N20 日计数：32,417例。严格口径依照入院记录计数；没有将按ID/年份去重的支架/手术分组当成总入院结局。

墨西哥主温度产品为已覆盖2018–2022的ERA5-Land州日温度：小时2米温度先转为摄氏度并取日均值，再按州界内网格点的 `cos(latitude)` 权重计算面积加权均值。输入文件名标有2015–2023，实际逐日记录止于2022-12-31；本分析窗口内32州日均完整。NASA POWER的同协变量T2M在州府坐标处提取，作为补充产品单独报告。两种温度产品的空间代表性不同，不能按同一加权方法解释。

巴西严格N20从本地去标识日聚合中按 `diagnosis='N20'` 且 `CAR_INT` 为01或02构造，总数321,201例；其中急诊226,402例、择期94,799例。广义N20–N23日计数逐卫生区与原面板匹配。巴西州温度严格复现 `fit_same_caliber.py`：先对州内卫生区的日均温度按人口加权。

六个规格均在同一个2018–2022州日样本、状态分层和0–7日滞后模型中先重现宽口径，再估计严格N20；国家之间不合并RR。累积效应为温度滞后系数之和，标准误使用完整协方差。结果同时给出历史可比的正态近似区间和按州聚类数 `G-1` 自由度的t近似区间。州聚类数有限，t区间是敏感性近似，不改变主模型定义。

州级结局计数、输入与拟合样本数、零事件分层及模型公式见CSV；每个模型的完整系数协方差见 `mexico_strict_full_covariance.json`。所有临床输出均为州/日期聚合，不包含个人标识。
"""
    (OUT / "mexico_strict_phenotype_note_zh.md").write_text(note, encoding="utf-8")
    (OUT / "mexico_strict_manifest.json").write_text(
        json.dumps(
            {
                "script": str(Path(__file__).resolve()),
                "window": "2018-01-01 through 2022-12-31",
                "lag_days": "0-7; filter window and check complete daily calendar before lagging",
                "models": list(FULL_MODELS),
                "outputs": [
                    "mexico_strict_model_results.csv",
                    "mexico_strict_model_summary.csv",
                    "mexico_strict_endpoint_counts_by_state.csv",
                    "mexico_strict_source_coverage.json",
                    "mexico_strict_full_covariance.json",
                    "mexico_strict_phenotype_note_zh.md",
                ],
                "source_coverage": SOURCE_COVERAGE,
            },
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )
    persist()
    print(f"COMPLETE models={len(FULL_MODELS)} output={OUT}", flush=True)


if __name__ == "__main__":
    main()
