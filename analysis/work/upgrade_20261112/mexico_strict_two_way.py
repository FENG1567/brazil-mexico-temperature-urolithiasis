"""Run only the prespecified Brazil and Mexico strict-N20 two-way fits.

Reuses the verified aggregate loaders from mexico_strict_phenotype.py, while
writing the sensitivity results to separate artifacts. No row-level IDs are
loaded or printed.
"""
import json
import sys
import time
from pathlib import Path

# Match the dependency path used by the existing analysis scripts.
sys.path.insert(0, "F:/analysis/climate_biliary_20261102/work/pylibs")

import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import norm, t as student_t

import mexico_strict_phenotype as strict


OUT = strict.OUT
RESULTS_PATH = OUT / "mexico_strict_two_way_model_results.csv"
SUMMARY_PATH = OUT / "mexico_strict_two_way_model_summary.csv"
COVARIANCE_PATH = OUT / "mexico_strict_two_way_full_covariance.json"
MANIFEST_PATH = OUT / "mexico_strict_two_way_manifest.json"
EXPECTED = {
    "Brazil": {
        "states": 27,
        "model_n": 45_436,
        "RR_per_C": 1.0130309362413275,
    },
    "Mexico": {
        "states": 32,
        "model_n": 39_359,
        "RR_per_C": 1.0264757508632651,
    },
}


def state_temperature_frame(state_daily, exposure_frame):
    exposure = exposure_frame[["state", "date", "temperature"]].copy()
    exposure["state"] = exposure.state.astype(str).str.zfill(2)
    exposure["date"] = pd.to_datetime(exposure.date).dt.normalize()
    if exposure.duplicated(["state", "date"]).any():
        raise ValueError("Duplicate exposure state-days")

    frame = state_daily[
        ["state", "date", "stratum", "strict_N20"]
    ].copy()
    frame["state"] = frame.state.astype(str).str.zfill(2)
    frame["date"] = pd.to_datetime(frame.date).dt.normalize()
    if frame.duplicated(["state", "date"]).any():
        raise ValueError("Duplicate outcome state-days")
    if not frame.date.between("2018-01-01", "2022-12-31").all():
        raise ValueError("Outcome rows outside the prespecified 2018–2022 window")

    joined = frame.merge(
        exposure, on=["state", "date"], how="left", validate="one_to_one"
    )
    if joined.temperature.isna().any():
        raise ValueError("Missing temperature on matched state-days")
    joined = joined.sort_values(["state", "date"]).reset_index(drop=True)
    diffs = joined.groupby("state", observed=True).date.diff().dropna().dt.days
    if not diffs.eq(1).all():
        raise ValueError("Nonconsecutive state calendars; refusing to create lags")

    for lag in range(8):
        joined[f"t{lag}"] = joined.groupby("state", sort=False).temperature.shift(lag)
    joined["date_id"] = joined.date.dt.strftime("%Y%m%d")
    return joined


def fit_two_way(country, product, frame):
    outcome = "strict_N20"
    lag_terms = [f"t{lag}" for lag in range(8)]
    model_frame = frame[
        ["state", "date", "date_id", "stratum", *lag_terms, outcome]
    ].dropna(subset=lag_terms).copy()
    if model_frame.empty:
        raise ValueError(f"Empty two-way sample for {country}")

    expected = EXPECTED[country]
    if int(model_frame.state.nunique()) != expected["states"]:
        raise ValueError(f"Unexpected state-cluster count for {country}")
    if int(model_frame.date_id.nunique()) != 1_819:
        raise ValueError(f"Unexpected date-cluster count for {country}")

    model_key = f"{country.lower()}_{product}_strict_N20_lag07_two_way"
    formula = outcome + " ~ " + " + ".join(lag_terms) + " | stratum"
    strata_events = model_frame.groupby("stratum", observed=True)[outcome].sum()
    strata_sizes = model_frame.groupby("stratum", observed=True).size()
    zero_levels = strata_events.index[strata_events.eq(0)]
    zero_fe_rows = int(strata_sizes.reindex(zero_levels).sum())
    input_n = len(model_frame)
    input_events = int(model_frame[outcome].sum())
    state_clusters = int(model_frame.state.nunique())
    date_clusters = int(model_frame.date_id.nunique())

    print(
        f"START {model_key} input_n={input_n} input_events={input_events} "
        f"state_clusters={state_clusters} date_clusters={date_clusters}",
        flush=True,
    )
    started = time.time()
    model = pf.fepois(
        formula, data=model_frame, vcov={"CRV1": "state + date_id"}
    )
    coef = model.coef()
    names = [str(name) for name in coef.index]
    values = coef.to_numpy(dtype=float)
    covariance = np.asarray(model._vcov, dtype=float)
    if covariance.shape != (len(names), len(names)):
        raise ValueError(f"Unexpected covariance dimensions: {covariance.shape}")
    if not np.allclose(covariance, covariance.T, rtol=0, atol=1e-10):
        raise ValueError("Two-way covariance matrix is not symmetric")

    weights = np.array([1.0 if name in lag_terms else 0.0 for name in names])
    if not weights.any():
        raise ValueError(f"No lag coefficients found: {names}")
    beta = float(weights @ values)
    variance = float(weights @ covariance @ weights)
    if variance < 0 and abs(variance) < 1e-12:
        variance = 0.0
    if variance < 0:
        raise ValueError(f"Negative cumulative variance: {variance}")
    se = float(np.sqrt(variance))
    z = beta / se if se > 0 else float("inf")
    p_normal = float(2 * norm.sf(abs(z))) if np.isfinite(z) else 0.0

    df = min(state_clusters, date_clusters) - 1
    if df <= 0:
        raise ValueError("Insufficient two-way clusters")
    tcrit = float(student_t.ppf(0.975, df))
    p_t = float(2 * student_t.sf(abs(z), df)) if np.isfinite(z) else 0.0
    normal_ci = [float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se))]
    t_ci = [float(np.exp(beta - tcrit * se)), float(np.exp(beta + tcrit * se))]
    rr = float(np.exp(beta))
    model_n = int(model._N)
    dropped_n = input_n - model_n
    expected_rr = expected["RR_per_C"]
    if model_n != expected["model_n"]:
        raise ValueError(
            f"{country} fitted N changed: expected {expected['model_n']}, got {model_n}"
        )
    if not np.isclose(rr, expected_rr, rtol=0, atol=1e-8):
        raise ValueError(
            f"{country} two-way point estimate does not match the saved one-way fit: "
            f"{rr} versus {expected_rr}"
        )

    summary = {
        "model": model_key,
        "country": country,
        "exposure_product": product,
        "phenotype": "strict_N20",
        "diagnosis_definition": (
            "all admissions with principal diagnosis prefix N20"
            if country == "Mexico"
            else "all admissions with diagnosis N20; CAR_INT 01 and 02"
        ),
        "formula": formula,
        "fixed_effects": "state x calendar year-month x day-of-week stratum",
        "cluster": "CRV1 by state + date_id (YYYYMMDD)",
        "state_clusters": state_clusters,
        "date_clusters": date_clusters,
        "t_df_min_clusters_minus_1": df,
        "input_n": input_n,
        "model_n": model_n,
        "rows_dropped": dropped_n,
        "input_strata": int(model_frame.stratum.nunique()),
        "zero_fe_strata": int(len(zero_levels)),
        "zero_fe_rows_expected_dropped": zero_fe_rows,
        "other_rows_dropped_by_fit": dropped_n - zero_fe_rows,
        "events_full_window": int(frame[outcome].sum()),
        "events_after_lag_window": input_events,
        "model_events": input_events if dropped_n == zero_fe_rows else None,
        "one_way_reference_RR_per_C": expected_rr,
        "point_estimate_matches_one_way": True,
        "seconds": round(time.time() - started, 2),
    }
    result = {
        **summary,
        "contrast": "cumulative_temperature_lag_0_7_per_C",
        "beta_sum": beta,
        "se_full_covariance": se,
        "RR_per_C": rr,
        "normal_CI_low": normal_ci[0],
        "normal_CI_high": normal_ci[1],
        "normal_p": p_normal,
        "t_min_clusters_minus_1_CI_low": t_ci[0],
        "t_min_clusters_minus_1_CI_high": t_ci[1],
        "t_min_clusters_minus_1_p": p_t,
        "contrast_terms": json.dumps(lag_terms),
    }
    full_model = {
        **summary,
        "coefficient_names": names,
        "coefficients": values.tolist(),
        "covariance": covariance.tolist(),
        "cumulative_temperature_contrast": {
            "weights": weights.tolist(),
            "beta_sum": beta,
            "variance": variance,
            "se": se,
            "RR_per_C": rr,
            "normal_CI": normal_ci,
            "normal_p": p_normal,
            "t_min_clusters_minus_1_CI": t_ci,
            "t_min_clusters_minus_1_p": p_t,
        },
    }
    print(
        f"RESULT {model_key} RR={rr:.9f} normalCI={normal_ci} "
        f"tCI={t_ci} df={df} n={model_n} "
        f"state_clusters={state_clusters} date_clusters={date_clusters}",
        flush=True,
    )
    return result, summary, full_model


def main():
    target_paths = [RESULTS_PATH, SUMMARY_PATH, COVARIANCE_PATH, MANIFEST_PATH]
    existing = [path.name for path in target_paths if path.exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite two-way artifacts: {existing}")

    # Reuse the already-audited source loaders. They do not fit models.
    mexico, mexico_era5, _mexico_nasa = strict.load_mexico()
    brazil = strict.load_brazil()

    brazil_exposure = brazil[["state", "date", "tmean"]].rename(
        columns={"tmean": "temperature"}
    )
    brazil_frame = state_temperature_frame(brazil, brazil_exposure)
    mexico_frame = state_temperature_frame(mexico, mexico_era5)

    results = []
    summaries = []
    models = {}
    for country, product, frame in [
        ("Brazil", "population_weighted_panel_temperature", brazil_frame),
        ("Mexico", "ERA5_Land_area_weighted_state_daily", mexico_frame),
    ]:
        result, summary, full_model = fit_two_way(country, product, frame)
        results.append(result)
        summaries.append(summary)
        models[summary["model"]] = full_model

    pd.DataFrame(results).to_csv(RESULTS_PATH, index=False, encoding="utf-8-sig")
    pd.DataFrame(summaries).to_csv(SUMMARY_PATH, index=False, encoding="utf-8-sig")
    COVARIANCE_PATH.write_text(
        json.dumps(
            {
                "window": "2018-01-01 through 2022-12-31; filtered before lags",
                "contrast": "sum of temperature coefficients for lags 0-7",
                "cluster": "two-way CRV1 by state and YYYYMMDD date_id",
                "models": models,
            },
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "script": str(Path(__file__).resolve()),
                "source_script": str(Path(strict.__file__).resolve()),
                "window": "2018-01-01 through 2022-12-31",
                "lag_days": "0-7; sample window is filtered before lags",
                "date_cluster_id": "YYYYMMDD",
                "models": [item["model"] for item in results],
                "outputs": [path.name for path in target_paths],
                "point_estimates_match_one_way": all(
                    item["point_estimate_matches_one_way"] for item in results
                ),
            },
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ),
        encoding="utf-8",
    )
    print(f"COMPLETE models={len(models)} output={OUT}", flush=True)


if __name__ == "__main__":
    main()
