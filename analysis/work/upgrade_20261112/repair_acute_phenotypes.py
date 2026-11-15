"""Reproduce the broad N20-N23 acute phenotype and run strict N20 repair fits.

This script reads the frozen source exposure panel and the de-identified daily
clinical aggregate. It writes only acute_repair_* artifacts under outputs.
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
from scipy.stats import norm


BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "work/upgrade_20261112"
OUT = BASE / "outputs/upgrade_20261112"
SOURCE_WORK = Path("F:/analysis/climate_biliary_20261102/work")
PANEL_PATH = SOURCE_WORK / "panel_composite_humidity.parquet"
AGGREGATE_PATH = WORK / "clinical_code_day_aggregate.parquet"
OUT.mkdir(parents=True, exist_ok=True)

GROUP_LEVELS = ["neither", "heat_only", "smoke_only", "both"]
DIAGNOSES = ("N20", "N21", "N22", "N23")
DRAINAGE_CODES = [
    "0409010103", "0409010162", "0409010170",
    "0409010278", "0409010286", "0409010294",
]
DOUBLE_J_CODE = "0409010170"
EXPECTED_DRAINAGE_EVENTS = 47_918
EXPECTED_DOUBLE_J_EVENTS = 46_949

db = duckdb.connect()
db.execute("SET threads=8")

RESULTS = []
SUMMARIES = []
COUNTS = []
FULL_MODELS = {}


def sql_path(path):
    return Path(path).as_posix().replace("'", "''")


def persist():
    pd.DataFrame(RESULTS).to_csv(
        OUT / "acute_repair_model_results.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(SUMMARIES).to_csv(
        OUT / "acute_repair_model_summary.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(COUNTS).to_csv(
        OUT / "acute_repair_endpoint_counts.csv", index=False, encoding="utf-8-sig"
    )
    payload = {
        "panel_source": str(PANEL_PATH),
        "clinical_aggregate": str(AGGREGATE_PATH),
        "window": "2018-01-01 through 2022-12-31; window filtered before lagging",
        "cluster": "CRV1 by hrid",
        "models": FULL_MODELS,
    }
    (OUT / "acute_repair_full_covariance.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )


def fit_model(key, formula, frame, zero_fe_columns, count_column, tags):
    """Fit one CRV1 Poisson model and retain all coefficients and covariance."""
    started = time.time()
    input_n = len(frame)
    fe_sum = frame.groupby(zero_fe_columns, observed=True)[count_column].sum()
    fe_size = frame.groupby(zero_fe_columns, observed=True).size()
    zero_levels = fe_sum.index[fe_sum.eq(0)]
    zero_fe_rows = int(fe_size.reindex(zero_levels).sum())
    zero_fe_cells = int(len(zero_levels))
    print(
        f"START {key} input_n={input_n:,} zero_fe_cells={zero_fe_cells:,} "
        f"formula={formula}",
        flush=True,
    )
    model = pf.fepois(formula, data=frame, vcov={"CRV1": "hrid"})
    coef = model.coef()
    names = [str(name) for name in coef.index]
    values = coef.to_numpy(dtype=float)
    covariance = np.asarray(model._vcov, dtype=float)
    if covariance.shape != (len(names), len(names)):
        raise ValueError(f"Unexpected covariance dimensions for {key}: {covariance.shape}")
    model_n = int(model._N)
    events = float(frame[count_column].sum())
    fitted_events = (
        float(events) if input_n - model_n == zero_fe_rows else None
    )  # all-zero FE cells carry no events; avoid guessing if anything else is dropped
    summary = {
        "model": key,
        **tags,
        "formula": formula,
        "cluster": "CRV1:hrid",
        "input_n": input_n,
        "model_n": model_n,
        "dropped_n": input_n - model_n,
        "input_events": events,
        "model_events": fitted_events,
        "zero_fe_cells_removed": zero_fe_cells,
        "zero_fe_rows_expected_removed": zero_fe_rows,
        "other_rows_removed_by_fit": input_n - model_n - zero_fe_rows,
        "input_strata": int(frame["stratum"].nunique()),
        "seconds": round(time.time() - started, 2),
    }
    SUMMARIES.append(summary)
    FULL_MODELS[key] = {
        **summary,
        "coefficient_names": names,
        "coefficients": values.tolist(),
        "covariance": covariance.tolist(),
        "contrasts": {},
    }
    print(
        f"FITTED {key} model_n={model_n:,} zero_fe_cells={zero_fe_cells:,} "
        f"dropped_n={input_n - model_n:,} seconds={summary['seconds']}",
        flush=True,
    )
    return model, names, values, covariance, summary


def add_contrast(key, group, names, values, covariance, weights, contrast, summary,
                 measure, label):
    beta = float(weights @ values)
    variance = float(weights @ covariance @ weights)
    if variance < 0 and abs(variance) < 1e-12:
        variance = 0.0
    if variance < 0:
        raise ValueError(f"Negative contrast variance for {key}/{contrast}: {variance}")
    se = float(np.sqrt(variance))
    z = beta / se if se > 0 else float("inf")
    p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else 0.0
    row = {
        "model": key,
        "phenotype": label,
        "contrast": contrast,
        "measure": measure,
        "RR": float(np.exp(beta)),
        "CI_low": float(np.exp(beta - 1.96 * se)),
        "CI_high": float(np.exp(beta + 1.96 * se)),
        "p": p,
        "beta": beta,
        "se": se,
        "variance": variance,
        "contrast_terms": json.dumps(
            [name for name, weight in zip(names, weights) if weight != 0],
            ensure_ascii=False,
        ),
        **summary,
    }
    RESULTS.append(row)
    FULL_MODELS[key]["contrasts"][contrast] = {
        "weights": weights.tolist(),
        "beta": beta,
        "variance": variance,
        "se": se,
        "RR": row["RR"],
        "CI_low": row["CI_low"],
        "CI_high": row["CI_high"],
        "p": p,
    }
    print(
        f"RESULT {key} {contrast} RR={row['RR']:.9f} "
        f"CI=[{row['CI_low']:.9f},{row['CI_high']:.9f}] p={p:.9g}",
        flush=True,
    )


def fit_separate_categorical(frame, phenotype, outcome, label):
    key = f"{phenotype}_{outcome}_categorical"
    formula = f"{outcome} ~ C(group) | stratum"
    fitted = frame[["hrid", "stratum", "group", outcome]].copy()
    model, names, values, covariance, summary = fit_model(
        key, formula, fitted, ["stratum"], outcome,
        {"analysis": "separate_categorical", "outcome": outcome},
    )
    for group in GROUP_LEVELS[1:]:
        term = f"C(group)[T.{group}]"
        if term not in names:
            raise ValueError(f"Categorical coefficient missing: {term} in {key}: {names}")
        weights = np.zeros(len(names))
        weights[names.index(term)] = 1.0
        add_contrast(
            key, group, names, values, covariance, weights,
            f"{group}_vs_neither", summary, "rate_ratio_vs_neither", label,
        )
    persist()


def fit_joint_categorical(frame, phenotype, label):
    key = f"{phenotype}_urgent_elective_joint_categorical"
    long = frame[["hrid", "stratum", "group", "urgent", "elective"]].melt(
        id_vars=["hrid", "stratum", "group"],
        value_vars=["elective", "urgent"],
        var_name="outcome",
        value_name="count",
    )
    long["outcome"] = pd.Categorical(long.outcome, categories=["elective", "urgent"])
    formula = "count ~ C(outcome)*C(group) | stratum^outcome"
    model, names, values, covariance, summary = fit_model(
        key, formula, long, ["stratum", "outcome"], "count",
        {"analysis": "joint_categorical", "outcome": "urgent_vs_elective_effect_ratio"},
    )
    for group in GROUP_LEVELS[1:]:
        candidates = [
            i for i, name in enumerate(names)
            if "C(outcome)[T.urgent]" in name
            and f"C(group)[T.{group}]" in name
            and ":" in name
        ]
        if len(candidates) != 1:
            raise ValueError(f"Expected one joint interaction for {group}: {names}")
        weights = np.zeros(len(names))
        weights[candidates[0]] = 1.0
        add_contrast(
            key, group, names, values, covariance, weights,
            f"{group}_urgent_RR_over_elective_RR", summary,
            "urgent_RR_over_elective_RR", label,
        )
    persist()


def fit_lag_temperature(frame, phenotype, outcome_values, label, max_lag):
    key = f"{phenotype}_urgent_temperature_lag_0_{max_lag}"
    terms = [f"t{lag}" for lag in range(max_lag + 1)]
    model_frame = frame[["hrid", "stratum", *terms]].copy()
    model_frame["urgent"] = outcome_values.to_numpy()
    model_frame = model_frame.dropna(subset=terms).copy()
    formula = "urgent ~ " + " + ".join(terms) + " | stratum"
    model, names, values, covariance, summary = fit_model(
        key, formula, model_frame, ["stratum"], "urgent",
        {
            "analysis": "continuous_temperature_lag",
            "outcome": "urgent",
            "max_lag_days": max_lag,
        },
    )
    weights = np.array([1.0 if name in terms else 0.0 for name in names])
    if not weights.any():
        raise ValueError(f"No temperature lag coefficients in {key}: {names}")
    add_contrast(
        key, None, names, values, covariance, weights,
        f"cumulative_temperature_lag_0_{max_lag}_per_C", summary,
        "cumulative_RR_per_C", label,
    )
    persist()


def fit_code_sensitivity(frame, outcome, name, label, include_same_day_rh):
    key = name
    lag_terms = [f"t{lag}" for lag in range(8)] + [f"f{lag}" for lag in range(8)]
    terms = lag_terms + (["rh"] if include_same_day_rh else [])
    model_frame = frame[["hrid", "stratum", *terms]].copy()
    model_frame["count"] = outcome.to_numpy()
    model_frame = model_frame.dropna(subset=terms).copy()
    formula = "count ~ " + " + ".join(terms) + " | stratum"
    model, names, values, covariance, summary = fit_model(
        key, formula, model_frame, ["stratum"], "count",
        {
            "analysis": "prespecified_code_sensitivity",
            "outcome": label,
            "temperature_lag_days": "0-7",
            "fire_exposure": "fire_pm25 / 10, lags 0-7",
            "same_day_rh_adjusted": bool(include_same_day_rh),
        },
    )
    heat_terms = [f"t{lag}" for lag in range(8)]
    weights = np.array([1.0 if term in heat_terms else 0.0 for term in names])
    if not weights.any():
        raise ValueError(f"No temperature lag coefficients in {key}: {names}")
    add_contrast(
        key, None, names, values, covariance, weights,
        "cumulative_temperature_lag_0_7_per_C", summary,
        "cumulative_RR_per_C", label,
    )
    persist()


def main():
    panel = db.sql(
        "SELECT * FROM read_parquet('" + sql_path(PANEL_PATH) + "') "
        "WHERE YEAR(date) BETWEEN 2018 AND 2022 AND fire_pm25 IS NOT NULL"
    ).df()
    aggregate_path = sql_path(AGGREGATE_PATH)
    meta = db.sql(
        f"""
        SELECT COUNT(*) AS rows,
               COUNT(DISTINCT IDENT) AS ident_levels,
               COUNT(DISTINCT diagnosis) AS diagnosis_levels,
               COUNT(DISTINCT CAR_INT) AS admission_type_levels,
               SUM(admissions) AS all_admissions,
               MIN(date) AS min_date, MAX(date) AS max_date
        FROM read_parquet('{aggregate_path}')
        WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
        """
    ).df().iloc[0].to_dict()
    daily = db.sql(
        f"""
        SELECT hrid, date,
          SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') THEN admissions ELSE 0 END) AS broad_total,
          SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') AND CAR_INT='02' THEN admissions ELSE 0 END) AS broad_urgent,
          SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') AND CAR_INT='01' THEN admissions ELSE 0 END) AS broad_elective,
          SUM(CASE WHEN diagnosis='N20' THEN admissions ELSE 0 END) AS strict_total,
          SUM(CASE WHEN diagnosis='N20' AND CAR_INT='02' THEN admissions ELSE 0 END) AS strict_urgent,
          SUM(CASE WHEN diagnosis='N20' AND CAR_INT='01' THEN admissions ELSE 0 END) AS strict_elective,
          SUM(CASE WHEN diagnosis='N20' AND CAR_INT='02'
                    AND PROC_REA IN ({','.join(repr(c) for c in DRAINAGE_CODES)})
                   THEN admissions ELSE 0 END) AS drainage,
          SUM(CASE WHEN diagnosis='N20' AND CAR_INT='02'
                    AND PROC_REA='{DOUBLE_J_CODE}'
                   THEN admissions ELSE 0 END) AS double_j_only
        FROM read_parquet('{aggregate_path}')
        WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
        GROUP BY hrid, date
        """
    ).df()

    panel["date"] = pd.to_datetime(panel.date)
    panel["hrid"] = panel.hrid.astype(str)
    daily["date"] = pd.to_datetime(daily.date)
    daily["hrid"] = daily.hrid.astype(str)
    if panel.duplicated(["hrid", "date"]).any():
        raise ValueError("Source panel has duplicate health-region dates")
    if daily.duplicated(["hrid", "date"]).any():
        raise ValueError("Clinical daily aggregate has duplicate health-region dates")
    if not set(daily.hrid).issubset(set(panel.hrid)):
        # The clinical aggregate is sparse by design; every event region must map.
        raise ValueError("Unexpected health-region mapping discrepancy")
    if not panel.fire_pm25.notna().all():
        raise ValueError("Categorical input unexpectedly has missing fire exposure")
    panel = panel.merge(
        daily, on=["hrid", "date"], how="left", validate="one_to_one"
    )
    endpoint_columns = [
        "broad_total", "broad_urgent", "broad_elective",
        "strict_total", "strict_urgent", "strict_elective",
        "drainage", "double_j_only",
    ]
    panel[endpoint_columns] = panel[endpoint_columns].fillna(0).astype("int64")
    panel["N20_N23"] = panel.N20_N23.astype("int64")
    mismatch = panel.broad_total.ne(panel.N20_N23)
    if mismatch.any():
        raise ValueError(
            f"Broad N20-N23 aggregate mismatches panel on {int(mismatch.sum())} days"
        )
    if set(panel.heat.dropna().unique()) - {0, 1} or set(panel.smoke.dropna().unique()) - {0, 1}:
        raise ValueError("Heat/smoke exposure fields are not binary")
    if panel[["heat", "smoke", "stratum", "rh", "tmean"]].isna().any().any():
        raise ValueError("Required panel fields contain missing values")

    panel["group"] = np.select(
        [
            (panel.heat.eq(0) & panel.smoke.eq(0)),
            (panel.heat.eq(1) & panel.smoke.eq(0)),
            (panel.heat.eq(0) & panel.smoke.eq(1)),
            (panel.heat.eq(1) & panel.smoke.eq(1)),
        ],
        GROUP_LEVELS,
        default="other",
    )
    if not panel.group.isin(GROUP_LEVELS).all():
        raise ValueError("At least one exposure day is outside the four intended groups")
    panel["group"] = pd.Categorical(panel.group, categories=GROUP_LEVELS)
    panel = panel.sort_values(["hrid", "date"]).reset_index(drop=True)
    if not panel.groupby("hrid", observed=True).date.diff().dropna().dt.days.eq(1).all():
        raise ValueError("Calendar has gaps; refusing to create chronological lags")

    # Confirm source-code scope and source counts before fitting any association.
    idents = db.sql(
        f"SELECT DISTINCT IDENT FROM read_parquet('{aggregate_path}') ORDER BY IDENT"
    ).df().IDENT.astype(str).tolist()
    admission_types = db.sql(
        f"SELECT DISTINCT CAR_INT FROM read_parquet('{aggregate_path}') ORDER BY CAR_INT"
    ).df().CAR_INT.astype(str).tolist()
    diagnoses = db.sql(
        f"SELECT DISTINCT diagnosis FROM read_parquet('{aggregate_path}') ORDER BY diagnosis"
    ).df().diagnosis.astype(str).tolist()
    if idents != ["1"]:
        raise ValueError(f"Expected aggregate IDENT='1' only; found {idents}")
    if not set(diagnoses).issubset(DIAGNOSES) or not set(admission_types).issubset({"01", "02"}):
        raise ValueError(f"Unexpected clinical codes: diagnosis={diagnoses}, CAR_INT={admission_types}")
    broad_urgent_events = int(panel.broad_urgent.sum())
    broad_elective_events = int(panel.broad_elective.sum())
    strict_urgent_events = int(panel.strict_urgent.sum())
    strict_elective_events = int(panel.strict_elective.sum())
    drainage_events = int(panel.drainage.sum())
    double_j_events = int(panel.double_j_only.sum())
    if drainage_events != EXPECTED_DRAINAGE_EVENTS:
        raise ValueError(
            f"Drainage whitelist event count {drainage_events} != expected {EXPECTED_DRAINAGE_EVENTS}"
        )
    if double_j_events != EXPECTED_DOUBLE_J_EVENTS:
        raise ValueError(
            f"Double-J event count {double_j_events} != expected {EXPECTED_DOUBLE_J_EVENTS}"
        )

    for endpoint, column in [
        ("broad_N20_N23_total", "broad_total"),
        ("broad_N20_N23_urgent", "broad_urgent"),
        ("broad_N20_N23_elective", "broad_elective"),
        ("strict_N20_total", "strict_total"),
        ("strict_N20_urgent", "strict_urgent"),
        ("strict_N20_elective", "strict_elective"),
        ("primary_drainage_whitelist", "drainage"),
        ("double_j_only_0409010170", "double_j_only"),
    ]:
        by_group = panel.groupby("group", observed=True)[column].sum().reindex(GROUP_LEVELS, fill_value=0)
        days = panel.groupby("group", observed=True).size().reindex(GROUP_LEVELS, fill_value=0)
        for group in GROUP_LEVELS:
            COUNTS.append({
                "endpoint": endpoint,
                "exposure_group": group,
                "input_days": int(days[group]),
                "admissions": int(by_group[group]),
            })
    persist()

    print(
        "INPUTS " + json.dumps({
            "panel_rows": len(panel),
            "aggregate": {key: (str(value) if isinstance(value, (pd.Timestamp,)) else value) for key, value in meta.items()},
            "broad_urgent": broad_urgent_events,
            "broad_elective": broad_elective_events,
            "strict_N20_urgent": strict_urgent_events,
            "strict_N20_elective": strict_elective_events,
            "drainage_whitelist": drainage_events,
            "double_j_only": double_j_events,
            "panel_rows_matching_broad_total": int((~mismatch).sum()),
        }, ensure_ascii=False, default=str),
        flush=True,
    )

    # Reproduce the old categorical comparisons on the original stratum/cluster.
    for phenotype, label, total, urgent, elective in [
        ("broad_N20_N23", "broad N20-N23", "broad_total", "broad_urgent", "broad_elective"),
        ("strict_N20", "strict N20", "strict_total", "strict_urgent", "strict_elective"),
    ]:
        for outcome in [total, urgent, elective]:
            fit_separate_categorical(panel, phenotype, outcome, label)
        fit_joint_categorical(
            panel.rename(columns={urgent: "urgent", elective: "elective"}),
            phenotype,
            label,
        )

    # The historical temperature models filter to 2018-2022 before creating
    # region-specific chronological lags. No pre-window days are carried in.
    panel["fire10"] = panel.fire_pm25 / 10.0
    for lag in range(22):
        panel[f"t{lag}"] = panel.groupby("hrid", sort=False).tmean.shift(lag)
        panel[f"f{lag}"] = panel.groupby("hrid", sort=False).fire10.shift(lag)
    for max_lag in [7, 14, 21]:
        fit_lag_temperature(
            panel, "broad_N20_N23", panel.broad_urgent,
            "broad N20-N23 urgent", max_lag,
        )
        fit_lag_temperature(
            panel, "strict_N20", panel.strict_urgent,
            "strict N20 urgent", max_lag,
        )

    # Prospective addendum sensitivities for the primary procedure endpoint:
    # first add same-day RH to the six-code endpoint; then narrow to the double-J
    # code alone and retain the original lagged temperature/fire specification.
    fit_code_sensitivity(
        panel, panel.drainage,
        "sensitivity_drainage_whitelist_lag07_rh_same_day",
        "urgent N20 with drainage whitelist",
        include_same_day_rh=True,
    )
    fit_code_sensitivity(
        panel, panel.double_j_only,
        "sensitivity_double_j_0409010170_lag07_no_rh",
        "urgent N20 with procedure code 0409010170 only",
        include_same_day_rh=False,
    )

    note = """# 急诊结局诊断标签修复与复现说明

冻结稿第74行的“urgent N20”标签与原构建口径不一致。只读源脚本 `build_acute_n20.py` 以主诊断前三位筛选 N20、N21、N22、N23，因此原分类结局和温度滞后结局应称为“急诊 N20–N23”宽口径；本次另按 `diagnosis == N20` 估计“急诊 N20”严格口径。冻结文件和 DOCX 未修改。

源日面板为2018–2022年821,700个卫生区日，按原分层 `stratum` 拟合条件Poisson并按卫生区做 CRV1 聚类。宽口径包括急诊和择期分别模型及堆叠联合模型；联合模型中热/烟组的系数表示该组相对“皆无”组的急诊/择期率比之比。温度滞后先限制到2018–2022年窗口，再按卫生区、日期排序生成滞后；累积效应为0至指定日滞后温度系数之和，标准误使用完整协方差矩阵。

两项预先确定的临床代码敏感性均报告：六代码引流/支架组加同日相对湿度；以及仅限0409010170代码、沿用温度和火源PM2.5滞后0–7模型、不加相对湿度。它们检验预先指定的协变量和代码定义变化，不用于根据显著性重新定义主终点。所有计数来自无个人标识的聚合日数据。

模型输入行数、拟合行数、零事件固定效应层数、事件数及完整系数协方差见 `acute_repair_model_summary.csv` 和 `acute_repair_full_covariance.json`。分类组日数与住院数见 `acute_repair_endpoint_counts.csv`。
"""
    (OUT / "acute_repair_phenotype_note_zh.md").write_text(note, encoding="utf-8")
    persist()
    manifest = {
        "script": str(Path(__file__).resolve()),
        "panel_source": str(PANEL_PATH),
        "clinical_aggregate": str(AGGREGATE_PATH),
        "window": "2018-01-01 through 2022-12-31",
        "ident_levels": idents,
        "diagnosis_levels": diagnoses,
        "admission_type_levels": admission_types,
        "panel_input_rows": len(panel),
        "broad_urgent_admissions": broad_urgent_events,
        "broad_elective_admissions": broad_elective_events,
        "strict_N20_urgent_admissions": strict_urgent_events,
        "strict_N20_elective_admissions": strict_elective_events,
        "drainage_whitelist_codes": DRAINAGE_CODES,
        "drainage_whitelist_admissions": drainage_events,
        "double_j_code": DOUBLE_J_CODE,
        "double_j_admissions": double_j_events,
        "models": list(FULL_MODELS),
    }
    (OUT / "acute_repair_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"COMPLETE models={len(FULL_MODELS)} output={OUT}", flush=True)


if __name__ == "__main__":
    main()
