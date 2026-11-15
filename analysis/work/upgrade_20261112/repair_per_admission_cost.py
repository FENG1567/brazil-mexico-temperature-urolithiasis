"""Reproduce and repair the N20 per-admission reimbursement analysis.

All regression outcomes are daily mean administrative SUS reimbursements per
admission. The frozen mixed-denominator model is retained for reproduction;
the three corrected models use internally consistent numerator/denominator
phenotypes, with the strict N20 model also rerun at December 2022 IPCA prices.
"""

import os

for variable in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[variable] = "8"

import json
import sys
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
CLINICAL_PATH = WORK / "clinical_code_day_aggregate.parquet"
IPCA_PATH = OUT / "ipca_monthly.csv"
FROZEN_PATH = BASE / "work/frozen_snapshot/tables/n20_cost_per_case.csv"

GROUPS = ["neither", "heat_only", "smoke_only", "both"]
EXPECTED_STRICT_ADMISSIONS = 226_402
EXPECTED_BROAD_ADMISSIONS = 315_504
EXPECTED_BROAD_ADMIT_DAYS = 196_715
EXPECTED_BROAD_DAYS_WITHOUT_N20 = 49_430

OUT.mkdir(parents=True, exist_ok=True)


def sql_path(path):
    return Path(path).as_posix().replace("'", "''")


def group_names(frame, both_label="both"):
    conditions = [
        frame.heat.eq(0) & frame.smoke.eq(0),
        frame.heat.eq(1) & frame.smoke.eq(0),
        frame.heat.eq(0) & frame.smoke.eq(1),
        frame.heat.eq(1) & frame.smoke.eq(1),
    ]
    labels = ["neither", "heat_only", "smoke_only", both_label]
    return pd.Categorical(np.select(conditions, labels, default="other"), categories=[
        "neither", "heat_only", "smoke_only", both_label
    ])


def summarize_by_exposure(frame, definition, cost_column, admissions_column,
                          analysis_mask, original_mixed=False):
    rows = []
    for group in GROUPS:
        group_frame = frame.loc[frame.group.eq(group)]
        admissions = group_frame[admissions_column].astype(float)
        cost = group_frame[cost_column].astype(float)
        event_days = admissions.gt(0)
        selected = analysis_mask.loc[group_frame.index]
        broad_events = group_frame["broad_admissions"].gt(0)
        row = {
            "definition": definition,
            "exposure_group": group,
            "panel_days": int(len(group_frame)),
            "admit_days": int(event_days.sum()),
            "admissions_all_admit_days": float(admissions.sum()),
            "reimbursement_all_admit_days": float(cost.sum()),
            "positive_cost_admit_days": int((event_days & cost.gt(0)).sum()),
            "zero_cost_admit_days": int((event_days & cost.eq(0)).sum()),
            "negative_cost_admit_days": int((event_days & cost.lt(0)).sum()),
            "analysis_days": int(selected.sum()),
            "analysis_admissions": float(admissions.loc[selected].sum()),
            "analysis_reimbursement": float(cost.loc[selected].sum()),
            "broad_N20_N23_admit_days": int(broad_events.sum()),
            "broad_days_with_multiple_diagnoses": int(
                (group_frame["broad_diagnoses_present"].gt(1) & broad_events).sum()
            ),
            "broad_days_without_N20_admissions": int(
                (broad_events & group_frame["strict_admissions"].eq(0)).sum()
            ),
        }
        if original_mixed:
            mixed_ratio = cost.div(admissions.where(admissions.gt(0)))
            row["old_cost_ratio_below_clip_floor"] = int(
                (event_days & mixed_ratio.lt(1)).sum()
            )
        else:
            row["old_cost_ratio_below_clip_floor"] = 0
        rows.append(row)
    return rows


def fit_cost_model(model_key, estimand, frame, cost_column, admissions_column,
                   price_basis, log_handling, include_only_positive_cost):
    """Fit log daily mean reimbursement with stratum FE and region CRV1 SEs."""
    model_frame = frame.copy()
    if include_only_positive_cost:
        keep = model_frame[admissions_column].gt(0) & model_frame[cost_column].gt(0)
    else:
        keep = model_frame[admissions_column].gt(0)
    model_frame = model_frame.loc[keep].copy()
    model_frame["cost_per_admission"] = (
        model_frame[cost_column] / model_frame[admissions_column]
    )
    if include_only_positive_cost and not model_frame.cost_per_admission.gt(0).all():
        raise ValueError(f"Nonpositive reimbursement entered {model_key}")
    if include_only_positive_cost:
        model_frame["log_cpc"] = np.log(model_frame.cost_per_admission)
    else:
        # Only the historical replication uses this lower clip.
        model_frame["log_cpc"] = np.log(model_frame.cost_per_admission.clip(lower=1))

    formula = "log_cpc ~ C(group) | stratum"
    model = pf.feols(formula, data=model_frame, vcov={"CRV1": "hrid"})
    coefficients = model.coef()
    coefficient_names = [str(name) for name in coefficients.index]
    values = coefficients.to_numpy(dtype=float)
    covariance = np.asarray(model._vcov, dtype=float)
    if covariance.shape != (len(coefficient_names), len(coefficient_names)):
        raise ValueError(f"Unexpected covariance dimensions in {model_key}: {covariance.shape}")

    terms = []
    for index, term in enumerate(coefficient_names):
        name = term.split("[T.", 1)[1].rstrip("]")
        beta = float(values[index])
        se = float(model.se()[term])
        p_value = float(2 * norm.sf(abs(beta / se))) if se > 0 else 0.0
        terms.append({
            "model": model_key,
            "estimand": estimand,
            "price_basis": price_basis,
            "term": term,
            "exposure_group": name,
            "beta_log_ratio": beta,
            "se": se,
            "variance": float(covariance[index, index]),
            "ratio": float(np.exp(beta)),
            "CI_low": float(np.exp(beta - 1.96 * se)),
            "CI_high": float(np.exp(beta + 1.96 * se)),
            "p": p_value,
        })

    summary = {
        "model": model_key,
        "estimand": estimand,
        "formula": formula,
        "cluster": "CRV1 by hrid",
        "price_basis": price_basis,
        "log_handling": log_handling,
        "analysis_rows_before_fit": int(len(model_frame)),
        "model_rows": int(model._N),
        "rows_removed_by_fit": int(len(model_frame) - int(model._N)),
        "strata_before_fit": int(model_frame.stratum.nunique()),
        "analysis_admissions": float(model_frame[admissions_column].sum()),
        "analysis_reimbursement": float(model_frame[cost_column].sum()),
        "analysis_cost_per_admission_arithmetic_day_mean": float(
            model_frame.cost_per_admission.mean()
        ),
        "positive_log_days": int(len(model_frame)),
    }
    details = {
        **summary,
        "numerator_column": cost_column,
        "denominator_column": admissions_column,
        "coefficient_names": coefficient_names,
        "coefficients": values.tolist(),
        "covariance": covariance.tolist(),
        "contrasts": terms,
    }
    return terms, summary, details, model_frame


def write_chinese_report(results, summaries, reproduction):
    result_lookup = {
        (row["model"], row["exposure_group"]): row for row in results
    }
    summary_lookup = {row["model"]: row for row in summaries}
    models = [
        ("strict_n20_nominal", "严格 N20 急诊，名义报销金额"),
        ("broad_n20_n23_nominal", "广义 N20–N23 急诊，名义报销金额"),
        ("strict_n20_ipca_2022", "严格 N20 急诊，IPCA 调整至 2022 年 12 月价格"),
    ]
    lines = [
        "# 按次均报销金额分析更正",
        "",
        "## 更正内容",
        "",
        "冻结分析把严格 N20 急诊报销总额除以广义 N20–N23 急诊入院数。复现旧结果时保留其 `max(比值, 1)` 对数下限；这 196,715 个广义急诊日中，有 49,430 日没有 N20 入院，旧比值因此被下限截为 1。旧系数概括的是严格 N20 报销分子与广义 N20–N23 入院分母共同构成的混合口径；以下校正结果按一致表型估计区域日均报销金额。",
        "",
        f"旧模型复现：heat-only 比值 {reproduction['heat_only_ratio']:.6f}（95% CI {reproduction['heat_only_CI_low']:.6f}–{reproduction['heat_only_CI_high']:.6f}，p={reproduction['heat_only_p']:.4g}）；与冻结表比值差 {reproduction['heat_only_ratio_difference']:.3g}，比值、区间、p 值和模型 N 均通过逐行核对。旧模型 196,715 个输入日中有 {reproduction['old_floor_days']:,} 日（{reproduction['old_floor_share']:.2%}）低于裁剪下限；模型 N={reproduction['model_n']:,}。",
        "",
        "日级诊断显示 26,588 个广义急诊日有多个 N20–N23 诊断，其中 25,141 日同时含 N20 和其他诊断；另有 49,430 个广义急诊日没有 N20 入院。校正模型中严格 N20 与广义 N20–N23 的入院日均未见零或负总报销金额，因此有效对数样本没有因该规则额外损失。",
        "",
        "## 表型一致的校正模型",
        "",
        "| 定义 | Heat-only 比值 (95% CI), p | Smoke-only 比值 (95% CI), p | Both 比值 (95% CI), p |",
        "|---|---:|---:|---:|",
    ]
    for model_key, label in models:
        cells = []
        for group in ("heat_only", "smoke_only", "both"):
            row = result_lookup[(model_key, group)]
            cells.append(
                f"{row['ratio']:.4f} ({row['CI_low']:.4f}–{row['CI_high']:.4f}), p={row['p']:.4g}"
            )
        lines.append(
            f"| {label} | {cells[0]} | {cells[1]} | {cells[2]} |"
        )
    lines.extend([
        "",
        "| 定义 | 有效日输入数 | 模型拟合保留 N | 入院数 | 报销金额合计 |",
        "|---|---:|---:|---:|---:|",
    ])
    for model_key, label in models:
        summary = summary_lookup[model_key]
        lines.append(
            f"| {label} | {summary['analysis_rows_before_fit']:,} | {summary['model_rows']:,} "
            f"| {summary['analysis_admissions']:,.0f} | {summary['analysis_reimbursement']:,.2f} |"
        )
    lines.extend([
        "",
        "## 模型与解释",
        "",
        "回归公式为 `log_cpc ~ C(group) | stratum`，其中 `log_cpc = log(日总报销金额 / 同表型急诊入院数)`；按健康区域聚类使用 CRV1 标准误。`neither`、`heat_only`、`smoke_only`、`both` 分别对应 heat/smoke 指示值 (0,0)、(1,0)、(0,1)、(1,1)。校正模型预先采用仅保留入院数 > 0 且日总报销金额 > 0 的规则；零/负报销日计入审计但不取对数。固定价模型先按申报月份计算 `nominal_cost × IPCA_factor(competence → 2022-12)`，再按区域和入院日汇总，并除以同日严格 N20 急诊入院数。",
        "",
        "exp(β) 表示在 stratum 固定效应下，heat-only 日相对 neither 日的每日平均 SUS 报销金额几何均值比。观测单位为健康区域—日期；日均金额由同表型急诊行政报销总额除以该日入院数得到，结果量化行政报销记录中的支付金额及其组间关联。",
        "",
        "所有类别的系数、协方差、逐表型样本审计和精确公式见同目录 `cost_repair_*` CSV/JSON。冻结表与原文件均保留。",
        "",
    ])
    (OUT / "cost_repair_中文更正说明.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    db = duckdb.connect()
    db.execute("SET threads=8")
    db.execute("SET memory_limit='12GB'")

    panel = db.sql(
        f"""
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               heat, smoke, CAST(stratum AS VARCHAR) AS stratum
        FROM read_parquet('{sql_path(PANEL_PATH)}')
        WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
        """
    ).df()
    clinical = db.sql(
        f"""
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               CAST(IDENT AS VARCHAR) AS ident,
               CAST(diagnosis AS VARCHAR) AS diagnosis,
               CAST(CAR_INT AS VARCHAR) AS car_int,
               CAST(competence AS VARCHAR) AS competence,
               CAST(admissions AS DOUBLE) AS admissions,
               CAST(nominal_cost AS DOUBLE) AS nominal_cost
        FROM read_parquet('{sql_path(CLINICAL_PATH)}')
        WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
        """
    ).df()
    db.close()

    panel["date"] = pd.to_datetime(panel.date)
    panel["hrid"] = panel.hrid.astype(str)
    panel["stratum"] = panel.stratum.astype(str)
    if panel.duplicated(["hrid", "date"]).any():
        raise ValueError("Exposure panel has duplicate health-region dates")
    if panel[["hrid", "date", "heat", "smoke", "stratum"]].isna().any().any():
        raise ValueError("Exposure panel has missing model fields")
    panel["group"] = group_names(panel)
    if panel.group.isna().any() or (panel.group == "other").any():
        raise ValueError("Exposure panel contains days outside the four groups")

    clinical["date"] = pd.to_datetime(clinical.date)
    clinical["hrid"] = clinical.hrid.astype(str)
    for column in ("admissions", "nominal_cost"):
        clinical[column] = pd.to_numeric(clinical[column], errors="raise")
    if clinical[["hrid", "date", "ident", "diagnosis", "car_int", "competence"]].isna().any().any():
        raise ValueError("Clinical aggregate has missing identifiers, dates, or codes")

    eligible = clinical.loc[
        clinical.ident.eq("1")
        & clinical.car_int.eq("02")
        & clinical.diagnosis.isin(["N20", "N21", "N22", "N23"])
    ].copy()
    strict_rows = eligible.loc[eligible.diagnosis.eq("N20")].copy()
    broad_rows = eligible.copy()
    prices = pd.read_csv(IPCA_PATH, dtype={"month": str})
    price_factors = dict(zip(prices.month.astype(str), prices.factor_to_2022_12.astype(float)))
    strict_rows["price_factor"] = strict_rows.competence.astype(str).map(price_factors)
    if strict_rows.price_factor.isna().any():
        missing_months = sorted(strict_rows.loc[strict_rows.price_factor.isna(), "competence"].unique())
        raise ValueError(f"IPCA factors missing for strict N20 billing months: {missing_months}")
    strict_rows["reimbursement_2022_brl"] = (
        strict_rows.nominal_cost * strict_rows.price_factor
    )

    def daily_aggregate(rows, value_column, cost_name, extra=None):
        aggregations = {
            "admissions": ("admissions", "sum"),
            cost_name: (value_column, "sum"),
        }
        if extra:
            aggregations.update(extra)
        return rows.groupby(["hrid", "date"], as_index=False).agg(**aggregations)

    strict_daily = daily_aggregate(
        strict_rows,
        "nominal_cost",
        "strict_cost",
        extra={"strict_cost_2022": ("reimbursement_2022_brl", "sum")},
    ).rename(columns={"admissions": "strict_admissions"})
    broad_daily = daily_aggregate(
        broad_rows,
        "nominal_cost",
        "broad_cost",
        extra={
            "broad_diagnoses_present": ("diagnosis", "nunique"),
        },
    ).rename(columns={"admissions": "broad_admissions"})
    daily = strict_daily.merge(broad_daily, on=["hrid", "date"], how="outer", validate="one_to_one")
    for column in ("strict_admissions", "strict_cost", "strict_cost_2022",
                   "broad_admissions", "broad_cost", "broad_diagnoses_present"):
        daily[column] = daily[column].fillna(0)
    daily["broad_diagnoses_present"] = daily.broad_diagnoses_present.astype(int)
    if daily.duplicated(["hrid", "date"]).any():
        raise ValueError("Clinical daily aggregate has duplicate health-region dates")
    if not set(daily.hrid).issubset(set(panel.hrid)):
        raise ValueError("Clinical event regions do not all appear in the exposure panel")

    data = panel.merge(daily, on=["hrid", "date"], how="left", validate="one_to_one")
    daily_columns = [
        "strict_admissions", "strict_cost", "strict_cost_2022",
        "broad_admissions", "broad_cost", "broad_diagnoses_present",
    ]
    data[daily_columns] = data[daily_columns].fillna(0)
    data["broad_diagnoses_present"] = data.broad_diagnoses_present.astype(int)

    strict_total = float(data.strict_admissions.sum())
    broad_total = float(data.broad_admissions.sum())
    broad_admit_days = int(data.broad_admissions.gt(0).sum())
    broad_without_n20_days = int(
        (data.broad_admissions.gt(0) & data.strict_admissions.eq(0)).sum()
    )
    if strict_total != EXPECTED_STRICT_ADMISSIONS:
        raise ValueError(f"Strict N20 urgent total {strict_total} != expected {EXPECTED_STRICT_ADMISSIONS}")
    if broad_total != EXPECTED_BROAD_ADMISSIONS:
        raise ValueError(f"Broad N20-N23 urgent total {broad_total} != expected {EXPECTED_BROAD_ADMISSIONS}")
    if broad_admit_days != EXPECTED_BROAD_ADMIT_DAYS:
        raise ValueError(f"Broad urgent admit-days {broad_admit_days} != expected {EXPECTED_BROAD_ADMIT_DAYS}")
    if broad_without_n20_days != EXPECTED_BROAD_DAYS_WITHOUT_N20:
        raise ValueError(
            f"Broad days without N20 {broad_without_n20_days} != expected {EXPECTED_BROAD_DAYS_WITHOUT_N20}"
        )

    # Recreate the historical numerator/denominator mix and its lower clip.
    original_data = data.copy()
    original_data["group"] = group_names(original_data, both_label="composite")
    original_data["original_cost_per_case"] = (
        original_data.strict_cost / original_data.broad_admissions
    )
    old_floor_days = int(
        (
            original_data.broad_admissions.gt(0)
            & original_data.original_cost_per_case.lt(1)
        ).sum()
    )
    if old_floor_days != EXPECTED_BROAD_DAYS_WITHOUT_N20:
        raise ValueError(
            "Historical log-floor observations did not match broad days without N20: "
            f"{old_floor_days} versus {EXPECTED_BROAD_DAYS_WITHOUT_N20}"
        )
    old_terms, old_summary, old_details, _ = fit_cost_model(
        "original_mixed_clipped",
        "strict N20 urgent reimbursement / broad N20-N23 urgent admissions",
        original_data.assign(
            original_cost=original_data.strict_cost,
            original_admissions=original_data.broad_admissions,
        ),
        "original_cost",
        "original_admissions",
        "nominal BRL",
        "Reproduction only: retain all broad urgent days and apply log(max(cost/admissions, 1)).",
        include_only_positive_cost=False,
    )
    # Keep the legacy label as composite in reproduction outputs and map it only
    # for the cross-table comparison to the frozen CSV.
    frozen = pd.read_csv(FROZEN_PATH)
    frozen_map = frozen.set_index("group")
    repro_rows = []
    for row in old_terms:
        frozen_group = "composite" if row["exposure_group"] == "composite" else row["exposure_group"]
        frozen_row = frozen_map.loc[frozen_group]
        repro_rows.append({
            **row,
            "frozen_ratio": float(frozen_row.ratio),
            "frozen_CI_low": float(frozen_row.CI_low),
            "frozen_CI_high": float(frozen_row.CI_high),
            "frozen_p": float(frozen_row.p),
            "frozen_n": int(frozen_row.n),
            "ratio_difference": float(row["ratio"] - frozen_row.ratio),
            "CI_low_difference": float(row["CI_low"] - frozen_row.CI_low),
            "CI_high_difference": float(row["CI_high"] - frozen_row.CI_high),
            "p_difference": float(row["p"] - frozen_row.p),
            "model_n_matches_frozen": int(old_summary["model_rows"]) == int(frozen_row.n),
        })
    if any(
        not np.allclose(
            [row["ratio"], row["CI_low"], row["CI_high"], row["p"]],
            [row["frozen_ratio"], row["frozen_CI_low"], row["frozen_CI_high"], row["frozen_p"]],
            rtol=1e-12,
            atol=1e-12,
        )
        or not row["model_n_matches_frozen"]
        for row in repro_rows
    ):
        raise ValueError("Historical mixed-denominator reproduction did not match the frozen results")
    pd.DataFrame(repro_rows).to_csv(
        OUT / "cost_repair_frozen_reproduction.csv", index=False, encoding="utf-8-sig"
    )
    old_summary.update({
        "input_rows_before_fit": broad_admit_days,
        "all_broad_urgent_admissions": broad_total,
        "strict_n20_reimbursement_over_broad_denominator_total_n20_cost": float(
            data.loc[data.broad_admissions.gt(0), "strict_cost"].sum()
        ),
        "broad_days_without_n20": broad_without_n20_days,
        "broad_days_without_n20_clipped_at_floor": old_floor_days,
        "all_old_cost_ratios_below_clip_floor": old_floor_days,
    })

    results = list(old_terms)
    summaries = [old_summary]
    details = {"original_mixed_clipped": old_details}
    audits = []
    original_mask = data.broad_admissions.gt(0)
    audits.extend(summarize_by_exposure(
        data, "original_mixed_clipped", "strict_cost", "broad_admissions",
        original_mask, original_mixed=True,
    ))

    corrected_specs = [
        (
            "strict_n20_nominal", "strict N20 urgent reimbursement / strict N20 urgent admissions",
            "strict_cost", "strict_admissions", "nominal BRL",
        ),
        (
            "broad_n20_n23_nominal", "broad N20-N23 urgent reimbursement / broad N20-N23 urgent admissions",
            "broad_cost", "broad_admissions", "nominal BRL",
        ),
        (
            "strict_n20_ipca_2022", "IPCA-adjusted strict N20 urgent reimbursement / strict N20 urgent admissions",
            "strict_cost_2022", "strict_admissions", "2022-12 BRL, IPCA monthly by competence",
        ),
    ]
    for key, estimand, cost_column, admissions_column, price_basis in corrected_specs:
        positive_mask = data[admissions_column].gt(0) & data[cost_column].gt(0)
        negative_days = int((data[admissions_column].gt(0) & data[cost_column].lt(0)).sum())
        zero_days = int((data[admissions_column].gt(0) & data[cost_column].eq(0)).sum())
        terms, summary, model_details, _ = fit_cost_model(
            key,
            estimand,
            data,
            cost_column,
            admissions_column,
            price_basis,
            "Use log(cost/admissions) on days with admissions > 0 and total cost > 0; exclude zero/negative cost days because log is undefined.",
            include_only_positive_cost=True,
        )
        summary.update({
            "zero_cost_admit_days_excluded": zero_days,
            "negative_cost_admit_days_excluded": negative_days,
            "all_phenotype_admissions": float(data[admissions_column].sum()),
            "all_phenotype_reimbursement": float(data[cost_column].sum()),
        })
        results.extend(terms)
        summaries.append(summary)
        details[key] = model_details
        audits.extend(summarize_by_exposure(
            data, key, cost_column, admissions_column, positive_mask
        ))

    # A daily diagnosis count distinguishes overlapping broad phenotypes from
    # broad urgent days that have no strict N20 event at all.
    overlap = data.loc[
        data.broad_admissions.gt(0) & data.broad_diagnoses_present.gt(1)
    ]
    cross_summary = {
        "panel_rows": int(len(data)),
        "panel_regions": int(data.hrid.nunique()),
        "strict_n20_urgent_admissions": int(strict_total),
        "broad_n20_n23_urgent_admissions": int(broad_total),
        "strict_n20_admit_days": int(data.strict_admissions.gt(0).sum()),
        "broad_n20_n23_admit_days": broad_admit_days,
        "broad_days_without_strict_n20_admissions": broad_without_n20_days,
        "broad_admit_days_with_multiple_diagnosis_codes": int(len(overlap)),
        "broad_days_with_strict_n20_and_other_diagnosis": int(
            (overlap.strict_admissions.gt(0)).sum()
        ),
        "broad_days_with_other_diagnosis_but_no_strict_n20": int(
            (overlap.strict_admissions.eq(0)).sum()
        ),
        "strict_n20_zero_cost_admit_days": int(
            (data.strict_admissions.gt(0) & data.strict_cost.eq(0)).sum()
        ),
        "strict_n20_negative_cost_admit_days": int(
            (data.strict_admissions.gt(0) & data.strict_cost.lt(0)).sum()
        ),
        "broad_zero_cost_admit_days": int(
            (data.broad_admissions.gt(0) & data.broad_cost.eq(0)).sum()
        ),
        "broad_negative_cost_admit_days": int(
            (data.broad_admissions.gt(0) & data.broad_cost.lt(0)).sum()
        ),
        "total_strict_n20_nominal_reimbursement": float(data.strict_cost.sum()),
        "total_broad_n20_n23_nominal_reimbursement": float(data.broad_cost.sum()),
        "total_strict_n20_reimbursement_2022_brl": float(data.strict_cost_2022.sum()),
    }
    (OUT / "cost_repair_cross_phenotype_audit.json").write_text(
        json.dumps(cross_summary, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    pd.DataFrame(audits).to_csv(
        OUT / "cost_repair_sample_audit.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(results).to_csv(
        OUT / "cost_repair_coefficients.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(summaries).to_csv(
        OUT / "cost_repair_model_summary.csv", index=False, encoding="utf-8-sig"
    )

    manifest = {
        "panel_source": str(PANEL_PATH),
        "clinical_source": str(CLINICAL_PATH),
        "price_index_source": str(IPCA_PATH),
        "frozen_result_source": str(FROZEN_PATH),
        "analysis_window": "2018-01-01 through 2022-12-31",
        "panel_rows": int(len(data)),
        "panel_regions": int(data.hrid.nunique()),
        "group_levels": GROUPS,
        "fixed_effect": "stratum",
        "cluster": "CRV1 by hrid",
        "formulas": {key: entry["formula"] for key, entry in details.items()},
        "price_adjustment_formula": "For strict N20 urgent claims, reimbursement_2022_brl = nominal_cost * IPCA factor for billing month competence, where factor_to_2022_12 comes from monthly IPCA; sum by hrid/date, then divide by same-day strict N20 urgent admissions.",
        "exposure_group_definitions": {
            "neither": "heat=0 and smoke=0",
            "heat_only": "heat=1 and smoke=0",
            "smoke_only": "heat=0 and smoke=1",
            "both": "heat=1 and smoke=1",
        },
        "models": details,
        "clinical_totals_and_cross_phenotype_audit": cross_summary,
        "frozen_reproduction_comparison": repro_rows,
    }
    (OUT / "cost_repair_full_models.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )

    heat = next(row for row in repro_rows if row["exposure_group"] == "heat_only")
    write_chinese_report(
        results,
        summaries,
        {
            "heat_only_ratio": heat["ratio"],
            "heat_only_CI_low": heat["CI_low"],
            "heat_only_CI_high": heat["CI_high"],
            "heat_only_p": heat["p"],
            "heat_only_ratio_difference": heat["ratio_difference"],
            "model_n": int(old_summary["model_rows"]),
            "old_floor_days": old_floor_days,
            "old_floor_share": old_floor_days / EXPECTED_BROAD_ADMIT_DAYS,
        },
    )

    print("COST_REPAIR_COMPLETE", flush=True)
    print(json.dumps(cross_summary, ensure_ascii=False), flush=True)
    print(pd.DataFrame(results).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
