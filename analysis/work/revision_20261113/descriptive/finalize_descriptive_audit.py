from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "revision_20261113" / "descriptive"
MEXICO_SEX_LABELS = {
    "1": "male", "2": "female", "3": "intersex", "8": "ignored", "9": "unspecified",
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def selected_counts(rows: list[dict[str, str]], diagnosis_codes: set[str], dimension: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        if row["diagnosis"] not in diagnosis_codes or row["dimension"] != dimension:
            continue
        category = row["category"]
        counts[category] = counts.get(category, 0) + int(float(row["admissions"]))
    return counts


def label_source_codes(counts: dict[str, int], labels: dict[str, str]) -> dict[str, int]:
    return {
        f"{code} ({labels.get(code, 'unmapped source code')})": count
        for code, count in counts.items()
    }


def publication_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    _, brazil = read_csv(OUT / "brazil_demographics_marginal.csv")
    _, mexico = read_csv(OUT / "mexico_demographics_marginal.csv")
    broad_codes = {"N20", "N21", "N22", "N23"}
    brazil_n = sum(v for k, v in selected_counts(brazil, broad_codes, "sex_code").items())
    mexico_sex_raw = selected_counts(mexico, broad_codes, "sex_code")
    mexico_n = sum(mexico_sex_raw.values())

    specifications = [
        ("Brazil", "2008–2022", brazil_n, "sex_code", "sex_code", "Raw SEXO source codes retained."),
        ("Brazil", "2008–2022", brazil_n, "age_group", "age_group", "Age groups derived from DATASUS COD_IDADE/IDADE; COD_IDADE=0 remains unknown."),
        ("Mexico", "2018–2022", mexico_n, "sex_code", "sex_code", "Official SAEH SEXO codebook label; source code retained."),
        ("Mexico", "2018–2022", mexico_n, "age_category", "age_category", "Official SAEH CVEEDAD units applied; year bands assigned only when CVEEDAD=3."),
    ]
    output: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "Brazil": {"scope": "N20–N23, 2008–2022", "denominator_admissions": brazil_n},
        "Mexico": {"scope": "valid-state N20–N23, 2018–2022", "denominator_admissions": mexico_n},
    }
    for country, period, denominator, source_dimension, output_dimension, interpretation in specifications:
        source_rows = brazil if country == "Brazil" else mexico
        counts = selected_counts(source_rows, broad_codes, source_dimension)
        if country == "Mexico" and source_dimension == "sex_code":
            counts = label_source_codes(counts, MEXICO_SEX_LABELS)
        ordered_categories = sorted(counts, key=lambda value: (int(value) if value.isdigit() else 10_000, value))
        if output_dimension == "age_group":
            order = {"0–14": 0, "15–44": 1, "45–64": 2, "65+": 3, "unknown_or_invalid": 4}
            ordered_categories = sorted(counts, key=lambda value: order.get(value, 99))
        elif output_dimension == "age_category":
            ordered_categories = sorted(counts, key=lambda value: (
                0 if value.startswith("0–14 years") else
                1 if value.startswith("15–44 years") else
                2 if value.startswith("45–64 years") else
                3 if value.startswith("65+ years") else
                4 if value.startswith("hours ") else
                5 if value.startswith("days ") else
                6 if value.startswith("months ") else
                7 if value.startswith("unspecified age") else
                8 if value.startswith("unreported/invalid EDAD") else
                9 if value.startswith("unmapped CVEEDAD=") else 10,
                value,
            ))
        for category in ordered_categories:
            count = counts[category]
            output.append({
                "country": country,
                "period": period,
                "diagnosis_scope": "N20–N23",
                "denominator_admissions": denominator,
                "dimension": output_dimension,
                "category_or_source_code": category,
                "admissions": count,
                "share_percent": round(100 * count / denominator, 3) if denominator else None,
                "interpretation": interpretation,
            })

    brazil_age = selected_counts(brazil, broad_codes, "age_group")
    brazil_sex = selected_counts(brazil, broad_codes, "sex_code")
    mexico_sex = selected_counts(mexico, broad_codes, "sex_code")
    mexico_age = selected_counts(mexico, broad_codes, "age_category")
    mexico_cveedad = selected_counts(mexico, broad_codes, "cveedad_code")
    summary["Brazil"].update({"sex_code_counts": brazil_sex, "age_group_counts": brazil_age})
    summary["Mexico"].update({
        "sex_code_counts": mexico_sex_raw,
        "sex_code_labels": MEXICO_SEX_LABELS,
        "labelled_sex_counts": label_source_codes(mexico_sex_raw, MEXICO_SEX_LABELS),
        "age_category_counts": mexico_age,
        "CVEEDAD_code_counts": mexico_cveedad,
    })
    return output, summary


def write_publication_table(rows: list[dict[str, object]]) -> None:
    fields = ["country", "period", "diagnosis_scope", "denominator_admissions", "dimension", "category_or_source_code", "admissions", "share_percent", "interpretation"]
    write_csv(OUT / "demographic_summary_table.csv", fields, rows)
    sections = []
    for country in ("Brazil", "Mexico"):
        country_rows = [row for row in rows if row["country"] == country]
        denominator = country_rows[0]["denominator_admissions"] if country_rows else 0
        period = country_rows[0]["period"] if country_rows else ""
        lines = [f"### {country}: N20–N23, {period} (N={int(denominator):,})", "", "| Dimension | Category / source code | Admissions | % |", "|---|---|---:|---:|"]
        for row in country_rows:
            lines.append(f"| {row['dimension']} | {row['category_or_source_code']} | {int(row['admissions']):,} | {float(row['share_percent']):.3f} |")
        sections.append("\n".join(lines))
    notes = (
        "## Notes\n\n"
        "Brazil percentages use the 2008–2022 broad N20–N23 total. Brazilian SEXO values are shown as source codes. "
        "Mexico percentages use valid-state N20–N23 admissions from 2018–2022. Mexican SEXO labels follow the official SAEH "
        "codebook while source codes remain in the detailed files. Year-based age groups use EDAD only when CVEEDAD=3; "
        "CVEEDAD=4 and 5 are listed as unmapped because the consulted codebooks do not define them. EDAD=999 and unusual "
        "values are not assigned to an age band. The detailed joint and marginal files retain year and diagnosis strata.\n"
    )
    text = "# Aggregate age and sex summary\n\n" + "\n\n".join(sections) + "\n\n" + notes
    (OUT / "demographic_summary_table.md").write_text(text, encoding="utf-8")


def update_flow() -> None:
    path = OUT / "sample_flow_summary.csv"
    fields, rows = read_csv(path)
    rows = [row for row in rows if row.get("stage") != "SP_any_removal_full_window"]
    row = {
        "country": "Brazil", "stage": "SP_any_removal_full_window", "period": "2018–2022",
        "unit": "procedure-linked admission", "source_or_input_n": "821700", "valid_geography_n": "450",
        "exposure_complete_n": "821700", "lag0_7_complete_n": "818550", "events_broad_or_named": "71833",
        "events_urgent": "", "events_elective": "", "events_after_fit": "71684", "rows_after_fit": "158873",
        "note": "71,833 full-window events; continuous lag-0–7 model has 71,684 events and 158,873 fitted rows; categorical full-window model has 71,833 events and 159,591 fitted rows.",
    }
    insert_at = next((i + 1 for i, existing in enumerate(rows) if existing.get("stage") == "SP_any_drainage_full_window"), len(rows))
    rows.insert(insert_at, row)
    write_csv(path, fields, rows)


def update_methods() -> None:
    path = OUT / "methods_definitions.csv"
    fields, rows = read_csv(path)
    for row in rows:
        if row.get("domain") == "Brazil age cache":
            row["definition"] = "The legacy K80 age cache uses 0–14, 15–44, 45–64 and 65+ groups. A separate targeted source audit aggregates DATASUS COD_IDADE/IDADE and SEXO by year and diagnosis."
            row["limits"] = "The legacy age cache is K80-specific and has no sex field; source-level age/sex counts are in the joint/marginal tables and demographic_summary_table.csv. SEXO values remain source codes."
        elif row.get("domain") == "Brazil ICU/death/bed-days":
            row["period"] = "Broad severity aggregates: source 2008–2023, reported 2008–2022; strict urgent N20 model: 2018–2022"
            row["definition"] = "The broad daily builder counts death when MORTE=1 and ICU when UTI_MES_TO>0. The separate strict urgent N20 model requires DIAG_PRINC prefix N20 and CAR_INT=02; n20_severe is ICU (UTI_MES_TO>0 or UTI_INT_TO>0) or death (MORTE=1). It compares heat-only, smoke-only and both days with neither in conditional Poisson models with health-region×year-month×weekday strata and CRV1 health-region clustering."
            row["source_or_field"] = "build_daily_severity.py; build_beddays.py; build_n20_severe.py; fit_n20_severe.py; MORTE; UTI_MES_TO; UTI_INT_TO; CAR_INT"
            row["limits"] = "The strict urgent N20 severe result has 4,967 outcome events and 20,463 retained model area-day rows; model n is not a case count. Broad aggregate and strict-N20 ICU definitions are different and remain separately labelled."
        elif row.get("domain") == "Negative controls: Brazil":
            row["period"] = "Source extraction 2008–2023; fitted exposure-complete panel 2008–2022"
            row["definition"] = "Principal diagnosis prefixes H25/H26 for cataract and S02, S12, S22, S32, S42, S52, S62, S72, S82 and S92 for fracture. Outcome-specific conditional Poisson models use same-day tmean per 1°C and fire PM2.5 per 10 μg/m³ jointly, with health-region×year-month×weekday strata and CRV1 health-region clustering; no lag terms."
            row["source_or_field"] = "build_negative_controls.py; fit_negative_controls.py; DIAG_PRINC; tmean; fire_pm25"
            row["limits"] = "The 335,928 cataract and 2,356,578 fracture values are retained model area-day rows, not admission counts. Fracture may be heat-sensitive and is not a clean negative control; a null cataract estimate does not prove absence of confounding."
    present = {row.get("domain") for row in rows}
    extra = [
        {
            "domain": "Brazil source age and sex audit", "period": "2008–2022",
            "definition": "Filtered RD aggregates use DIAG_PRINC, DT_INTER, MUNIC_RES, SEXO, COD_IDADE and IDADE. COD_IDADE=2/3 sub-year units are grouped under 0–14; code 4 uses IDADE in years; code 5 is grouped under 65+; code 0 is retained as unknown. Age and sex are grouped by year and diagnosis.",
            "source_or_field": "aggregate_case_demographics.py; DATASUS RD fields SEXO, COD_IDADE, IDADE",
            "limits": "Source counts reconcile exactly to the cached annual diagnosis panel. SEXO values remain source codes; no person or admission IDs were exported.",
        },
        {
            "domain": "Mexico source age and sex audit", "period": "2018–2022",
            "definition": "Valid-state N20–N23 Egresos admissions are aggregated by SEXO, CVEEDAD and age category. Official SAEH mappings are CVEEDAD 0=hours, 1=days, 2=months, 3=years and 9=unspecified; age bands are assigned from EDAD only when CVEEDAD=3. Official SEXO labels are 1=male, 2=female, 3=intersex, 8=ignored and 9=unspecified.",
            "source_or_field": "aggregate_case_demographics.py; Egresos fields ENTIDAD, AFECPRIN, SEXO, EDAD, CVEEDAD; SAEH codebooks INSP/RIISP V11 and V14; 2014 historical GII-B002-SAEH-02-02 v2.2 PDF",
            "limits": "Counts reconcile exactly to the valid-state baseline cache. CVEEDAD=4/5 are present in the raw data but not defined by the consulted codebooks, so they remain unmapped; EDAD=999 and unusual values receive no year-based age group. Source codes are retained and no person or admission IDs were exported.",
        },
    ]
    for replacement in extra:
        if replacement["domain"] in present:
            for row in rows:
                if row.get("domain") == replacement["domain"]:
                    row.update(replacement)
        else:
            rows.append(replacement)
    write_csv(path, fields, rows)


def update_json(summary: dict[str, object], brazil_audit: dict[str, object], mexico_audit: dict[str, object]) -> None:
    path = OUT / "descriptive_results.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    brazil = summary["Brazil"]
    mexico = summary["Mexico"]
    data["brazil_primary"]["sex_case_distribution"] = {
        "scope": brazil["scope"], "counts_by_source_code": brazil["sex_code_counts"],
        "source_code_labels": "not assigned; raw SEXO codes retained",
    }
    data["brazil_primary"]["age_group_source_distribution"] = {
        "scope": brazil["scope"], "counts": brazil["age_group_counts"],
        "source_mapping": "DATASUS COD_IDADE/IDADE; unknown/invalid codes retained when present",
    }
    data["mexico_baseline_2018_2022"]["age_sex"] = {
        "scope": mexico["scope"], "denominator_admissions": mexico["denominator_admissions"],
        "sex_code_counts": mexico["sex_code_counts"],
        "sex_code_labels": mexico["sex_code_labels"],
        "age_category_counts": mexico["age_category_counts"],
        "CVEEDAD_code_counts": mexico["CVEEDAD_code_counts"],
        "age_note": "Official SAEH codebooks map CVEEDAD 0=hours, 1=days, 2=months, 3=years, and 9=unspecified. Year-based age bands use EDAD only for CVEEDAD=3; codes 4/5 remain unmapped.",
        "sex_note": "Official SAEH codebook labels SEXO 1=male, 2=female, 3=intersex, 8=ignored, and 9=unspecified; source codes remain in the detailed files.",
        "age_codebook_url": "https://riisp.insp.mx/nada/index.php/catalog/18/variable/V11",
        "historical_age_codebook_url": "https://riisp.insp.mx/nada/files_nesstar/SAEH_SSA/Referencias/NOTAS_TECNICAS/GII-B002-SAEH-02-02.pdf",
        "sex_codebook_url": "https://riisp.insp.mx/nada/index.php/catalog/18/variable/V14",
    }
    data["demographic_source_audit"] = {
        "status": "complete", "aggregate_only": True, "person_ids_read_or_written": False,
        "Brazil": brazil_audit, "Mexico": mexico_audit,
        "summary_table_csv": "demographic_summary_table.csv", "summary_table_markdown": "demographic_summary_table.md",
    }
    data["sample_flow"]["brazil_sp_removal_full_window_events"] = 71833
    data["sample_flow"]["brazil_sp_removal_lag0_7_events"] = 71684
    data["sample_flow"]["brazil_sp_removal_lag0_7_fitted_rows"] = 158873
    data["sample_flow"]["brazil_sp_removal_categorical_fitted_rows"] = 159591
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def update_summary(summary: dict[str, object]) -> None:
    path = OUT / "descriptive_summary_zh.md"
    text = path.read_text(encoding="utf-8")
    brazil = summary["Brazil"]
    mexico = summary["Mexico"]
    brazil_sex = "、".join(f"代码{k} {v:,}" for k, v in brazil["sex_code_counts"].items())
    brazil_age = "、".join(f"{k} {v:,}" for k, v in brazil["age_group_counts"].items())
    mexico_sex = "、".join(f"代码{k}（{mexico['sex_code_labels'].get(k, '未映射代码')}）{v:,}" for k, v in mexico["sex_code_counts"].items())
    mexico_age = "、".join(f"{k} {v:,}" for k, v in mexico["age_category_counts"].items())
    age_line = (
        f"定向源核查显示：巴西2008–2022宽N20–N23中性别源代码计数为{brazil_sex}；年龄组为{brazil_age}。"
        f"墨西哥2018–2022有效州宽N20–N23中SEXO计数为{mexico_sex}；按官方CVEEDAD单位分类后的年龄类别计数为{mexico_age}。"
        "墨西哥仅在CVEEDAD=3时将EDAD按年分组；代码4/5未由所查代码本定义，保持未映射，EDAD=999及异常值不归入年龄组。"
    )
    old_age_line = "巴西当前年龄汇总缓存有0–14、15–44、45–64、65+年龄组；原始病例性别字段仍在定向源文件核查中。墨西哥病例年龄/性别也待源文件字段核查。"
    text = text.replace(old_age_line, age_line)
    old_flow = "2018–2022巴西面板全窗为821,700个区日；构造完整0–7日滞后后为818,550个区日。SP引流全窗计数84,764，滞后完整计数84,572。"
    new_flow = "2018–2022巴西面板全窗为821,700个区日；构造完整0–7日滞后后为818,550个区日。SP引流全窗计数84,764、滞后完整计数84,572；SP取石/碎石相关记录全窗71,833、滞后完整模型71,684（15.9万拟合区日行）。"
    text = text.replace(old_flow, new_flow)
    for stale_marker in (
        "\n\n原始过滤后病例源已完成定向核查：",
        "\n\n年龄/性别源病例字段核查仍在进行",
        "\n\n## 年龄与性别源字段核查",
    ):
        if stale_marker in text:
            text = text.split(stale_marker, 1)[0]
    audit_note = (
        "\n\n## 年龄与性别源字段核查\n\n" + age_line + "\n\n"
        "墨西哥SEXO标签依据官方SAEH变量记录；CVEEDAD映射同时核对了2014年历史代码本验证规则表（v2.2）。"
        "原始代码4和5分别保留为未映射类别，不根据数值EDAD猜测其单位。\n\n"
        "可直接用于补充材料的汇总表见`demographic_summary_table.md`和`demographic_summary_table.csv`；"
        "按年度、诊断和联合分层的聚合计数见`*_demographics_joint.csv`与`*_demographics_marginal.csv`。\n"
    )
    text = text.rstrip() + audit_note
    path.write_text(text, encoding="utf-8")


def finalize() -> dict[str, object]:
    rows, summary = publication_rows()
    write_publication_table(rows)
    update_flow()
    update_methods()
    audit = json.loads((OUT / "demographic_source_audit.json").read_text(encoding="utf-8"))
    audit["N20_N23_summary"] = summary
    (OUT / "demographic_source_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    update_json(summary, audit["Brazil"], audit["Mexico"])
    update_summary(summary)
    return summary


if __name__ == "__main__":
    print(json.dumps(finalize(), ensure_ascii=False, indent=2))
