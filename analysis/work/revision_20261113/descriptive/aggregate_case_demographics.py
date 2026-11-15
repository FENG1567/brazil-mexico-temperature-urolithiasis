from __future__ import annotations

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
BASE = Path(r"F:\analysis\climate_biliary_20261102\work")
RD = Path(r"D:\Brazil_Acute_Inflammation_Climate_2003_2023\hybrid_work_v1\private\health_partitions\rd")
XW = Path(r"D:\Brazil_Acute_Inflammation_Climate_2003_2023\hybrid_work_v1\reference_step03\municipality_sih6_ibge7_health_region_2023_v2023_1.csv")
MX_FILES = [
    Path(r"D:\数据库集合\_transfer\mx_stage\mx2018_egresos.csv"),
    Path(r"D:\数据库集合\_transfer\mx_stage\mx2019_egresos.txt"),
    Path(r"D:\数据库集合\_transfer\mx_stage\mx2020_egresos.txt"),
    Path(r"D:\数据库集合\_transfer\mx_stage\mx2021_egresos.txt"),
    Path(r"D:\数据库集合\_transfer\mx_stage\mx2022_egresos.txt"),
]
OUT = ROOT / "outputs" / "revision_20261113" / "descriptive"
WORK = ROOT / "work" / "revision_20261113" / "descriptive"
STATUS = WORK / "demographic_source_status.json"

OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def log(message: str) -> None:
    print(f"{stamp()} {message}", flush=True)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def write_csv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")


def safe_str(series: pd.Series, missing: str) -> pd.Series:
    s = series.astype("string").str.strip()
    return s.mask(s.isna() | s.eq(""), missing).astype(str)


def marginal_rows(frame: pd.DataFrame, country: str, category_cols: list[str]) -> pd.DataFrame:
    output: list[pd.DataFrame] = []
    denominator = frame.groupby(["year", "diagnosis"], as_index=False).admissions.sum().rename(columns={"admissions": "diagnosis_year_total"})
    for dimension in category_cols:
        g = frame.groupby(["year", "diagnosis", dimension], as_index=False, dropna=False).admissions.sum()
        g["country"] = country
        g["dimension"] = dimension
        g = g.rename(columns={dimension: "category"})
        g = g.merge(denominator, on=["year", "diagnosis"], how="left", validate="many_to_one")
        g["share_within_diagnosis_year_pct"] = 100 * g.admissions / g.diagnosis_year_total
        output.append(g[["country", "year", "diagnosis", "dimension", "category", "admissions", "diagnosis_year_total", "share_within_diagnosis_year_pct"]])
    return pd.concat(output, ignore_index=True)


def brazil_age_group(cod_idade: pd.Series, idade: pd.Series) -> pd.Series:
    code = pd.to_numeric(cod_idade, errors="coerce")
    age = pd.to_numeric(idade, errors="coerce")
    out = pd.Series("unknown_or_invalid", index=cod_idade.index, dtype="object")
    year_age = code.eq(4) & age.ge(0)
    out.loc[year_age & age.lt(15)] = "0–14"
    out.loc[year_age & age.ge(15) & age.lt(45)] = "15–44"
    out.loc[year_age & age.ge(45) & age.lt(65)] = "45–64"
    out.loc[year_age & age.ge(65)] = "65+"
    # Observed SIH codes 2 and 3 are sub-year units; code 0 is ignored.
    out.loc[code.isin([2, 3])] = "0–14"
    out.loc[code.eq(5)] = "65+"
    return out


def mexico_age_category(age_series: pd.Series, cveedad_series: pd.Series, year_series: pd.Series) -> pd.Series:
    code = safe_str(cveedad_series, "missing_or_blank")
    age = pd.to_numeric(age_series, errors="coerce")
    year = pd.to_numeric(year_series, errors="coerce")
    out = pd.Series("unknown or inconsistent age code/value", index=age_series.index, dtype="object")
    old_vintage = year.le(2019)
    new_vintage = year.ge(2020)
    revised_code_values = year.ge(2019)

    # Legacy SAEH v2.2: 0=hours, 1=days, 2=months, 3=years, 9=unknown.
    # The 2019 source has both vocabularies. Code 3 means legacy years or
    # modern days. Assign a band only when both meanings imply the same band,
    # or EDAD exceeds the documented maximum of 29 days.
    unit_masks = {
        "hours": old_vintage & code.eq("0"),
        "modern_hours": new_vintage & code.eq("2"),
        "days": (old_vintage & code.eq("1")) | (new_vintage & code.eq("3")),
        "months": (old_vintage & code.eq("2")) | (revised_code_values & code.eq("4")),
        "years": (old_vintage & code.eq("3")) | (revised_code_values & code.eq("5")),
    }
    limits = {"hours": (1, 23), "modern_hours": (1, 23), "days": (1, 29), "months": (1, 11), "years": (1, 120)}
    for unit, mask in unit_masks.items():
        lower, upper = limits[unit]
        valid = mask & age.ge(lower) & age.le(upper)
        if unit != "years":
            out.loc[valid] = "0–14 years"
            continue
        out.loc[valid & age.lt(15)] = "0–14 years"
        out.loc[valid & age.ge(15) & age.lt(45)] = "15–44 years"
        out.loc[valid & age.ge(45) & age.lt(65)] = "45–64 years"
        out.loc[valid & age.ge(65)] = "65+ years"

    # Code 9 is the documented unknown value in both consulted vintages.
    out.loc[code.eq("9")] = "unspecified age (CVEEDAD=9)"
    # Code 0 changed meaning: hours in v2.2, not specified in v3.7.
    out.loc[new_vintage & code.eq("0") & age.eq(999)] = "unreported age (CVEEDAD=0)"
    out.loc[old_vintage & code.eq("0") & age.eq(999)] = "unreported age (CVEEDAD=0)"
    out.loc[year.eq(2019) & code.eq('3') & age.between(15,29)] = 'ambiguous 2019 age unit (years or days; EDAD 15–29)'
    out.loc[code.eq("missing_or_blank")] = "missing CVEEDAD code"
    return out


def age_value_band(series: pd.Series) -> pd.Series:
    age = pd.to_numeric(series, errors="coerce")
    out = pd.Series("missing_or_non_numeric", index=series.index, dtype="object")
    out.loc[age.ge(0) & age.lt(15)] = "0–14 numeric EDAD"
    out.loc[age.ge(15) & age.lt(45)] = "15–44 numeric EDAD"
    out.loc[age.ge(45) & age.lt(65)] = "45–64 numeric EDAD"
    out.loc[age.ge(65)] = "65+ numeric EDAD"
    return out


def run_brazil(start: float) -> tuple[pd.DataFrame, dict[str, object]]:
    import duckdb

    files = list(RD.glob("*/*/*.parquet"))
    if not files:
        raise FileNotFoundError("No Brazil RD Parquet partitions matched the existing builder glob")
    total_bytes = sum(p.stat().st_size for p in files)
    log(f"Brazil filtered RD scan started: files={len(files)} bytes={total_bytes}; query reads diagnosis/date/sex/age and crosswalk columns only")
    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    con.execute("PRAGMA memory_limit='16GB'")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"PRAGMA temp_directory='{(WORK / 'demographic_duckdb_tmp').as_posix()}'")
    rd_glob = (RD / "*" / "*" / "*.parquet").as_posix()
    xw_path = XW.as_posix()
    query = f"""
    WITH selected AS (
      SELECT
        CAST(substr(replace(CAST(r.DT_INTER AS VARCHAR), '-', ''), 1, 4) AS INTEGER) AS year,
        substr(upper(trim(CAST(r.DIAG_PRINC AS VARCHAR))), 1, 3) AS diagnosis,
        COALESCE(NULLIF(trim(CAST(r.SEXO AS VARCHAR)), ''), 'missing_or_blank') AS sex_code,
        COALESCE(NULLIF(trim(CAST(r.COD_IDADE AS VARCHAR)), ''), 'missing_or_blank') AS cod_idade,
        TRY_CAST(r.IDADE AS INTEGER) AS idade
      FROM read_parquet('{rd_glob}', union_by_name=true) r
      JOIN read_csv_auto('{xw_path}') x
        ON CAST(r.MUNIC_RES AS VARCHAR)=CAST(x.municipality_sih6 AS VARCHAR)
      WHERE substr(replace(CAST(r.DT_INTER AS VARCHAR), '-', ''), 1, 4) BETWEEN '2008' AND '2022'
        AND substr(upper(trim(CAST(r.DIAG_PRINC AS VARCHAR))), 1, 3)
            IN ('K80','K81','K85','N20','N21','N22','N23')
    ), binned AS (
      SELECT year, diagnosis, sex_code, cod_idade,
        CASE
          WHEN cod_idade='4' AND idade BETWEEN 0 AND 14 THEN '0–14'
          WHEN cod_idade='4' AND idade BETWEEN 15 AND 44 THEN '15–44'
          WHEN cod_idade='4' AND idade BETWEEN 45 AND 64 THEN '45–64'
          WHEN cod_idade='4' AND idade >= 65 THEN '65+'
          WHEN cod_idade IN ('2','3') THEN '0–14'
          WHEN cod_idade='5' THEN '65+'
          ELSE 'unknown_or_invalid'
        END AS age_group
      FROM selected
    )
    SELECT year, diagnosis, sex_code, cod_idade, age_group, COUNT(*)::BIGINT AS admissions
    FROM binned
    GROUP BY 1,2,3,4,5
    ORDER BY 1,2,3,4,5
    """
    t0 = time.perf_counter()
    frame = con.execute(query).df()
    scan_seconds = time.perf_counter() - t0
    con.close()
    frame["country"] = "Brazil"
    write_csv(frame, "brazil_demographics_joint.csv")
    marginal = marginal_rows(frame, "Brazil", ["sex_code", "age_group", "cod_idade"])
    write_csv(marginal, "brazil_demographics_marginal.csv")
    totals = frame.groupby(["year", "diagnosis"], as_index=False).admissions.sum()
    primary_map = {"K80": "K80_events", "K81": "K81_events", "K85": "K85_events", "N20": "N20_N23_events", "N21": "N20_N23_events", "N22": "N20_N23_events", "N23": "N20_N23_events"}
    panel = pd.read_csv(OUT / "brazil_primary_annual.csv")
    panel_rows = []
    for _, row in totals.iterrows():
        outcome = primary_map[row.diagnosis]
        target = int(panel.loc[panel.year.eq(int(row.year)), outcome].iloc[0])
        if row.diagnosis in ["N21", "N22", "N23"]:
            continue
        if row.diagnosis == "N20":
            source_group = int(totals.loc[(totals.year.eq(row.year)) & totals.diagnosis.isin(["N20", "N21", "N22", "N23"]), "admissions"].sum())
            target_value = int(panel.loc[panel.year.eq(int(row.year)), "N20_N23_events"].iloc[0])
            panel_rows.append({"country": "Brazil", "year": int(row.year), "diagnosis": "N20–N23", "source_admissions": source_group, "panel_admissions": target_value, "difference": source_group - target_value})
        else:
            panel_rows.append({"country": "Brazil", "year": int(row.year), "diagnosis": row.diagnosis, "source_admissions": int(row.admissions), "panel_admissions": target, "difference": int(row.admissions) - target})
    reconciliation = pd.DataFrame(panel_rows)
    write_csv(reconciliation, "brazil_demographics_source_reconciliation.csv")
    details = {
        "file_count": len(files), "bytes": total_bytes, "scan_seconds": round(scan_seconds, 2),
        "aggregate_joint_rows": int(len(frame)), "target_admissions": int(frame.admissions.sum()),
        "sex_codes": sorted(frame.sex_code.astype(str).unique().tolist()),
        "sex_code_labels": {"1": "male", "3": "female"},
        "sex_note": "SIH-RD source code 1=male and 3=female; codes are retained in detailed tables.",
        "age_codes": sorted(frame.cod_idade.astype(str).unique().tolist()),
        "age_code_interpretation": {"0": "ignored/unknown", "2": "sub-year unit", "3": "sub-year unit", "4": "years in IDADE", "5": "100 years or older"},
        "age_note": "COD_IDADE=4 uses IDADE in completed years; observed codes 2/3 are sub-year units and grouped under one year; code 5 is 100+; code 0 is retained as unknown.",
        "age_groups": sorted(frame.age_group.astype(str).unique().tolist()),
        "panel_reconciliation_max_abs_difference": int(reconciliation.difference.abs().max()) if len(reconciliation) else None,
        "source_fields_present": ["SEXO", "COD_IDADE", "IDADE", "DIAG_PRINC", "DT_INTER", "MUNIC_RES"],
        "ids_read_or_written": False,
    }
    log(f"Brazil aggregate complete: admissions={details['target_admissions']} elapsed_seconds={details['scan_seconds']} reconciliation_max_abs={details['panel_reconciliation_max_abs_difference']}")
    return frame, details


def run_mexico() -> tuple[pd.DataFrame, dict[str, object]]:
    required = ["INGRE", "ENTIDAD", "AFECPRIN", "SEXO", "EDAD", "CVEEDAD"]
    missing_files = [str(p.name) for p in MX_FILES if not p.exists()]
    if missing_files:
        raise FileNotFoundError(f"Missing Mexico Egresos source files: {missing_files}")
    sizes = {p.name: int(p.stat().st_size) for p in MX_FILES}
    total_bytes = sum(sizes.values())
    log(f"Mexico Egresos scan started: files={len(MX_FILES)} bytes={total_bytes}; reads only {required}; no ID/row-level output")
    parts: list[pd.DataFrame] = []
    source_rows = 0
    target_rows = 0
    chunk_size = 500_000
    for path in MX_FILES:
        year_hint = int(path.name[2:6])
        delim = "," if year_hint in (2018, 2019) else "|"
        t_file = time.perf_counter()
        part_rows = 0
        part_target = 0
        for chunk_number, chunk in enumerate(pd.read_csv(
            path, sep=delim, encoding="utf-8-sig", usecols=required,
            dtype="string", chunksize=chunk_size, on_bad_lines="error",
        ), start=1):
            source_rows += len(chunk)
            part_rows += len(chunk)
            dx = safe_str(chunk["AFECPRIN"], "missing").str.upper().str[:3]
            state = pd.to_numeric(chunk["ENTIDAD"], errors="coerce")
            date = pd.to_datetime(chunk["INGRE"], errors="coerce", format="mixed")
            years = date.dt.year
            keep = dx.isin(["N20", "N21", "N22", "N23"]) & state.between(1, 32) & years.between(2018, 2022)
            if not keep.any():
                if chunk_number % 4 == 0:
                    log(f"Mexico {year_hint}: parsed_chunks={chunk_number} source_rows={part_rows} target_rows={part_target}")
                continue
            selected = chunk.loc[keep, ["SEXO", "EDAD", "CVEEDAD"]].copy()
            selected["year"] = years.loc[keep].astype(int).to_numpy()
            selected["diagnosis"] = dx.loc[keep].to_numpy()
            selected["sex_code"] = safe_str(selected.SEXO, "missing_or_blank").to_numpy()
            selected["cveedad_code"] = safe_str(selected.CVEEDAD, "missing_or_blank").to_numpy()
            selected['source_year'] = year_hint
            selected["age_category"] = mexico_age_category(selected.EDAD, selected.CVEEDAD, selected.source_year).to_numpy()
            counts = selected.groupby(["year", "source_year", "diagnosis", "sex_code", "cveedad_code", "age_category"], as_index=False, dropna=False).size()
            counts = counts.rename(columns={"size": "admissions"})
            parts.append(counts)
            part_target += int(len(selected))
            target_rows += int(len(selected))
            if chunk_number % 4 == 0:
                log(f"Mexico {year_hint}: parsed_chunks={chunk_number} source_rows={part_rows} target_rows={part_target}")
        log(f"Mexico {year_hint} complete: source_rows={part_rows} valid_state_target_admissions={part_target} elapsed_seconds={round(time.perf_counter()-t_file,2)}")

    raw = pd.concat(parts, ignore_index=True)
    group_cols = ["year", "source_year", "diagnosis", "sex_code", "cveedad_code", "age_category"]
    frame = raw.groupby(group_cols, as_index=False, dropna=False).admissions.sum()
    frame["country"] = "Mexico"
    write_csv(frame, "mexico_demographics_joint.csv")
    code_year = frame.groupby(["year", "source_year", "cveedad_code"], as_index=False, dropna=False).admissions.sum()
    code_year["country"] = "Mexico"
    code_year["age_codebook_vintage"] = code_year.source_year.map(
        lambda y: "legacy 2014 SAEH v2.2" if y == 2018 else
        "mixed 2019 source; code3 EDAD15–29 unit ambiguous" if y == 2019 else
        "GIIS-B002-03-07 v3.7 (December 2020)"
    )
    code_year = code_year[["country", "year", "source_year", "cveedad_code", "admissions", "age_codebook_vintage"]]
    write_csv(code_year, "mexico_cveedad_code_year_counts.csv")
    marginal = marginal_rows(frame, "Mexico", ["sex_code", "cveedad_code", "age_category"])
    write_csv(marginal, "mexico_demographics_marginal.csv")
    valid_totals = frame.groupby(["year", "diagnosis"], as_index=False).admissions.sum()
    baseline_total = int(valid_totals.admissions.sum())
    cached = pd.read_csv(OUT / "mexico_annual_diagnosis.csv")
    cached_totals = cached.groupby(["year", "diagnosis"], as_index=False).admissions.sum()
    rec = valid_totals.merge(cached_totals, on=["year", "diagnosis"], how="outer", suffixes=("_raw", "_cache"), validate="one_to_one").fillna(0)
    rec["difference"] = rec.admissions_raw.astype(int) - rec.admissions_cache.astype(int)
    rec["country"] = "Mexico"
    write_csv(rec[["country", "year", "diagnosis", "admissions_raw", "admissions_cache", "difference"]], "mexico_demographics_source_reconciliation.csv")
    details = {
        "files": sizes, "file_count": len(MX_FILES), "bytes": total_bytes,
        "source_rows_scanned": source_rows, "valid_state_target_admissions": target_rows,
        "aggregate_joint_rows": int(len(frame)), "baseline_admissions": baseline_total,
        "sex_codes": sorted(frame.sex_code.astype(str).unique().tolist()),
        "cveedad_codes": sorted(frame.cveedad_code.astype(str).unique().tolist()),
        "age_categories": sorted(frame.age_category.astype(str).unique().tolist()),
        "cveedad_code_year_counts": code_year.to_dict(orient="records"),
        "cache_reconciliation_max_abs_difference": int(rec.difference.abs().max()) if len(rec) else None,
        "source_fields_present": required,
        "ids_read_or_written": False,
        "age_note": "Source-file vintage mapping: 2018 legacy units; 2019 code4/5 months/years and legacy code1/2; 2020–2022 modern units code2=hours, code3=days, code4=months, code5=years. The mixed 2019 source code3 EDAD15–29 cannot distinguish years from days and remains explicitly ambiguous. Code3 EDAD1–14 belongs in 0–14 under either unit; EDAD30–120 can only be valid years. All mappings follow source-file year, not admission year. EDAD999, inconsistent values and undocumented codes remain unknown; validation bounds hours1–23, days1–29, months1–11, years1–120.",
        "sex_note": "Official INSP/RIISP SAEH codebook maps SEXO 1=male, 2=female, 3=intersex, 8=ignored, 9=unspecified; detailed files retain source codes.",
        "age_codebook_url": "https://riisp.insp.mx/nada/index.php/catalog/18/variable/V11",
        "historical_age_codebook_url": "https://riisp.insp.mx/nada/files_nesstar/SAEH_SSA/Referencias/NOTAS_TECNICAS/GII-B002-SAEH-02-02.pdf",
        "current_age_codebook_url": "https://gobi.salud.gob.mx/gobi/guias/egresos_urgencias_lesiones/egresos/GIIS-B002-03-07_20210220.pdf?v=2022.01.13",
        "sex_codebook_url": "https://riisp.insp.mx/nada/index.php/catalog/18/variable/V14",
    }
    log(f"Mexico aggregate complete: source_rows={source_rows} admissions={baseline_total} reconciliation_max_abs={details['cache_reconciliation_max_abs_difference']}")
    return frame, details


def main() -> None:
    t0 = time.perf_counter()
    write_json(STATUS, {"status": "running", "pid": os.getpid(), "started_at_utc": stamp(), "outputs_are_aggregate_only": True})
    try:
        brazil_frame, brazil_meta = run_brazil(t0)
        mexico_frame, mexico_meta = run_mexico()
        brazil_total = int(brazil_frame.admissions.sum())
        mexico_total = int(mexico_frame.admissions.sum())
        if int(brazil_meta["panel_reconciliation_max_abs_difference"]) != 0:
            raise ValueError("Brazil demographic-source counts do not reconcile to the cached annual primary panel")
        if int(mexico_meta["cache_reconciliation_max_abs_difference"]) != 0:
            raise ValueError("Mexico Egresos target counts do not reconcile to the audited 2018–2022 cache")
        summary = {
            "status": "complete", "started_at_utc": json.loads(STATUS.read_text(encoding="utf-8"))["started_at_utc"],
            "completed_at_utc": stamp(), "pid": os.getpid(), "elapsed_seconds": round(time.perf_counter()-t0, 2),
            "aggregate_only": True, "person_ids_read_or_written": False,
            "Brazil": brazil_meta, "Mexico": mexico_meta,
            "brazil_admissions_total": brazil_total, "mexico_admissions_total": mexico_total,
        }
        write_json(OUT / "demographic_source_audit.json", summary)
        write_json(STATUS, summary)
        # Extend the existing Chinese explanation without replacing its earlier descriptive sections.
        md_path = OUT / "descriptive_summary_zh.md"
        md = md_path.read_text(encoding="utf-8")
        md = md.replace(
            "年龄/性别源病例字段核查仍在进行；如能由已有原始或过滤源完成聚合，将追加去标识汇总表并更新本说明。",
            f"""原始过滤后病例源已完成定向核查：Brazil RD共{brazil_meta['file_count']:,}个Parquet文件（{brazil_meta['bytes']:,}字节），Mexico Egresos 2018–2022共{mexico_meta['file_count']}个文件（{mexico_meta['bytes']:,}字节）。只选择诊断、入院日、地理有效标志、性别和年龄字段并直接聚合，没有读取或写出个人/住院ID。Brazil按DATASUS `COD_IDADE`/`IDADE`生成年龄组；Mexico保留`CVEEDAD`代码并按`EDAD`数值分组，不擅自解释年龄代码。Brazil聚合病例数与既有主面板逐年诊断计数最大差为0；Mexico有效州N20–N23总数与既有基线缓存最大差为0。

`brazil_demographics_joint.csv`、`brazil_demographics_marginal.csv`、`mexico_demographics_joint.csv`和`mexico_demographics_marginal.csv`提供去标识年度、诊断、性别代码与年龄组计数；`demographic_source_reconciliation.csv`和`demographic_source_audit.json`记录来源核对及读取规模。"""
        )
        md_path.write_text(md, encoding="utf-8")
        log(f"Demographic audit complete: Brazil={brazil_total} Mexico={mexico_total} total_elapsed={summary['elapsed_seconds']} seconds")
    except Exception as exc:
        write_json(STATUS, {"status": "failed", "pid": os.getpid(), "updated_at_utc": stamp(), "error": repr(exc)})
        raise


if __name__ == "__main__":
    sys.path.insert(0, str(BASE / "pylibs"))
    main()
