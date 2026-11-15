"""Count procedure-code feasibility in strict N20 Mexico admissions.

Admission identifiers are used only in memory to recover source-file year and
join annual procedure tables by ID + source-file year. Admission date remains
the exposure date; outputs contain aggregate counts and no identifiers.
"""

from __future__ import annotations

import json
import csv
import hashlib
import shutil
import sys
from pathlib import Path

sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
import pandas as pd
import audit_mexico_source_partitions as partition_audit


BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs" / "upgrade_20261112"
COHORT_PATH = Path(
    r"F:\analysis\climate_biliary_20261102\work\mexico_n20_id_proced.parquet"
)
PROCEDURE_ROOT = Path(r"D:\数据库集合\_transfer\mx_stage")
CHUNK_SIZE = 250_000

EXPECTED_COHORT_BY_YEAR = {
    2018: 7_446,
    2019: 8_552,
    2020: 3_602,
    2021: 5_367,
    2022: 7_450,
}

PROCEDURE_FILES = {
    2018: (PROCEDURE_ROOT / "mx2018_proced.csv", ","),
    2019: (PROCEDURE_ROOT / "mx2019_proced.txt", ","),
    2020: (PROCEDURE_ROOT / "mx2020_proced.txt", "|"),
    2021: (PROCEDURE_ROOT / "mx2021_proced.txt", "|"),
    2022: (PROCEDURE_ROOT / "mx2022_proced.txt", "|"),
}

SOURCE_FILE_YEAR_FILES = {
    year: (
        partition_audit.ROOT
        / (f"mx{year}_egresos." + ("csv" if year == 2018 else "txt")),
        "," if year in (2018, 2019) else "|",
    )
    for year in range(2018, 2023)
}

EXPECTED_CROSS_PARTITION_COUNTS = {
    (2019, 2018): 20,
    (2020, 2019): 20,
    (2021, 2020): 5,
    (2022, 2021): 13,
}
EXPECTED_STRICT_COHORT_TOTAL = 32_417

PRE_REPAIR_OUTPUT_NAMES = (
    "mexico_procedure_feasibility.json",
    "mexico_procedure_feasibility.md",
    "mexico_specific_daily_counts.csv",
    "mexico_specific_model_results.csv",
    "mexico_specific_full_covariance.json",
    "mexico_specific_codebook.csv",
    "mexico_specific_provenance.json",
)
PRE_REPAIR_DIR = OUT / "pre_partition_linkage_repair"

# These exact CIE-9-MC catalogue values were fixed before counting.
# 598X is the literal catalogue code/title for ureteral catheterization.
CODE_GROUPS = {
    "nephrostomy": {
        "label": "Nephrostomy",
        "codes": {"5502", "5503", "5504"},
    },
    "ureteral_catheterization": {
        "label": "Ureteral catheterization (not stent-specific)",
        "codes": {"598X"},
    },
    "urinary_stone_fragmentation_lithotripsy": {
        "label": "Urinary stone fragmentation/lithotripsy",
        "codes": {"5504", "5995", "9851"},
    },
    "urinary_drainage_device_exchange_irrigation_removal": {
        "label": "Urinary drainage-device exchange/irrigation/removal",
        "codes": {"5593", "9645", "9646", "9761", "9762", "9764"},
    },
}

REQUIRED_PROCEDURE_COLUMNS = {"ID", "NUMPROMED", "PROMED", "TIPO"}


def preserve_pre_repair_outputs() -> None:
    """Keep prior aggregate outputs with the superseded partition rule disclosed."""
    manifest_path = PRE_REPAIR_DIR / "preservation_manifest.json"
    if manifest_path.exists():
        return
    PRE_REPAIR_DIR.mkdir(parents=True, exist_ok=True)
    preserved = []
    for name in PRE_REPAIR_OUTPUT_NAMES:
        source = OUT / name
        destination = PRE_REPAIR_DIR / name
        if source.is_file() and not destination.exists():
            shutil.copy2(source, destination)
        if destination.is_file():
            digest = hashlib.sha256()
            with destination.open("rb") as archived:
                for block in iter(lambda: archived.read(1024 * 1024), b""):
                    digest.update(block)
            preserved.append(
                {
                    "file": name,
                    "bytes": destination.stat().st_size,
                    "sha256": digest.hexdigest(),
                }
            )
    manifest = {
        "purpose": "Preserve aggregate procedure outputs produced before source-partition repair.",
        "superseded_link_rule": (
            "Procedure source file selected by admission date year."
        ),
        "audit_finding": (
            "The strict N20 source audit found 58 of 32,417 admissions in an annual "
            "egresos partition whose file year differs from admission year: "
            "20 source-2019/admission-2018, 20 source-2020/admission-2019, "
            "5 source-2021/admission-2020, and 13 source-2022/admission-2021. "
            "Prior linked counts and models may therefore be affected."
        ),
        "preserved_outputs": preserved,
        "raw_admission_identifiers_or_records_copied": False,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def recover_source_file_years(cohort: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Bind cohort admissions to egresos partitions in memory using the audit signature."""
    required = {"ID", "INGRE", "ENTIDAD", "AFECPRIN", "year"}
    missing = required.difference(cohort.columns)
    if missing:
        raise ValueError("Cohort lacks fields needed for source-partition recovery")
    exposure_dates = pd.to_datetime(cohort["date"], errors="coerce").dt.strftime(
        "%Y-%m-%d"
    )
    source_dates = cohort["INGRE"].astype("string").str.strip().str[:10]
    if exposure_dates.isna().any() or not exposure_dates.eq(source_dates).all():
        raise ValueError("Admission exposure dates disagree with the audit admission date")

    signature_to_admission: dict[tuple[str, str, str, str], tuple[str, int]] = {}
    for row in cohort[["ID", "INGRE", "ENTIDAD", "AFECPRIN", "year"]].itertuples(
        index=False, name=None
    ):
        identifier, raw_date, state, diagnosis, admission_year = row
        signature = partition_audit.signature(identifier, raw_date, state, diagnosis)
        key = (str(identifier).strip(), int(admission_year))
        if signature in signature_to_admission:
            raise ValueError("Cohort source signatures are not unique after deduplication")
        signature_to_admission[signature] = key

    source_year_by_admission: dict[tuple[str, int], int] = {}
    source_rows_by_year = {}
    for source_year, (path, delimiter) in SOURCE_FILE_YEAR_FILES.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing annual egresos source for {source_year}")
        records = blank_lines = footer_count = 0
        footer_total = None
        with path.open("r", encoding="utf-8-sig", newline="") as source:
            reader = csv.reader(source, delimiter=delimiter, strict=True)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError(f"Empty annual egresos source for {source_year}") from exc
            required_source_columns = {"ID", "INGRE", "ENTIDAD", "AFECPRIN"}
            if not required_source_columns.issubset(header):
                raise ValueError(f"Annual egresos source {path.name} lacks required columns")
            indexes = {name: header.index(name) for name in required_source_columns}
            for row in reader:
                if not row or (len(row) == 1 and not row[0].strip()):
                    blank_lines += 1
                    continue
                trailer = (
                    partition_audit.FOOTER.fullmatch(row[0].strip())
                    if len(row) == 1
                    else None
                )
                if trailer:
                    if footer_count:
                        raise ValueError(f"Duplicate SQL-export trailer in {source_year}")
                    footer_total = int(trailer.group(1))
                    footer_count += 1
                    continue
                if len(row) != len(header) or footer_count:
                    raise ValueError(
                        f"Unexpected annual egresos structure in {source_year} "
                        f"at line {reader.line_num}"
                    )
                records += 1
                signature = partition_audit.signature(
                    row[indexes["ID"]],
                    row[indexes["INGRE"]],
                    row[indexes["ENTIDAD"]],
                    row[indexes["AFECPRIN"]],
                )
                admission = signature_to_admission.get(signature)
                if admission is not None:
                    if admission in source_year_by_admission:
                        raise ValueError(
                            "A cohort admission matched more than one source partition"
                        )
                    source_year_by_admission[admission] = source_year
        if footer_total is not None and footer_total != records:
            raise ValueError(
                f"SQL-export trailer total mismatch in {source_year}: "
                f"{footer_total} != {records}"
            )
        source_rows_by_year[source_year] = {
            "file": path.name,
            "source_rows_parsed": records,
            "blank_lines_excluded": blank_lines,
            "sql_export_trailers_excluded": footer_count,
            "trailer_record_total": footer_total,
        }

    unmatched = len(signature_to_admission) - len(source_year_by_admission)
    if unmatched:
        raise ValueError(
            f"Source-partition recovery left {unmatched} strict cohort admissions unmatched"
        )
    if len(source_year_by_admission) != EXPECTED_STRICT_COHORT_TOTAL:
        raise ValueError(
            "Source-partition recovery did not reconcile all 32,417 strict N20 admissions"
        )

    cross_counts: dict[tuple[int, int], int] = {}
    for (identifier, admission_year), source_year in source_year_by_admission.items():
        if source_year != admission_year:
            pair = (source_year, admission_year)
            cross_counts[pair] = cross_counts.get(pair, 0) + 1
    if cross_counts != EXPECTED_CROSS_PARTITION_COUNTS:
        raise ValueError(
            "Recovered source/admission-year mismatches differ from the audited "
            f"partition counts: {cross_counts}"
        )

    recovered = cohort.copy()
    recovered["source_file_year"] = [
        source_year_by_admission[(str(identifier).strip(), int(admission_year))]
        for identifier, admission_year in zip(recovered["ID"], recovered["year"])
    ]
    if recovered.duplicated(["ID", "source_file_year"]).any():
        raise ValueError("ID + source_file_year is not unique for the strict cohort")

    summary = {
        "strict_admissions_reconciled": len(source_year_by_admission),
        "unmatched_admissions": 0,
        "ambiguous_admissions": 0,
        "matching_signature": (
            "ID + ISO admission date + ENTIDAD + AFECPRIN, in memory only"
        ),
        "procedure_link_key": "ID + source_file_year",
        "exposure_date_field": "date (admission date)",
        "cross_year_counts": [
            {
                "source_file_year": source_year,
                "admission_year": admission_year,
                "admissions": count,
            }
            for (source_year, admission_year), count in sorted(cross_counts.items())
        ],
        "source_rows_by_year": source_rows_by_year,
    }
    return recovered, summary


def load_cohort() -> tuple[pd.DataFrame, dict[int, set[str]], dict]:
    db = duckdb.connect()
    db.execute('SET threads=8')
    source = db.execute('SELECT * FROM read_parquet(?)', [str(COHORT_PATH)]).df()
    required = {"ID", "INGRE", "ENTIDAD", "AFECPRIN", "date"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("Cohort source is missing required columns")

    dates = pd.to_datetime(source["date"], errors="coerce")
    years = dates.dt.year
    diagnoses = source["AFECPRIN"].astype("string").str[:3]
    state_numbers = pd.to_numeric(source["ENTIDAD"], errors="coerce")
    strict_in_window = diagnoses.eq("N20") & years.between(2018, 2022)
    strict = source.loc[
        strict_in_window, ["ID", "INGRE", "ENTIDAD", "AFECPRIN", "date"]
    ].copy()
    strict["year"] = years.loc[strict_in_window].astype(int).to_numpy()
    strict["state_number"] = state_numbers.loc[strict_in_window].to_numpy()
    valid = strict.loc[strict["state_number"].between(1, 32)].copy()
    valid["ID"] = valid["ID"].astype("string").str.strip()

    missing_id = valid["ID"].isna() | valid["ID"].eq("")
    if missing_id.any():
        raise ValueError("Valid-state strict N20 cohort contains missing IDs")

    duplicate_state_conflicts = (
        valid.groupby(["ID", "year"], observed=True)["state_number"].nunique().gt(1)
    )
    if duplicate_state_conflicts.any():
        raise ValueError("Duplicate ID-year cohort rows disagree on state")

    before_dedup = len(valid)
    cohort = valid.drop_duplicates(["ID", "year"]).copy()
    cohort_counts = {
        int(year): int(count)
        for year, count in cohort.groupby("year", observed=True).size().items()
    }
    if cohort_counts != EXPECTED_COHORT_BY_YEAR:
        raise ValueError(
            "Unexpected strict N20 valid-state cohort counts; refusing to count "
            "procedure codes"
        )

    cohort, partition_audit_summary = recover_source_file_years(cohort)
    ids_by_source_year = {
        int(source_year): set(group["ID"].astype(str))
        for source_year, group in cohort.groupby("source_file_year", observed=True)
    }
    audit = {
        "strict_n20_records_in_window_before_state_filter": int(strict.shape[0]),
        "strict_n20_records_with_valid_state_before_deduplication": int(before_dedup),
        "strict_n20_unique_id_year_admissions": int(len(cohort)),
        "valid_state_duplicate_id_year_rows_removed": int(before_dedup - len(cohort)),
        "unique_admissions_by_year": cohort_counts,
        "valid_state_definition": "ENTIDAD numeric value from 1 through 32 inclusive",
        "admission_key": "ID + admission year",
        "procedure_link_key": "ID + source_file_year",
        "source_partition_reconciliation": partition_audit_summary,
    }
    return cohort, ids_by_source_year, audit


def validate_procedure_header(path: Path, delimiter: str) -> None:
    header = pd.read_csv(
        path,
        sep=delimiter,
        encoding="utf-8-sig",
        nrows=0,
        dtype="string",
        on_bad_lines="error",
        engine="c",
    )
    if not REQUIRED_PROCEDURE_COLUMNS.issubset(header.columns):
        raise ValueError(
            f"Procedure source {path.name} lacks required columns"
        )


def scan_procedure_tables(
    ids_by_source_year: dict[int, set[str]],
) -> tuple[dict, dict, dict]:
    seen_rows: set[tuple[str, int]] = set()
    seen_coded_admissions: set[tuple[str, int]] = set()
    admissions_with_missing_code_row: set[tuple[str, int]] = set()
    any_surgery: set[tuple[str, int]] = set()
    any_matches = {key: set() for key in CODE_GROUPS}
    first_ordinal_matches = {key: set() for key in CODE_GROUPS}
    per_year_source = {}

    for source_year, (path, delimiter) in PROCEDURE_FILES.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing procedure source for {source_year}")
        validate_procedure_header(path, delimiter)

        rows_parsed = 0
        linked_rows = 0
        missing_code_rows = 0
        reader = pd.read_csv(
            path,
            sep=delimiter,
            encoding="utf-8-sig",
            usecols=sorted(REQUIRED_PROCEDURE_COLUMNS),
            dtype="string",
            keep_default_na=False,
            chunksize=CHUNK_SIZE,
            on_bad_lines="error",
            engine="c",
        )
        for chunk in reader:
            rows_parsed += len(chunk)
            ids = chunk["ID"].astype("string").str.strip()
            selected_mask = ids.isin(ids_by_source_year.get(source_year, set()))
            if not selected_mask.any():
                continue

            selected = chunk.loc[
                selected_mask, ["NUMPROMED", "PROMED", "TIPO"]
            ].copy()
            selected_ids = ids.loc[selected_mask]
            keys = [
                (str(identifier), source_year)
                for identifier in selected_ids.tolist()
            ]
            linked_rows += len(selected)
            seen_rows.update(keys)

            tipos = selected["TIPO"].astype("string").str.strip().str.upper()
            for index in selected.index[tipos.eq("Q")]:
                any_surgery.add(keys[selected.index.get_loc(index)])

            codes = selected["PROMED"].astype("string").str.strip().str.upper()
            first_ordinal = pd.to_numeric(
                selected["NUMPROMED"].astype("string").str.strip(),
                errors="coerce",
            ).eq(1)
            missing_code = codes.isna() | codes.eq("")
            missing_code_rows += int(missing_code.sum())
            for index in selected.index[missing_code]:
                admissions_with_missing_code_row.add(keys[selected.index.get_loc(index)])

            nonmissing_indices = selected.index[~missing_code]
            for index in nonmissing_indices:
                position = selected.index.get_loc(index)
                key = keys[position]
                seen_coded_admissions.add(key)
                code = codes.loc[index]
                for group_key, group in CODE_GROUPS.items():
                    if code in group["codes"]:
                        any_matches[group_key].add(key)
                        if bool(first_ordinal.loc[index]):
                            first_ordinal_matches[group_key].add(key)

        per_year_source[source_year] = {
            "file": path.name,
            "delimiter": delimiter,
            "cohort_admissions_assigned_to_source_file": int(
                len(ids_by_source_year.get(source_year, set()))
            ),
            "procedure_rows_parsed": int(rows_parsed),
            "procedure_rows_linked_to_cohort": int(linked_rows),
            "linked_admissions_with_any_procedure_row": int(
                sum(1 for key in seen_rows if key[1] == source_year)
            ),
            "linked_admissions_with_any_nonblank_PROMED": int(
                sum(1 for key in seen_coded_admissions if key[1] == source_year)
            ),
            "linked_admissions_with_at_least_one_blank_PROMED_row": int(
                sum(
                    1
                    for key in admissions_with_missing_code_row
                    if key[1] == source_year
                )
            ),
            "cohort_admissions_with_no_nonblank_PROMED": int(
                len(ids_by_source_year.get(source_year, set()))
                - sum(
                    1
                    for key in seen_coded_admissions
                    if key[1] == source_year
                )
            ),
            "linked_procedure_rows_missing_PROMED": int(missing_code_rows),
            "admissions_with_any_TIPO_Q_row": int(
                sum(1 for key in any_surgery if key[1] == source_year)
            ),
        }

    source_totals = {
        "procedure_rows_parsed": int(
            sum(item["procedure_rows_parsed"] for item in per_year_source.values())
        ),
        "procedure_rows_linked_to_cohort": int(
            sum(item["procedure_rows_linked_to_cohort"] for item in per_year_source.values())
        ),
        "linked_procedure_rows_missing_PROMED": int(
            sum(item["linked_procedure_rows_missing_PROMED"] for item in per_year_source.values())
        ),
    }
    return (
        {
            "linked_admissions_with_any_procedure_row": int(len(seen_rows)),
            "linked_admissions_with_any_nonblank_PROMED": int(
                len(seen_coded_admissions)
            ),
            "cohort_admissions_with_no_nonblank_PROMED": int(
                sum(len(ids) for ids in ids_by_source_year.values())
                - len(seen_coded_admissions)
            ),
            "linked_admissions_with_at_least_one_blank_PROMED_row": int(
                len(admissions_with_missing_code_row)
            ),
            "admissions_with_any_TIPO_Q_row": int(len(any_surgery)),
            **source_totals,
        },
        per_year_source,
        {
            "any_recorded_matching_procedure": any_matches,
            "first_ordinal_NUMPROMED_1": first_ordinal_matches,
            "any_surgery": any_surgery,
        },
    )


def make_group_counts(
    cohort: pd.DataFrame,
    cohort_counts: dict[int, int],
    match_sets: dict,
) -> list[dict]:
    rows = []
    years = sorted(cohort_counts)
    admission_year_by_key = {
        (str(identifier).strip(), int(source_year)): int(admission_year)
        for identifier, source_year, admission_year in zip(
            cohort["ID"], cohort["source_file_year"], cohort["year"]
        )
    }
    for group_key in CODE_GROUPS:
        if not match_sets["any_recorded_matching_procedure"][group_key].issubset(
            admission_year_by_key
        ):
            raise ValueError("Procedure match set contains an admission outside the cohort")
        if not match_sets["first_ordinal_NUMPROMED_1"][group_key].issubset(
            match_sets["any_recorded_matching_procedure"][group_key]
        ):
            raise ValueError("First-ordinal matches exceed any-code matches")

    def count_in_admission_year(keys: set[tuple[str, int]], year: int) -> int:
        return sum(1 for key in keys if admission_year_by_key[key] == year)

    for group_key, group in CODE_GROUPS.items():
        for year in years:
            denominator = cohort_counts[year]
            any_count = count_in_admission_year(
                match_sets["any_recorded_matching_procedure"][group_key], year
            )
            first_count = count_in_admission_year(
                match_sets["first_ordinal_NUMPROMED_1"][group_key], year
            )
            if first_count > any_count:
                raise ValueError("First-ordinal matches exceed any-code matches")
            rows.append(
                {
                    "group_key": group_key,
                    "group_label": group["label"],
                    "year": year,
                    "strict_n20_admissions": denominator,
                    "admissions_any_matching_code": int(any_count),
                    "admissions_first_ordinal_num_promed_1": int(first_count),
                }
            )
        denominator = sum(cohort_counts.values())
        total_any = len(match_sets["any_recorded_matching_procedure"][group_key])
        total_first = len(match_sets["first_ordinal_NUMPROMED_1"][group_key])
        rows.append(
            {
                "group_key": group_key,
                "group_label": group["label"],
                "year": "total",
                "strict_n20_admissions": denominator,
                "admissions_any_matching_code": int(total_any),
                "admissions_first_ordinal_num_promed_1": int(total_first),
            }
        )
    return rows


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def build_report(result: dict) -> str:
    coverage_rows = []
    for source_year in sorted(
        int(key) for key in result["procedure_sources_by_file_year"]
    ):
        source = result["procedure_sources_by_file_year"][str(source_year)]
        coverage_rows.append(
            [
                str(source_year),
                f"{source['cohort_admissions_assigned_to_source_file']:,}",
                f"{source['procedure_rows_parsed']:,}",
                f"{source['procedure_rows_linked_to_cohort']:,}",
                f"{source['linked_admissions_with_any_procedure_row']:,}",
                f"{source['linked_admissions_with_any_nonblank_PROMED']:,}",
                f"{source['cohort_admissions_with_no_nonblank_PROMED']:,}",
                f"{source['linked_procedure_rows_missing_PROMED']:,}",
                f"{source['admissions_with_any_TIPO_Q_row']:,}",
            ]
        )
    coverage_rows.append(
        [
            "Total",
            f"{result['cohort']['strict_n20_unique_id_year_admissions']:,}",
            f"{result['coverage']['procedure_rows_parsed']:,}",
            f"{result['coverage']['procedure_rows_linked_to_cohort']:,}",
            f"{result['coverage']['linked_admissions_with_any_procedure_row']:,}",
            f"{result['coverage']['linked_admissions_with_any_nonblank_PROMED']:,}",
            f"{result['coverage']['cohort_admissions_with_no_nonblank_PROMED']:,}",
            f"{result['coverage']['linked_procedure_rows_missing_PROMED']:,}",
            f"{result['coverage']['admissions_with_any_TIPO_Q_row']:,}",
        ]
    )

    group_rows = []
    for row in result["procedure_group_counts_by_year"]:
        group_rows.append(
            [
                row["group_label"],
                str(row["year"]),
                f"{row['admissions_any_matching_code']:,}",
                f"{row['admissions_first_ordinal_num_promed_1']:,}",
            ]
        )

    group_definitions = [
        f"- {group['label']}: {', '.join(sorted(group['codes']))}."
        for group in CODE_GROUPS.values()
    ]
    group_table = markdown_table(
        [
            "Code group",
            "Admission year",
            "Admissions with any matching code",
            "Admissions with matching code at NUMPROMED=1",
        ],
        group_rows,
    )
    coverage_table = markdown_table(
        [
            "Procedure source file year",
            "Strict N20 admissions assigned to source file",
            "Procedure-table rows parsed",
            "Linked procedure rows",
            "Admissions with a linked procedure row",
            "Admissions with any nonblank PROMED",
            "Admissions with no nonblank PROMED",
            "Linked rows with blank PROMED",
            "Admissions with any TIPO=Q row",
        ],
        coverage_rows,
    )
    total = result["cohort"]["strict_n20_unique_id_year_admissions"]
    linked = result["coverage"]["linked_admissions_with_any_procedure_row"]
    coded = result["coverage"]["linked_admissions_with_any_nonblank_PROMED"]
    no_code = result["coverage"]["cohort_admissions_with_no_nonblank_PROMED"]
    missing_admissions = result["coverage"][
        "linked_admissions_with_at_least_one_blank_PROMED_row"
    ]
    missing_code_rows = result["coverage"]["linked_procedure_rows_missing_PROMED"]
    partition = result["cohort"]["source_partition_reconciliation"]
    cross_rows = [
        [
            str(row["source_file_year"]),
            str(row["admission_year"]),
            f"{row['admissions']:,}",
        ]
        for row in partition["cross_year_counts"]
    ]
    cross_table = markdown_table(
        ["Source file year", "Admission year", "Admissions"], cross_rows
    )

    return f"""# Mexico procedure-code feasibility

## Cohort and linkage

The cohort contains {total:,} unique 2018–2022 admissions with a principal diagnosis beginning N20 and a valid state code (ENTIDAD 1–32), deduplicated by ID and admission year. Admission identifiers were used only in memory. The annual egresos source file was recovered with the audit signature ID + ISO admission date + ENTIDAD + AFECPRIN; all {partition['strict_admissions_reconciled']:,} strict admissions matched exactly one source partition, with no unmatched or ambiguous admission. Procedure records are linked by ID + source file year. Admission date remains the exposure date.

The previous source-selection rule chose the procedure file by admission year. That rule differs from the audited source file year for 58 admissions, so pre-repair linked counts and models may be affected; those aggregate outputs are preserved in `pre_partition_linkage_repair` with the old rule documented. Annual procedure tables were parsed with strict row handling. The discharge-table field PROCED was not treated as a procedure code; the procedure code field used here is PROMED. Candidate code matches include the listed PROMED values regardless of TIPO because the reviewed catalogue classifies some of these procedure codes as therapeutic (T). Counts are unique admissions, not procedure rows.

The audited cross-partition admissions were:

{cross_table}

{coverage_table}

Across the full cohort, {linked:,} admissions ({100 * linked / total:.1f}%) had at least one linked procedure row; {coded:,} admissions had at least one nonblank PROMED code. {no_code:,} cohort admissions had no linked nonblank PROMED. There were {missing_code_rows:,} linked procedure rows with blank PROMED, affecting {missing_admissions:,} unique admissions. These metrics describe missingness in the linked annual procedure tables.

## Candidate code counts

{group_table}

Candidate groups were fixed from the reviewed government CIE-9-MC catalogue:

{chr(10).join(group_definitions)}

The two count columns are separate definitions: any linked procedure row bearing a group code, and a matching code whose NUMPROMED ordinal equals 1. The catalogue dictionary describes NUMPROMED as an ordinal from 1 to N; this report calls it the first ordinal and does not interpret it as a documented primary procedure. Group counts are not mutually exclusive: code 5504 occurs in both the nephrostomy and fragmentation/lithotripsy lists.

The reviewed catalogue labels 598X as ureteral catheterization and contained no specific ureteral-stent title. Catheterization is not counted as stent placement. Fragmentation/lithotripsy codes are not interpreted as stone extraction. The annual tables do not supply a verified operation date here, so the counts establish that a code was recorded during the linked admission, not when the procedure occurred.

The `TIPO=Q` measure is a separate admission flag for any row classified as surgery; it is not a code-group definition.

## Scope

This is a feasibility and count-evidence step only. No association model was fitted. Raw admission identifiers were not written to the report, JSON, or any intermediate output.
"""


def assemble_result(
    cohort: pd.DataFrame,
    cohort_audit: dict,
    coverage: dict,
    source_by_file_year: dict,
    match_sets: dict,
) -> dict:
    cohort_counts = {
        str(int(year)): int(count)
        for year, count in cohort.groupby("year", observed=True).size().items()
    }
    group_counts = make_group_counts(
        cohort,
        {int(year): count for year, count in cohort_counts.items()},
        match_sets,
    )

    admission_year_by_key = {
        (str(identifier).strip(), int(source_year)): int(admission_year)
        for identifier, source_year, admission_year in zip(
            cohort["ID"], cohort["source_file_year"], cohort["year"]
        )
    }
    if not match_sets["any_surgery"].issubset(admission_year_by_key):
        raise ValueError("TIPO=Q flags contain an admission outside the strict cohort")
    partition = cohort_audit["source_partition_reconciliation"]
    cross_partition_total = sum(
        int(row["admissions"]) for row in partition["cross_year_counts"]
    )
    if (
        partition["strict_admissions_reconciled"] != EXPECTED_STRICT_COHORT_TOTAL
        or partition["unmatched_admissions"] != 0
        or partition["ambiguous_admissions"] != 0
        or cross_partition_total != 58
    ):
        raise ValueError("Source-partition repair audit did not pass reconciliation")
    surgery_by_admission_year = {
        str(year): sum(
            admission_year_by_key[key] == year for key in match_sets["any_surgery"]
        )
        for year in sorted(int(value) for value in cohort_counts)
    }

    return {
        "title": "Mexico strict N20 procedure-code feasibility",
        "cohort": {
            **cohort_audit,
            "unique_admissions_by_year": cohort_counts,
            "unique_admissions_by_source_file_year": {
                str(year): int(count)
                for year, count in cohort.groupby(
                    "source_file_year", observed=True
                ).size().items()
            },
            "source_file": COHORT_PATH.name,
        },
        "procedure_sources_by_file_year": {
            str(year): source_by_file_year[year]
            for year in sorted(source_by_file_year)
        },
        "coverage": coverage,
        "surgery_type_q_admissions_by_admission_year": surgery_by_admission_year,
        "procedure_code_groups": {
            group_key: {
                "label": group["label"],
                "codes": sorted(group["codes"]),
            }
            for group_key, group in CODE_GROUPS.items()
        },
        "procedure_group_counts_by_year": group_counts,
        "count_definitions": {
            "any_matching_code": (
                "Unique strict N20 admissions with at least one linked PROMED value "
                "equal to a code in the group, regardless of TIPO"
            ),
            "first_ordinal": (
                "Unique strict N20 admissions with a matching PROMED value on a row "
                "whose NUMPROMED parses numerically to 1"
            ),
            "ordinal_interpretation": (
                "NUMPROMED is an ordinal 1-to-N; ordinal 1 is not described as "
                "a documented primary procedure"
            ),
        },
        "catalogue_note": (
            "Reviewed CIE-9-MC catalogue labels 598X as ureteral catheterization; "
            "no specific ureteral-stent title was identified. Code 5504 is shared "
            "by two candidate groups."
        ),
        "models_fitted": 0,
        "raw_admission_identifiers_written": False,
        "raw_records_written": False,
        "linkage_repair_audit": {
            "passed": True,
            "status": "PASS",
            "population": "strict N20 admissions with ENTIDAD 1–32",
            "strict_admissions": EXPECTED_STRICT_COHORT_TOTAL,
            "source_file_year_mapping_count": partition[
                "strict_admissions_reconciled"
            ],
            "unmatched_admissions": partition["unmatched_admissions"],
            "ambiguous_admissions": partition["ambiguous_admissions"],
            "strict_cross_partition_admissions": cross_partition_total,
            "cross_year_counts": partition["cross_year_counts"],
            "legacy_rule": "procedure source file selected by admission year",
            "corrected_procedure_link_key": "ID + source_file_year",
            "exposure_date_field": "date (admission date)",
            "source_year_mapping_written": False,
            "raw_admission_identifiers_or_records_written": False,
            "statement": (
                "The 58 source-file-year/admission-year mismatches were recovered "
                "to their audited annual procedure source partitions; exposure "
                "dates remain the admission dates."
            ),
        },
    }


def write_feasibility_outputs(result: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    json_path = OUT / "mexico_procedure_feasibility.json"
    markdown_path = OUT / "mexico_procedure_feasibility.md"
    repair_audit_path = OUT / "mexico_linkage_repair_audit.json"
    json_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    markdown_path.write_text(build_report(result), encoding="utf-8")
    repair_audit_path.write_text(
        json.dumps(
            result["linkage_repair_audit"],
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    preserve_pre_repair_outputs()
    cohort, ids_by_source_year, cohort_audit = load_cohort()
    coverage, source_by_file_year, match_sets = scan_procedure_tables(
        ids_by_source_year
    )
    result = assemble_result(
        cohort, cohort_audit, coverage, source_by_file_year, match_sets
    )
    write_feasibility_outputs(result)
    json_path = OUT / "mexico_procedure_feasibility.json"
    markdown_path = OUT / "mexico_procedure_feasibility.md"

    print(
        "COMPLETE "
        f"cohort={cohort_audit['strict_n20_unique_id_year_admissions']} "
        f"linked={coverage['linked_admissions_with_any_procedure_row']} "
        f"coded={coverage['linked_admissions_with_any_nonblank_PROMED']} "
        f"source_partition_reconciled={result['linkage_repair_audit']['source_file_year_mapping_count']} "
        f"outputs={markdown_path.name},{json_path.name}",
        flush=True,
    )


if __name__ == "__main__":
    main()
