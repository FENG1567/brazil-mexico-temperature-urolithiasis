"""Source-boundary sensitivity for Mexico admissions (2018–2022 dates).

Reads the archived 2023 DGIS discharge/procedure extracts, links source-year
2023 procedure rows in memory, and writes aggregate-only boundary outputs.
The exposure calendar remains 2018–2022.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import re
import sys
import zlib
from collections import Counter, defaultdict
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[variable] = "8"

sys.path.insert(0, "F:/analysis/climate_biliary_20261102/work/pylibs")

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import norm, t as student_t

import audit_mexico_source_partitions as partition_audit
import mexico_procedure_feasibility as feasibility
import mexico_strict_phenotype as country_inputs


BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "work" / "upgrade_20261112"
OUT = BASE / "outputs" / "upgrade_20261112"
OUT.mkdir(parents=True, exist_ok=True)

ARCHIVE_EVIDENCE = Path(
    r"D:\数据库集合\新增数据\Mexico_DGIS\_evidence\full_acquisition_2023_20261030T122231Z.json"
)
EGRESOS_2023 = Path(
    r"D:\数据库集合\新增数据\Mexico_DGIS\data_full\2023\Registros\Egresos.txt"
)
PROCEDURES_2023 = Path(
    r"D:\数据库集合\新增数据\Mexico_DGIS\data_full\2023\Registros\Procedimientos.txt"
)
BASELINE_CACHE = Path(
    r"F:\analysis\climate_biliary_20261102\work\mexico_n20_id_proced.parquet"
)
ORIGINAL_PANEL = Path(
    r"F:\analysis\climate_biliary_20261102\work\mexico_composite_state_day_panel.parquet"
)
OLD_DAILY = OUT / "mexico_specific_daily_counts.csv"
HOLIDAY_DATES = OUT / "holiday_calendar_dates.csv"
OLD_SOURCE_AUDIT = OUT / "mexico_specific_global_key_audit.json"
OLD_FEASIBILITY = OUT / "mexico_procedure_feasibility.json"
WINDOW_START = dt.date(2018, 1, 1)
WINDOW_END = dt.date(2022, 12, 31)
EXPECTED_STRICT_BASE = 32_417
EXPECTED_BROAD_BASE = 40_415
EXPECTED_STATE_DAYS = 58_432
EXPECTED_EGRESOS_COLUMNS = 49
EXPECTED_PROCEDURE_COLUMNS = 7
FOOTER = re.compile(r"^\((\d+) row\(s\) affected\)$")
ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}(?:[ T].*)?$")
TARGET_DIAGNOSES = ("N20", "N21", "N22", "N23")
GROUPS = ["nephrostomy", "ureteral_catheterization",
          "urinary_stone_fragmentation_lithotripsy"]
GROUP_CODES = {
    key: set(feasibility.CODE_GROUPS[key]["codes"])
    for key in GROUPS
}
STRICT_DAILY_COLUMNS = [
    "admissions", *GROUPS, "any_surgery",
]


def log(message: str, **values: object) -> None:
    print(message, json.dumps(values, ensure_ascii=False, sort_keys=True), flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_integrity(path: Path) -> dict:
    digest = hashlib.sha256()
    crc = 0
    byte_count = 0
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            byte_count += len(block)
            digest.update(block)
            crc = zlib.crc32(block, crc)
    return {
        "path": str(path),
        "bytes": byte_count,
        "sha256": digest.hexdigest(),
        "crc32": f"{crc & 0xFFFFFFFF:08x}",
    }


def state_key(value: object) -> str:
    value = str(value or "").strip()
    return str(int(value)) if value.isdigit() else value


def audit_signature(identifier: object, raw_date: object,
                    state: object, diagnosis: object) -> tuple[str, str, str, str]:
    return partition_audit.signature(identifier, raw_date, state, diagnosis)


def parse_iso_admission_date(raw: str) -> dt.date | None:
    raw = raw.strip()
    if not ISO_DATE.fullmatch(raw):
        return None
    try:
        return dt.date.fromisoformat(raw[:10])
    except ValueError:
        return None


def require_paths() -> None:
    paths = [
        ARCHIVE_EVIDENCE, EGRESOS_2023, PROCEDURES_2023, BASELINE_CACHE,
        ORIGINAL_PANEL, OLD_DAILY, HOLIDAY_DATES, OLD_SOURCE_AUDIT,
        OLD_FEASIBILITY,
    ]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required local boundary input is missing: " + ", ".join(missing))


def verify_archive_inputs() -> dict:
    evidence = json.loads(ARCHIVE_EVIDENCE.read_text(encoding="utf-8"))
    archive_record = next(
        (entry for entry in evidence
         if entry.get("year") == 2023 and entry.get("archive_members")),
        None,
    )
    if archive_record is None:
        raise ValueError("Archived 2023 acquisition evidence is missing its archive record")
    archive_path = Path(archive_record["path"])
    if not archive_path.is_file():
        raise FileNotFoundError("Archived 2023 ZIP is missing")
    archive_sha = sha256_file(archive_path)
    if archive_sha != archive_record.get("sha256"):
        raise ValueError("Archived 2023 ZIP hash differs from its acquisition evidence")
    members = {member["name"]: member for member in archive_record["archive_members"]}
    file_checks = []
    for member_name, path in (
        ("Registros/Egresos.txt", EGRESOS_2023),
        ("Registros/Procedimientos.txt", PROCEDURES_2023),
    ):
        expected = members.get(member_name)
        if expected is None:
            raise ValueError("Archived 2023 ZIP manifest lacks an expected source member")
        observed = file_integrity(path)
        if observed["bytes"] != int(expected["uncompressed_bytes"]):
            raise ValueError("Extracted 2023 source byte count differs from the ZIP manifest")
        if observed["crc32"].lower() != str(expected["crc32"]).lower():
            raise ValueError("Extracted 2023 source CRC differs from the ZIP manifest")
        file_checks.append({
            "archive_member": member_name,
            **observed,
            "archive_crc32": str(expected["crc32"]).lower(),
            "archive_uncompressed_bytes": int(expected["uncompressed_bytes"]),
        })
    return {
        "year": 2023,
        "archive_path": str(archive_path),
        "archive_sha256": archive_sha,
        "archive_sha256_expected": archive_record["sha256"],
        "original_url": archive_record.get("original_url"),
        "final_url": archive_record.get("final_url"),
        "provenance_url": archive_record.get("provenance_url"),
        "downloaded_at_utc": archive_record.get("downloaded_at_utc"),
        "metadata_version": "2024.05.31",
        "archive_capture_timestamp": "2025-09-15T17:03:11Z",
        "acquisition_status": archive_record.get("status"),
        "live_download_performed_for_this_analysis": False,
        "extracted_sources": file_checks,
    }


def load_baseline_signatures() -> tuple[set, set, dict, dict]:
    db = duckdb.connect()
    try:
        rows = db.execute(
            "SELECT ID, INGRE, ENTIDAD, AFECPRIN FROM read_parquet(?)",
            [str(BASELINE_CACHE)],
        ).fetchall()
    finally:
        db.close()
    broad = set()
    strict = set()
    identity_diagnoses: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    duplicate_keys = 0
    counts: Counter = Counter()
    for row in rows:
        key = audit_signature(*row)
        identifier, raw_date, state, diagnosis = key
        counts[key] += 1
        if (diagnosis.startswith(TARGET_DIAGNOSES)
                and WINDOW_START.isoformat() <= raw_date[:10] <= WINDOW_END.isoformat()
                and state.isdigit() and 1 <= int(state) <= 32):
            broad.add(key)
            identity_diagnoses[(identifier, raw_date[:10], state)].add(diagnosis)
            if diagnosis.startswith("N20"):
                strict.add(key)
    duplicate_keys = sum(value - 1 for value in counts.values() if value > 1)
    if len(broad) != EXPECTED_BROAD_BASE:
        raise ValueError(
            "Cached baseline broad N20-N23 signature count differs from the audited 40,415"
        )
    if len(strict) != EXPECTED_STRICT_BASE:
        raise ValueError(
            "Cached baseline strict N20 signature count differs from the audited 32,417"
        )
    return broad, strict, dict(identity_diagnoses), {
        "cache_path": str(BASELINE_CACHE),
        "cache_sha256": sha256_file(BASELINE_CACHE),
        "cached_rows": len(rows),
        "unique_cached_signatures": len(counts),
        "broad_n20_n23_in_window_valid_state_signatures": len(broad),
        "duplicate_exact_signature_rows": duplicate_keys,
        "strict_n20_in_window_valid_state_signatures": len(strict),
        "signature_definition": "ID + ISO INGRE date + normalized ENTIDAD + uppercase AFECPRIN",
        "identifiers_written": False,
    }


def scan_egresos() -> tuple[dict[tuple[str, str, str, str], dict], dict]:
    path = EGRESOS_2023
    seen_annual_ids: set[str] = set()
    candidates: dict[tuple[str, str, str, str], dict] = {}
    target_signatures_by_id: dict[str, set[tuple[str, str, str, str]]] = defaultdict(set)
    rows = missing_ids = duplicate_id_rows = blanks = trailers = iso_dates = 0
    malformed_iso_dates = invalid_states = target_rows = missing_target_ids = 0
    discharge_parseable = discharge_nonblank_noniso = discharge_before_ingre = 0
    discharge_conflicting_rows = 0
    discharge_conflicting_signatures = set()
    footer_total = None
    diagnosis_rows = valid_date_rows = in_window_rows = 0
    with path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="|", strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("2023 Egresos source is empty") from exc
        if len(header) != EXPECTED_EGRESOS_COLUMNS or len(set(header)) != len(header):
            raise ValueError("2023 Egresos header does not match the audited 49-column format")
        required = {"ID", "INGRE", "EGRESO", "ENTIDAD", "AFECPRIN"}
        if not required.issubset(header):
            raise ValueError("2023 Egresos source lacks required admission-key fields")
        idx = {name: header.index(name) for name in required}
        for row in reader:
            if not row or (len(row) == 1 and not row[0].strip()):
                blanks += 1
                continue
            trailer = FOOTER.fullmatch(row[0].strip()) if len(row) == 1 else None
            if trailer:
                if trailers:
                    raise ValueError("Duplicate SQL export trailer in 2023 Egresos")
                footer_total = int(trailer.group(1))
                trailers += 1
                continue
            if len(row) != len(header) or trailers:
                raise ValueError(
                    f"Unexpected 2023 Egresos row structure at line {reader.line_num}"
                )
            rows += 1
            identifier = row[idx["ID"]].strip()
            if not identifier:
                missing_ids += 1
            elif identifier in seen_annual_ids:
                duplicate_id_rows += 1
            else:
                seen_annual_ids.add(identifier)

            raw_date = row[idx["INGRE"]].strip()
            admission_date = parse_iso_admission_date(raw_date)
            iso_dates += int(bool(ISO_DATE.fullmatch(raw_date)))
            if ISO_DATE.fullmatch(raw_date) and admission_date is None:
                malformed_iso_dates += 1
            if admission_date is not None:
                valid_date_rows += 1
                if WINDOW_START <= admission_date <= WINDOW_END:
                    in_window_rows += 1
            diagnosis = row[idx["AFECPRIN"]].strip().upper()
            if not diagnosis.startswith(TARGET_DIAGNOSES):
                continue
            diagnosis_rows += 1
            if admission_date is None or not (WINDOW_START <= admission_date <= WINDOW_END):
                continue
            state_raw = row[idx["ENTIDAD"]].strip()
            if not state_raw.isdigit() or not 1 <= int(state_raw) <= 32:
                invalid_states += 1
                continue
            target_rows += 1
            if not identifier:
                missing_target_ids += 1
                continue
            state_number = int(state_raw)
            key = audit_signature(identifier, raw_date, state_raw, diagnosis)
            raw_discharge = row[idx["EGRESO"]].strip()
            discharge_date = parse_iso_admission_date(raw_discharge)
            if discharge_date is not None:
                discharge_parseable += 1
                discharge_before_ingre += int(discharge_date < admission_date)
            elif raw_discharge:
                discharge_nonblank_noniso += 1
            target_signatures_by_id[identifier].add(key)
            if key not in candidates:
                candidates[key] = {
                    "ID": identifier,
                    "source_file_year": 2023,
                    "date": admission_date,
                    "discharge_date": discharge_date,
                    "state": f"{state_number:02d}",
                    "diagnosis": diagnosis,
                    "strict_N20": int(diagnosis.startswith("N20")),
                }
            elif candidates[key]["discharge_date"] != discharge_date:
                discharge_conflicting_rows += 1
                discharge_conflicting_signatures.add(key)
    if footer_total is not None and footer_total != rows:
        raise ValueError("2023 Egresos SQL export trailer count did not reconcile")
    ambiguous_target_ids = sum(
        1 for signatures in target_signatures_by_id.values() if len(signatures) > 1
    )
    if missing_target_ids:
        raise ValueError(
            f"{missing_target_ids} in-window target admissions lack an ID for procedure linkage"
        )
    if ambiguous_target_ids:
        raise ValueError(
            f"{ambiguous_target_ids} in-window target IDs map to multiple admission signatures"
        )
    for key in discharge_conflicting_signatures:
        candidates.pop(key, None)
    if rows == 0 or len(seen_annual_ids) == 0:
        raise ValueError("2023 Egresos source contains no auditable admission rows")
    return candidates, {
        "file": str(path),
        "delimiter": "|",
        "header_columns": len(header),
        "required_fields_present": sorted(required),
        "records": rows,
        "unique_annual_IDs": len(seen_annual_ids),
        "missing_ID_records": missing_ids,
        "duplicate_rows_beyond_unique_annual_ID": duplicate_id_rows,
        "annual_ID_uniqueness_passed": duplicate_id_rows == 0 and missing_ids == 0,
        "blank_lines_excluded": blanks,
        "sql_export_trailers_excluded": trailers,
        "trailer_record_total": footer_total,
        "ISO_INGRE_records": iso_dates,
        "ISO_shaped_but_invalid_INGRE_dates": malformed_iso_dates,
        "valid_INGRE_dates": valid_date_rows,
        "records_with_INGRE_in_2018_2022": in_window_rows,
        "target_diagnosis_prefix_records_N20_N23": diagnosis_rows,
        "target_records_after_date_and_valid_state_filters": target_rows,
        "invalid_or_out_of_range_state_among_target_diagnosis_dates": invalid_states,
        "unique_target_admission_signatures": len(candidates),
        "exact_duplicate_target_rows_within_2023_source": target_rows - len(candidates),
        "target_IDs_with_conflicting_admission_signatures": ambiguous_target_ids,
        "target_admissions_with_parseable_EGRESO": discharge_parseable,
        "target_admissions_with_nonblank_nonISO_EGRESO": discharge_nonblank_noniso,
        "target_admissions_with_EGRESO_before_INGRE": discharge_before_ingre,
        "duplicate_target_signatures_with_conflicting_EGRESO_excluded": len(
            discharge_conflicting_signatures
        ),
        "source_rows_with_conflicting_EGRESO_for_duplicate_signature": discharge_conflicting_rows,
        "raw_identifiers_written": False,
        "filter": "principal AFECPRIN prefix N20-N23; ISO INGRE 2018-01-01 through 2022-12-31; ENTIDAD 1-32",
    }


def reconcile_baseline(candidates: dict, baseline_signatures: set,
                       baseline_strict_signatures: set,
                       baseline_identity_diagnoses: dict) -> tuple[dict, dict, dict]:
    accepted = {}
    exact_duplicates = 0
    identity_conflicts = 0
    baseline_internal_identity_conflicts = sum(
        len(diagnoses) > 1 for diagnoses in baseline_identity_diagnoses.values()
    )
    exact_broad_duplicates = 0
    exact_strict_duplicates = 0
    for signature, record in candidates.items():
        identifier, raw_date, state, diagnosis = signature
        identity = (identifier, raw_date[:10], state)
        prior_diagnoses = baseline_identity_diagnoses.get(identity, set())
        if signature in baseline_signatures and prior_diagnoses == {diagnosis}:
            exact_duplicates += 1
            exact_broad_duplicates += 1
            exact_strict_duplicates += int(signature in baseline_strict_signatures)
            continue
        if prior_diagnoses:
            identity_conflicts += 1
            continue
        accepted[signature] = record
    if baseline_internal_identity_conflicts:
        raise ValueError(
            f"{baseline_internal_identity_conflicts} baseline admissions have conflicting diagnoses"
        )
    if identity_conflicts:
        raise ValueError(
            f"{identity_conflicts} source-2023 admissions conflict with baseline diagnosis at ID/date/state"
        )
    accepted_by_year: dict[str, Counter] = defaultdict(Counter)
    discharge_year_counts: Counter = Counter()
    discharge_evaluated = discharge_passed = discharge_failed = discharge_unavailable = 0
    for record in accepted.values():
        admission_year = str(record["date"].year)
        accepted_by_year[admission_year]["broad_N20_N23"] += 1
        accepted_by_year[admission_year]["strict_N20"] += int(record["strict_N20"])
        discharge_date = record["discharge_date"]
        if discharge_date is None:
            discharge_unavailable += 1
            discharge_year_counts["missing_or_nonISO"] += 1
        else:
            discharge_year_counts[str(discharge_date.year)] += 1
            discharge_evaluated += 1
            if record["date"] <= discharge_date:
                discharge_passed += 1
            else:
                discharge_failed += 1
    duplicate_audit = {
        "compared_against": "audited N20-N23 baseline identifier cache",
        "baseline_signature_match_duplicates_excluded": exact_duplicates,
        "baseline_broad_signature_duplicates_excluded": exact_broad_duplicates,
        "baseline_strict_N20_signature_duplicates_excluded": exact_strict_duplicates,
        "same_ID_date_state_but_diagnosis_conflicts": identity_conflicts,
        "new_unique_admissions_retained": len(accepted),
        "duplicate_decision": (
            "Exact same admission signatures are excluded from additions. A conflicting "
            "principal diagnosis at an existing ID/date/state stops the pipeline for review."
        ),
        "identifiers_written": False,
    }
    discharge_audit = {
        "accepted_admissions_by_admission_year": {
            year: dict(counts) for year, counts in sorted(accepted_by_year.items())
        },
        "accepted_admissions_by_discharge_year": dict(sorted(discharge_year_counts.items())),
        "INGRE_le_EGRESO_validation": {
            "assessable": discharge_evaluated,
            "passed": discharge_passed,
            "failed": discharge_failed,
            "missing_or_nonISO_EGRESO": discharge_unavailable,
            "all_assessable_admissions_passed": discharge_failed == 0,
        },
    }
    return accepted, duplicate_audit, discharge_audit


def scan_procedures(target_ids: set[str]) -> tuple[dict, dict]:
    path = PROCEDURES_2023
    rows = linked_rows = missing_code_rows = blanks = trailers = 0
    footer_total = None
    seen_any_row = set()
    seen_nonblank_code = set()
    seen_blank_code = set()
    q_ids = set()
    group_ids = {key: set() for key in GROUPS}
    with path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="|", strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("2023 Procedimientos source is empty") from exc
        if len(header) != EXPECTED_PROCEDURE_COLUMNS or len(set(header)) != len(header):
            raise ValueError("2023 Procedimientos header does not match the audited 7-column format")
        required = set(feasibility.REQUIRED_PROCEDURE_COLUMNS)
        if not required.issubset(header):
            raise ValueError("2023 Procedimientos source lacks required linkage fields")
        idx = {name: header.index(name) for name in required}
        for row in reader:
            if not row or (len(row) == 1 and not row[0].strip()):
                blanks += 1
                continue
            trailer = FOOTER.fullmatch(row[0].strip()) if len(row) == 1 else None
            if trailer:
                if trailers:
                    raise ValueError("Duplicate SQL export trailer in 2023 Procedimientos")
                footer_total = int(trailer.group(1))
                trailers += 1
                continue
            if len(row) != len(header) or trailers:
                raise ValueError(
                    f"Unexpected 2023 Procedimientos row structure at line {reader.line_num}"
                )
            rows += 1
            identifier = row[idx["ID"]].strip()
            if identifier not in target_ids:
                continue
            linked_rows += 1
            seen_any_row.add(identifier)
            if row[idx["TIPO"]].strip().upper() == "Q":
                q_ids.add(identifier)
            code = row[idx["PROMED"]].strip().upper()
            if not code:
                missing_code_rows += 1
                seen_blank_code.add(identifier)
                continue
            seen_nonblank_code.add(identifier)
            for group_key, codes in GROUP_CODES.items():
                if code in codes:
                    group_ids[group_key].add(identifier)
    if footer_total is not None and footer_total != rows:
        raise ValueError("2023 Procedimientos SQL export trailer count did not reconcile")
    if not q_ids.issubset(target_ids):
        raise ValueError("TIPO=Q procedure links exceed the selected boundary cohort")
    if any(not values.issubset(target_ids) for values in group_ids.values()):
        raise ValueError("Procedure-code groups exceed the selected boundary cohort")
    summary = {
        "file": str(path),
        "delimiter": "|",
        "header_columns": len(header),
        "required_fields_present": sorted(required),
        "procedure_rows_parsed": rows,
        "procedure_rows_linked_to_boundary_cohort": linked_rows,
        "boundary_admissions_with_any_procedure_row": len(seen_any_row),
        "boundary_admissions_with_any_nonblank_PROMED": len(seen_nonblank_code),
        "boundary_admissions_with_no_nonblank_PROMED": len(target_ids - seen_nonblank_code),
        "boundary_admissions_with_any_blank_PROMED_row": len(seen_blank_code),
        "linked_procedure_rows_missing_PROMED": missing_code_rows,
        "boundary_admissions_with_any_TIPO_Q_row": len(q_ids),
        "blank_lines_excluded": blanks,
        "sql_export_trailers_excluded": trailers,
        "trailer_record_total": footer_total,
        "raw_identifiers_written": False,
    }
    return {
        "TIPO_Q": q_ids,
        **group_ids,
        "_any_row": seen_any_row,
        "_any_nonblank_PROMED": seen_nonblank_code,
    }, summary


def sha_for_existing_source_audit() -> dict:
    audit = json.loads(OLD_SOURCE_AUDIT.read_text(encoding="utf-8"))
    if len(audit.get("years", [])) != 5 or not audit.get("global_annual_ID_uniqueness_passed"):
        raise ValueError("Existing 2018-2022 annual-source audit is missing or failed")
    return {
        "source_scope": "annual egresos source files 2018-2022; admission dates filtered to 2018-2022",
        "annual_source_key_audit": str(OLD_SOURCE_AUDIT),
        "annual_source_key_audit_sha256": sha256_file(OLD_SOURCE_AUDIT),
        "source_partition_audit": str(OLD_FEASIBILITY),
        "source_partition_audit_sha256": sha256_file(OLD_FEASIBILITY),
        "annual_files": audit["sources"],
        "annual_id_uniqueness_passed": audit["global_annual_ID_uniqueness_passed"],
        "strict_n20_partition_reconciliation": audit["source_partition_audit"],
    }


def make_boundary_daily(accepted: dict, procedure_flags: dict,
                        daily_base: pd.DataFrame, panel: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    additions = []
    q_ids = procedure_flags["TIPO_Q"]
    for signature, record in accepted.items():
        identifier = record["ID"]
        is_strict = int(record["strict_N20"])
        additions.append({
            "state": record["state"],
            "date": pd.Timestamp(record["date"]),
            "added_2023_broad": 1,
            "added_2023_strict": is_strict,
            "added_2023_nephrostomy": int(is_strict and identifier in procedure_flags["nephrostomy"]),
            "added_2023_ureteral_catheterization": int(
                is_strict and identifier in procedure_flags["ureteral_catheterization"]
            ),
            "added_2023_urinary_stone_fragmentation_lithotripsy": int(
                is_strict and identifier in procedure_flags["urinary_stone_fragmentation_lithotripsy"]
            ),
            "added_2023_TIPO_Q": int(is_strict and identifier in q_ids),
        })
    add_cols = [
        "added_2023_broad", "added_2023_strict", "added_2023_nephrostomy",
        "added_2023_ureteral_catheterization",
        "added_2023_urinary_stone_fragmentation_lithotripsy", "added_2023_TIPO_Q",
    ]
    added = pd.DataFrame(additions)
    if added.empty:
        added_daily = pd.DataFrame(columns=["state", "date", *add_cols])
    else:
        added_daily = added.groupby(["state", "date"], as_index=False, observed=True)[add_cols].sum()
    base_daily = daily_base.copy()
    base_daily["date"] = pd.to_datetime(base_daily.date).dt.normalize()
    base_daily["state"] = base_daily.state.astype(str).str.zfill(2)
    panel = panel.copy()
    panel["date"] = pd.to_datetime(panel.date).dt.normalize()
    panel["state"] = panel.state.astype(str).str.zfill(2)
    full = panel[["state", "date", "broad_N20_N23", "strict_N20", "stratum"]].merge(
        base_daily[["state", "date", *STRICT_DAILY_COLUMNS]],
        on=["state", "date"], how="left", validate="one_to_one",
    )
    for column in STRICT_DAILY_COLUMNS:
        full[column] = full[column].fillna(0).astype("int64")
    if len(full) != EXPECTED_STATE_DAYS or full.duplicated(["state", "date"]).any():
        raise ValueError("Base Mexico boundary calendar differs from 32 states x 2018-2022")
    if int(full.admissions.sum()) != EXPECTED_STRICT_BASE:
        raise ValueError("Corrected strict daily base does not total 32,417")
    if not full.admissions.eq(full.strict_N20).all():
        raise ValueError("Corrected strict daily base differs from audited strict phenotype")
    if not full.any_surgery.between(0, full.admissions).all():
        raise ValueError("Base TIPO=Q counts do not partition strict admissions")
    base_counts = {
        "broad_N20_N23": int(full.broad_N20_N23.sum()),
        "strict_N20": int(full.admissions.sum()),
        **{column: int(full[column].sum()) for column in GROUPS},
        "TIPO_Q": int(full.any_surgery.sum()),
        "no_TIPO_Q": int((full.admissions - full.any_surgery).sum()),
    }
    if base_counts["broad_N20_N23"] != EXPECTED_BROAD_BASE:
        raise ValueError("Base broad Mexico counts do not total 40,415")
    if base_counts["TIPO_Q"] != 10_488 or base_counts["no_TIPO_Q"] != 21_929:
        raise ValueError("Base corrected TIPO=Q and no-Q counts differ from the audited totals")

    full = full.merge(added_daily, on=["state", "date"], how="left", validate="one_to_one")
    for column in add_cols:
        full[column] = full[column].fillna(0).astype("int64")
    full["broad_N20_N23"] += full["added_2023_broad"]
    full["admissions"] += full["added_2023_strict"]
    for group in GROUPS:
        full[group] += full[f"added_2023_{group}"]
    full["any_surgery"] += full["added_2023_TIPO_Q"]
    full["no_TIPO_Q"] = full.admissions - full.any_surgery
    if not full.any_surgery.between(0, full.admissions).all():
        raise ValueError("Augmented TIPO=Q counts do not partition strict admissions")
    for group in GROUPS:
        if not full[group].between(0, full.admissions).all():
            raise ValueError("Augmented procedure group exceeds strict admissions")
    full = full.rename(columns={
        "any_surgery": "TIPO_Q",
    })
    audit = {
        "base_counts_2018_2022_source_years": base_counts,
        "accepted_source_2023_additions": {
            "broad_N20_N23": int(full.added_2023_broad.sum()),
            "strict_N20": int(full.added_2023_strict.sum()),
            **{group: int(full[f"added_2023_{group}"].sum()) for group in GROUPS},
            "TIPO_Q": int(full.added_2023_TIPO_Q.sum()),
            "no_TIPO_Q": int(full.added_2023_strict.sum() - full.added_2023_TIPO_Q.sum()),
        },
        "augmented_counts_source_years_2018_2023": {
            "broad_N20_N23": int(full.broad_N20_N23.sum()),
            "strict_N20": int(full.admissions.sum()),
            **{group: int(full[group].sum()) for group in GROUPS},
            "TIPO_Q": int(full.TIPO_Q.sum()),
            "no_TIPO_Q": int(full.no_TIPO_Q.sum()),
        },
        "state_day_rows": len(full),
        "state_count": int(full.state.nunique()),
        "admission_date_window": "2018-01-01 through 2022-12-31",
        "exposure_date_window": "2018-01-01 through 2022-12-31",
    }
    return full, audit


def add_lags(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.sort_values(["state", "date"]).copy()
    for lag in range(8):
        for column in columns:
            prefix = "t" if column == "temperature" else "f"
            out[f"{prefix}{lag}"] = out.groupby("state", sort=False)[column].shift(lag)
    return out


def fit_two_way_lag_model(frame: pd.DataFrame, outcome: str, outcome_label: str,
                          full_events: int, models: dict, result_rows: list) -> None:
    terms = [f"t{lag}" for lag in range(8)]
    model_frame = frame[["state", "date", "date_id", "stratum", *terms, "holiday", outcome]].dropna(
        subset=terms
    ).copy()
    fe_sum = model_frame.groupby("stratum", observed=True)[outcome].sum()
    fe_size = model_frame.groupby("stratum", observed=True).size()
    zero_levels = fe_sum.index[fe_sum.eq(0)]
    expected_dropped = int(fe_size.reindex(zero_levels).sum())
    formula = outcome + " ~ " + " + ".join(terms + ["holiday"]) + " | stratum"
    model = pf.fepois(formula, data=model_frame, vcov={"CRV1": "state + date_id"})
    names = [str(name) for name in model.coef().index]
    beta = model.coef().to_numpy(dtype=float)
    covariance = np.asarray(model._vcov, dtype=float)
    if covariance.shape != (len(names), len(names)) or not np.isfinite(covariance).all():
        raise ValueError("Two-way boundary model returned an invalid covariance matrix")
    clusters_state = int(model_frame.state.nunique())
    clusters_date = int(model_frame.date_id.nunique())
    df = min(clusters_state, clusters_date) - 1
    if df <= 0:
        raise ValueError("Boundary model has insufficient two-way clusters")
    model_n = int(model._N)
    if model_n != len(model_frame) - expected_dropped:
        raise ValueError("Boundary model dropped rows beyond zero-event strata")
    terms_present = set(names)
    if not set(terms + ["holiday"]).issubset(terms_present):
        raise ValueError("Boundary model is missing a predeclared term")
    z = np.array([float(term in terms) for term in names])
    beta_sum = float(z @ beta)
    variance = float(z @ covariance @ z)
    if variance <= 0:
        raise ValueError("Boundary cumulative temperature contrast has nonpositive variance")
    se = float(np.sqrt(variance))
    z_score = beta_sum / se
    p_value = float(2 * student_t.sf(abs(z_score), df))
    row = {
        "model": outcome_label,
        "outcome": outcome,
        "contrast": "cumulative_temperature_lag_0_7_per_C",
        "RR": float(np.exp(beta_sum)),
        "CI_low": float(np.exp(beta_sum - student_t.ppf(0.975, df) * se)),
        "CI_high": float(np.exp(beta_sum + student_t.ppf(0.975, df) * se)),
        "p": p_value,
        "beta": beta_sum,
        "se_full_covariance": se,
        "q_BH_six_procedure_contrasts": None,
        "events_full_window": int(full_events),
        "events_after_lag_window": int(model_frame[outcome].sum()),
        "input_rows_after_lag": len(model_frame),
        "model_rows": model_n,
        "state_clusters": clusters_state,
        "date_clusters": clusters_date,
        "t_df": df,
        "formula": formula,
        "fixed_effects": "state x calendar year-month x day-of-week",
        "cluster": "state + calendar date CRV1",
        "inference": "t(min(G_state,G_date)-1)",
    }
    result_rows.append(row)
    models[outcome_label] = {
        "outcome": outcome,
        "terms": names,
        "coefficients": beta.tolist(),
        "vcov": covariance.tolist(),
        "formula": formula,
        "model_rows": model_n,
        "input_rows_after_lag": len(model_frame),
        "zero_event_strata": int(len(zero_levels)),
        "zero_event_rows_expected_dropped": expected_dropped,
        "events_full_window": int(full_events),
        "events_after_lag_window": int(model_frame[outcome].sum()),
        "state_clusters": clusters_state,
        "date_clusters": clusters_date,
        "t_df": df,
        "contrast": {
            "name": "cumulative_temperature_lag_0_7_per_C",
            "terms": terms,
            "beta": beta_sum,
            "variance": variance,
            "se": se,
            "RR": row["RR"],
            "CI_low": row["CI_low"],
            "CI_high": row["CI_high"],
            "p": p_value,
        },
    }
    log("MODEL_COMPLETE", model=outcome_label, events=full_events,
        rows=model_n, state_clusters=clusters_state, date_clusters=clusters_date)


def fit_procedure_lag_models(frame: pd.DataFrame, models: dict,
                             result_rows: list) -> None:
    lag_terms = [f"{prefix}{lag}" for lag in range(8) for prefix in ("t", "f")]
    for outcome in GROUPS:
        lag_valid = frame[[
            "state", "date", "date_id", "stratum", *lag_terms, "holiday", outcome
        ]].dropna(subset=lag_terms).copy()
        informative = lag_valid.groupby("stratum", observed=True)[outcome].transform("sum") > 0
        model_frame = lag_valid.loc[informative].copy()
        if model_frame.empty:
            raise ValueError("A predeclared boundary procedure model has no informative rows")
        formula = outcome + " ~ " + " + ".join(lag_terms + ["holiday"]) + " | stratum"
        model = pf.fepois(formula, data=model_frame, vcov={"CRV1": "state + date_id"})
        names = [str(name) for name in model.coef().index]
        beta = model.coef().to_numpy(dtype=float)
        covariance = np.asarray(model._vcov, dtype=float)
        if model._N != len(model_frame) or not np.isfinite(covariance).all():
            raise ValueError("Procedure model row count or covariance did not reconcile")
        state_clusters = int(model_frame.state.nunique())
        date_clusters = int(model_frame.date_id.nunique())
        df = min(state_clusters, date_clusters) - 1
        rows_for_model = []
        contrast_data = {}
        for prefix, label, scale_name in (
            ("t", "cumulative_temperature_lag_0_7_per_C", "temperature"),
            ("f", "cumulative_fire_pm25_lag_0_7_per_10ug_m3", "fire_pm25"),
        ):
            selected = {f"{prefix}{lag}" for lag in range(8)}
            if not selected.issubset(names):
                raise ValueError("Predeclared procedure distributed-lag terms are missing")
            weights = np.array([float(term in selected) for term in names])
            beta_sum = float(weights @ beta)
            variance = float(weights @ covariance @ weights)
            if variance <= 0:
                raise ValueError("Procedure distributed-lag contrast has nonpositive variance")
            se = float(np.sqrt(variance))
            p_value = float(2 * student_t.sf(abs(beta_sum / se), df))
            row = {
                "model": f"strict_N20_{outcome}",
                "outcome": outcome,
                "contrast": label,
                "RR": float(np.exp(beta_sum)),
                "CI_low": float(np.exp(beta_sum - student_t.ppf(0.975, df) * se)),
                "CI_high": float(np.exp(beta_sum + student_t.ppf(0.975, df) * se)),
                "p": p_value,
                "beta": beta_sum,
                "se_full_covariance": se,
                "q_BH_six_procedure_contrasts": None,
                "events_full_window": int(frame[outcome].sum()),
                "events_after_lag_window": int(lag_valid[outcome].sum()),
                "events_model_input": int(model_frame[outcome].sum()),
                "input_rows_after_lag": len(lag_valid),
                "informative_rows_after_lag": len(model_frame),
                "model_rows": int(model._N),
                "state_clusters": state_clusters,
                "date_clusters": date_clusters,
                "t_df": df,
                "formula": formula,
                "fixed_effects": "state x calendar year-month x day-of-week",
                "cluster": "state + calendar date CRV1",
                "inference": "t(min(G_state,G_date)-1)",
            }
            rows_for_model.append(row)
            result_rows.append(row)
            contrast_data[scale_name] = {
                "name": label,
                "terms": sorted(selected),
                "beta": beta_sum,
                "variance": variance,
                "se": se,
                "RR": row["RR"],
                "CI_low": row["CI_low"],
                "CI_high": row["CI_high"],
                "p": p_value,
            }
        models[f"strict_N20_{outcome}_temperature_fire_lag07_holiday"] = {
            "outcome": outcome,
            "terms": names,
            "coefficients": beta.tolist(),
            "vcov": covariance.tolist(),
            "formula": formula,
            "model_rows": int(model._N),
            "input_rows_after_lag": len(lag_valid),
            "informative_rows_after_lag": len(model_frame),
            "events_full_window": int(frame[outcome].sum()),
            "events_after_lag_window": int(lag_valid[outcome].sum()),
            "events_model_input": int(model_frame[outcome].sum()),
            "state_clusters": state_clusters,
            "date_clusters": date_clusters,
            "t_df": df,
            "contrasts": contrast_data,
        }
        log("MODEL_COMPLETE", model=f"strict_N20_{outcome}",
            events=int(frame[outcome].sum()), rows=int(model._N),
            state_clusters=state_clusters, date_clusters=date_clusters)


def fit_same_day_models(daily: pd.DataFrame, models: dict,
                        result_rows: list) -> None:
    db = duckdb.connect()
    try:
        panel = db.execute(
            "SELECT state, date, tmean, fire_pm25, stratum "
            "FROM read_parquet(?) WHERE date >= DATE '2018-01-01' "
            "AND date < DATE '2023-01-01'",
            [str(ORIGINAL_PANEL)],
        ).df()
    finally:
        db.close()
    panel["state"] = panel.state.astype(str).str.zfill(2)
    panel["date"] = pd.to_datetime(panel.date).dt.normalize()
    if len(panel) != EXPECTED_STATE_DAYS or panel.duplicated(["state", "date"]).any():
        raise ValueError("Original Mexico surgery model panel calendar changed")
    panel = panel.merge(
        daily[["state", "date", "TIPO_Q", "no_TIPO_Q"]],
        on=["state", "date"], how="left", validate="one_to_one",
    )
    panel[["TIPO_Q", "no_TIPO_Q"]] = panel[["TIPO_Q", "no_TIPO_Q"]].fillna(0).astype("int64")
    panel["fire10"] = panel.fire_pm25 / 10.0
    if panel[["tmean", "fire10", "stratum"]].isna().any().any():
        raise ValueError("Original same-day surgery covariate coverage changed")
    for outcome in ("TIPO_Q", "no_TIPO_Q"):
        formula = f"{outcome} ~ tmean + fire10 | stratum"
        fe_sum = panel.groupby("stratum", observed=True)[outcome].sum()
        fe_size = panel.groupby("stratum", observed=True).size()
        zero_levels = fe_sum.index[fe_sum.eq(0)]
        expected_dropped = int(fe_size.reindex(zero_levels).sum())
        model = pf.fepois(formula, data=panel, vcov={"CRV1": "state"})
        names = [str(name) for name in model.coef().index]
        beta = model.coef().to_numpy(dtype=float)
        covariance = np.asarray(model._vcov, dtype=float)
        if (model._N != len(panel) - expected_dropped
                or not np.isfinite(covariance).all()):
            raise ValueError("Same-day boundary model rows or covariance did not reconcile")
        if not {"tmean", "fire10"}.issubset(names):
            raise ValueError("Same-day surgery model is missing a predeclared term")
        for index, term in enumerate(names):
            coefficient = float(beta[index])
            se = float(np.sqrt(covariance[index, index]))
            if se <= 0:
                raise ValueError("Same-day surgery model has nonpositive standard error")
            result_rows.append({
                "model": f"strict_N20_{outcome}_same_day",
                "outcome": outcome,
                "contrast": (
                    "same_day_temperature_per_C" if term == "tmean"
                    else "same_day_fire_pm25_per_10ug_m3"
                ),
                "term": term,
                "RR": float(np.exp(coefficient)),
                "CI_low": float(np.exp(coefficient - 1.96 * se)),
                "CI_high": float(np.exp(coefficient + 1.96 * se)),
                "p": float(2 * norm.sf(abs(coefficient / se))),
                "beta": coefficient,
                "se_full_covariance": se,
                "q_BH_six_procedure_contrasts": None,
                "events_full_window": int(panel[outcome].sum()),
                "events_after_lag_window": int(panel[outcome].sum()),
                "input_rows_after_lag": len(panel),
                "model_rows": int(model._N),
                "input_rows": len(panel),
                "zero_event_strata": int(len(zero_levels)),
                "zero_event_rows_expected_dropped": expected_dropped,
                "state_clusters": int(panel.state.nunique()),
                "date_clusters": int(panel.date.nunique()),
                "t_df": None,
                "formula": formula,
                "fixed_effects": "state x calendar year-month x day-of-week stratum",
                "cluster": "state CRV1",
                "inference": "1.96 normal approximation",
            })
        models[f"strict_N20_{outcome}_same_day"] = {
            "outcome": outcome,
            "terms": names,
            "coefficients": beta.tolist(),
            "vcov": covariance.tolist(),
            "formula": formula,
            "model_rows": int(model._N),
            "input_rows": len(panel),
            "zero_event_strata": int(len(zero_levels)),
            "zero_event_rows_expected_dropped": expected_dropped,
            "events_full_window": int(panel[outcome].sum()),
            "state_clusters": int(panel.state.nunique()),
            "fixed_effects": "state x calendar year-month x day-of-week stratum",
            "cluster": "state CRV1",
            "inference": "1.96 normal approximation",
        }
        log("MODEL_COMPLETE", model=f"strict_N20_{outcome}_same_day",
            events=int(panel[outcome].sum()), rows=int(model._N),
            state_clusters=int(panel.state.nunique()))


def bh_adjust(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    q = np.empty(len(values), dtype=float)
    q[order] = np.minimum(adjusted, 1.0)
    return q


def write_outputs(source_provenance: dict, archive_provenance: dict,
                  base_cache_audit: dict, egresos_audit: dict,
                  duplicate_audit: dict, procedure_audit: dict,
                  daily_audit: dict, daily: pd.DataFrame,
                  results: list, models: dict) -> None:
    if len(models) != 7:
        raise ValueError(f"Exactly seven predeclared model fits were required, found {len(models)}")
    if len(results) != 12:
        raise ValueError(f"Expected twelve predeclared model contrasts, found {len(results)}")
    proc_mask = [
        row["model"].startswith("strict_N20_")
        and row["contrast"] in {
            "cumulative_temperature_lag_0_7_per_C",
            "cumulative_fire_pm25_lag_0_7_per_10ug_m3",
        }
        and row["outcome"] in GROUPS
        for row in results
    ]
    proc_rows = [row for row, keep in zip(results, proc_mask) if keep]
    if len(proc_rows) != 6:
        raise ValueError("The fixed six-contrast procedure multiplicity family changed")
    q_values = bh_adjust(np.array([row["p"] for row in proc_rows], dtype=float))
    for row, q_value in zip(proc_rows, q_values):
        row["q_BH_six_procedure_contrasts"] = float(q_value)
    results_path = OUT / "mexico_boundary_model_results.csv"
    covariance_path = OUT / "mexico_boundary_full_covariance.json"
    audit_path = OUT / "mexico_boundary_audit.json"
    report_path = OUT / "mexico_boundary_report_zh.md"
    daily_path = OUT / "mexico_boundary_daily_counts.csv"
    daily_columns = [
        "state", "date", "broad_N20_N23", "admissions", "nephrostomy",
        "ureteral_catheterization", "urinary_stone_fragmentation_lithotripsy",
        "TIPO_Q", "no_TIPO_Q", "added_2023_broad", "added_2023_strict",
        "added_2023_nephrostomy", "added_2023_ureteral_catheterization",
        "added_2023_urinary_stone_fragmentation_lithotripsy",
        "added_2023_TIPO_Q",
    ]
    daily = daily[daily_columns].sort_values(["state", "date"])
    daily["date"] = pd.to_datetime(daily.date).dt.strftime("%Y-%m-%d")
    daily.to_csv(daily_path, index=False, encoding="utf-8-sig")
    pd.DataFrame(results).to_csv(results_path, index=False, encoding="utf-8-sig")
    covariance_payload = {
        "source_years": [2018, 2019, 2020, 2021, 2022, 2023],
        "admission_date_window": "2018-01-01 through 2022-12-31",
        "exposure_window": "2018-01-01 through 2022-12-31",
        "predeclared_model_fits": 7,
        "models": models,
        "procedure_contrast_multiplicity": {
            "family_size": 6,
            "adjustment": "Benjamini-Hochberg",
            "contrasts": [
                f"{row['outcome']}:{row['contrast']}" for row in proc_rows
            ],
        },
    }
    covariance_path.write_text(
        json.dumps(covariance_payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    audit = {
        "status": "complete",
        "analysis": "Mexico source-boundary sensitivity",
        "source_periods": [2018, 2019, 2020, 2021, 2022, 2023],
        "admission_date_window": "2018-01-01 through 2022-12-31",
        "exposure_window": "2018-01-01 through 2022-12-31; no 2023 exposure dates",
        "principal_diagnosis_scope": "N20-N23; strict N20 outcomes separately identified",
        "annual_source_2023": archive_provenance,
        "prior_source_years_2018_2022": source_provenance,
        "baseline_ID_cache_audit": base_cache_audit,
        "egresos_2023_audit": egresos_audit,
        "cross_source_duplicate_audit": duplicate_audit,
        "procedures_2023_linkage_audit": procedure_audit,
        "daily_count_reconciliation": daily_audit,
        "model_fit_count": len(models),
        "model_keys": list(models),
        "model_specs": [
            "broad N20-N23 ERA5 temperature lags 0-7 + 38-date Mexico holiday indicator; state x year-month x weekday FE; state/date CRV1; t(min(Gstate,Gdate)-1)",
            "strict N20 same as broad model",
            "three strict N20 procedure-group models with ERA5 temperature and fire PM2.5 lags 0-7 + holiday; same FE, two-way clustering and t approximation; six BH-adjusted contrasts",
            "TIPO=Q and no recorded TIPO=Q same-day tmean + fire PM2.5/10; original 2018-2022 state-day panel, state CRV1, 1.96 normal intervals",
        ],
        "TIPO_Q_definition": (
            "Strict N20 admission with at least one linked TIPO=Q row in its audited "
            "source-file-year procedure partition; the 2023 increment is separately "
            "counted from the 2023 source-year procedure table."
        ),
        "no_TIPO_Q_definition": (
            "Strict N20 admission with no recorded TIPO=Q row in its audited "
            "source-file-year procedure partition."
        ),
        "procedure_group_codes": {
            group: sorted(GROUP_CODES[group]) for group in GROUPS
        },
        "procedure_code_semantics": {
            group: feasibility.CODE_GROUPS[group]["label"] for group in GROUPS
        },
        "procedure_flags_use_any_PROMED_row_regardless_of_TIPO": True,
        "operation_date_available": False,
        "raw_IDs_or_row_level_records_written": False,
        "municipality_flows_recomputed": False,
        "late_discharge_censoring": (
            "Admissions discharged after the available 2023 source year may still be absent."
        ),
        "source_version_difference": (
            "The archived 2023 SSA package carries metadata version 2024.05.31 and differs "
            "from the previously used annual source partitions; this is a source-boundary "
            "sensitivity, not proof of universal hospitalization completeness."
        ),
        "outputs": {
            path.name: sha256_file(path)
            for path in (results_path, covariance_path, daily_path)
        },
    }
    validation_checks = {
        "2023_annual_ID_uniqueness": bool(
            egresos_audit["annual_ID_uniqueness_passed"]
        ),
        "all_accepted_INGRE_dates_not_after_EGRESO": bool(
            egresos_audit["INGRE_le_EGRESO_validation"]["all_assessable_admissions_passed"]
            and egresos_audit["INGRE_le_EGRESO_validation"]["failed"] == 0
        ),
        "no_missing_or_nonISO_EGRESO_among_accepted": (
            egresos_audit["INGRE_le_EGRESO_validation"]["missing_or_nonISO_EGRESO"] == 0
        ),
        "no_ambiguous_target_IDs": (
            egresos_audit["target_IDs_with_conflicting_admission_signatures"] == 0
        ),
        "no_unresolved_cross_source_diagnosis_conflicts": (
            duplicate_audit["same_ID_date_state_but_diagnosis_conflicts"] == 0
        ),
        "baseline_broad_and_strict_signatures_reconciled": (
            base_cache_audit["broad_n20_n23_in_window_valid_state_signatures"]
            == EXPECTED_BROAD_BASE
            and base_cache_audit["strict_n20_in_window_valid_state_signatures"]
            == EXPECTED_STRICT_BASE
        ),
        "boundary_daily_calendar_reconciled": (
            daily_audit["state_day_rows"] == EXPECTED_STATE_DAYS
            and daily_audit["state_count"] == 32
            and daily_audit["exposure_date_window"] == "2018-01-01 through 2022-12-31"
        ),
        "exactly_seven_predeclared_model_fits": len(models) == 7,
        "all_twelve_predeclared_contrasts": len(results) == 12,
    }
    audit["validation_checks"] = validation_checks
    audit["passed"] = all(validation_checks.values())
    if not audit["passed"]:
        raise ValueError("At least one source, date, count, or model validation failed")
    audit_path.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_report(report_path, audit, results)
    audit["outputs"][report_path.name] = sha256_file(report_path)
    audit_path.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_report(path: Path, audit: dict, results: list) -> None:
    daily = audit["daily_count_reconciliation"]
    additions = daily["accepted_source_2023_additions"]
    totals = daily["augmented_counts_source_years_2018_2023"]
    base = daily["base_counts_2018_2022_source_years"]
    new_by_year = audit["egresos_2023_audit"].get("accepted_admissions_by_year", {})
    new_by_discharge_year = audit["egresos_2023_audit"].get(
        "accepted_admissions_by_discharge_year", {}
    )
    discharge_validation = audit["egresos_2023_audit"].get(
        "INGRE_le_EGRESO_validation", {}
    )
    duplicate = audit["cross_source_duplicate_audit"]
    rows = []
    for result in results:
        q = result.get("q_BH_six_procedure_contrasts")
        q_text = "" if q is None else f"{q:.4g}"
        rows.append(
            f"| {result['model']} | {result['contrast']} | {result['RR']:.4f} "
            f"| {result['CI_low']:.4f}–{result['CI_high']:.4f} | {result['p']:.4g} | {q_text} | "
            f"{result['events_full_window']:,} | {result['model_rows']:,} |"
        )
    count_rows = []
    for year in sorted(new_by_year, key=int):
        item = new_by_year[year]
        count_rows.append(
            f"| {year} | {item['broad_N20_N23']:,} | {item['strict_N20']:,} |"
        )
    count_table = "\n".join(count_rows) if count_rows else "| 2018–2022 | 0 | 0 |"
    discharge_rows = [
        f"| {year} | {count:,} |"
        for year, count in sorted(new_by_discharge_year.items())
    ]
    discharge_table = "\n".join(discharge_rows) if discharge_rows else "| 无法判定 | 0 |"
    model_table = "\n".join(rows)
    report = f"""# 墨西哥来源边界敏感性

本分析补入官方SSA年度出院来源文件2018–2023中，入院日期落在2018-01-01至2022-12-31的N20–N23记录。2023年来源来自历史存档的官方包（元数据版本2024.05.31）；本次未重新下载。按入院日期匹配暴露，温度、火源、节假日及所有模型日期均仍限于2018–2022。

## 新增记录与去重

| 入院年份 | 新增N20–N23 | 新增严格N20 |
|---|---:|---:|
{count_table}

2023来源按49列出院表和7列操作表严格解析；唯一ID检查为{audit['egresos_2023_audit']['unique_annual_IDs']:,}个年度ID、{audit['egresos_2023_audit']['duplicate_rows_beyond_unique_annual_ID']:,}条超出唯一ID的记录。该来源筛选到{audit['egresos_2023_audit']['unique_target_admission_signatures']:,}条不同目标入院签名。

与既有N20–N23缓存的ID、ISO入院日期、州和主要诊断签名比较后，排除完全重复入院{duplicate['baseline_signature_match_duplicates_excluded']:,}例，其中严格N20重复{duplicate['baseline_strict_N20_signature_duplicates_excluded']:,}例。ID/日期/州相同而主要诊断冲突时，脚本会停止供人工核对；本次完成时此类冲突为{duplicate['same_ID_date_state_but_diagnosis_conflicts']:,}例。保留的来源2023新增为N20–N23 {additions['broad_N20_N23']:,}例，其中严格N20 {additions['strict_N20']:,}例。

按可解析的EGRESO日期统计新增记录：

| EGRESO年份 | 新增住院 |
|---|---:|
{discharge_table}

在{discharge_validation.get('assessable', 0):,}条可解析记录中，INGRE不晚于EGRESO者为{discharge_validation.get('passed', 0):,}条，反向日期为{discharge_validation.get('failed', 0):,}条；缺失或非ISO EGRESO为{discharge_validation.get('missing_or_nonISO_EGRESO', 0):,}条。按入院日期筛选，未按出院年份筛选。

原2018–2022计数与补入来源2023后的计数如下：

| 结局 | 2018–2022基线 | 补入来源2023后 |
|---|---:|---:|
| N20–N23 | {base['broad_N20_N23']:,} | {totals['broad_N20_N23']:,} |
| 严格N20 | {base['strict_N20']:,} | {totals['strict_N20']:,} |
| 记录TIPO=Q | {base['TIPO_Q']:,} | {totals['TIPO_Q']:,} |
| 无记录TIPO=Q | {base['no_TIPO_Q']:,} | {totals['no_TIPO_Q']:,} |
| 肾造瘘代码组 | {base['nephrostomy']:,} | {totals['nephrostomy']:,} |
| 输尿管导管代码组 | {base['ureteral_catheterization']:,} | {totals['ureteral_catheterization']:,} |
| 碎石代码组 | {base['urinary_stone_fragmentation_lithotripsy']:,} | {totals['urinary_stone_fragmentation_lithotripsy']:,} |

三个操作结局按预先固定的PROMED代码白名单识别，任意NUMPROMED序号均计入且不要求TIPO为Q；分组可重叠。TIPO=Q是独立的手术分类记录标志，“无Q”表示未观察到链接的TIPO=Q行。年度操作表没有可验证的实际操作日期，因此暴露按入院日期赋值，不解释为操作当日效应。

## 预定七项模型

两项州日模型估计宽口径和严格N20的ERA5-Land温度lag 0–7，加入墨西哥38日公共日期指标；固定效应为州×年月×星期几，按州及日历日期双向聚类，使用min(G州,G日期)−1自由度t近似。三项操作模型为严格N20肾造瘘、输尿管导管和碎石代码组，分别放入温度及火源PM2.5 lag 0–7与同一日期指标，维持相同固定效应、聚类和t近似；六个对比一起做BH校正。另两项同日模型沿用原州日面板和州聚类，对TIPO=Q及无Q严格N20计数估计同日温度与火源，区间沿用1.96正态近似。

| 模型 | 对比 | RR | 95% CI | p | BH q（六个操作对比） | 全期事件 | 拟合行 |
|---|---|---:|---:|---:|---:|---:|---:|
{model_table}

## 解释范围

这项敏感性检验截至2023年出院来源的边界对既有入院期结果的影响。可用来源结束后才出院、未进入SSA出院记录或其他来源版本差异仍可能造成截尾；本分析不声称覆盖全部后续出院。现有市镇流量使用其原来源窗口，未在本步骤重算。去标识逐州逐日计数、七个模型的完整系数协方差、源文件校验和与审计数字均另存。
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    log("BOUNDARY_START", workspace=str(BASE))
    require_paths()
    archive_provenance = verify_archive_inputs()
    source_provenance = sha_for_existing_source_audit()
    baseline_signatures, baseline_strict_signatures, baseline_identity_diagnoses, base_cache_audit = load_baseline_signatures()
    log("INPUT_INTEGRITY_VERIFIED",
        egresos_bytes=archive_provenance["extracted_sources"][0]["bytes"],
        procedures_bytes=archive_provenance["extracted_sources"][1]["bytes"],
        baseline_signatures=len(baseline_signatures))

    candidates, egresos_audit = scan_egresos()
    accepted, duplicate_audit, discharge_audit = reconcile_baseline(
        candidates, baseline_signatures, baseline_strict_signatures,
        baseline_identity_diagnoses
    )
    egresos_audit["accepted_admissions_by_year"] = (
        discharge_audit["accepted_admissions_by_admission_year"]
    )
    egresos_audit["accepted_admissions_by_discharge_year"] = (
        discharge_audit["accepted_admissions_by_discharge_year"]
    )
    egresos_audit["INGRE_le_EGRESO_validation"] = (
        discharge_audit["INGRE_le_EGRESO_validation"]
    )
    egresos_audit["accepted_target_admissions"] = len(accepted)
    egresos_audit["accepted_strict_N20_admissions"] = sum(
        record["strict_N20"] for record in accepted.values()
    )
    log("EGRESOS_SCAN_COMPLETE",
        source_rows=egresos_audit["records"],
        unique_annual_ids=egresos_audit["unique_annual_IDs"],
        target_signatures=egresos_audit["unique_target_admission_signatures"],
        accepted=len(accepted),
        exact_duplicates=duplicate_audit["baseline_signature_match_duplicates_excluded"],
        conflicts=duplicate_audit["same_ID_date_state_but_diagnosis_conflicts"])

    target_ids = {record["ID"] for record in accepted.values()}
    if len(target_ids) != len(accepted):
        raise ValueError("Accepted boundary admissions do not have unique annual IDs")
    procedure_flags, procedure_audit = scan_procedures(target_ids)
    log("PROCEDURES_SCAN_COMPLETE",
        rows=procedure_audit["procedure_rows_parsed"],
        linked_rows=procedure_audit["procedure_rows_linked_to_boundary_cohort"],
        q_admissions=procedure_audit["boundary_admissions_with_any_TIPO_Q_row"])

    panel, era5, _ = country_inputs.load_mexico()
    panel["state"] = panel.state.astype(str).str.zfill(2)
    panel["date"] = pd.to_datetime(panel.date).dt.normalize()
    era5["state"] = era5.state.astype(str).str.zfill(2)
    era5["date"] = pd.to_datetime(era5.date).dt.normalize()
    base_panel = panel.merge(era5, on=["state", "date"], validate="one_to_one")
    db = duckdb.connect()
    try:
        fire = db.execute(
            "SELECT state, date, fire_pm25 FROM read_parquet(?) "
            "WHERE date >= DATE '2018-01-01' AND date < DATE '2023-01-01'",
            [str(ORIGINAL_PANEL)],
        ).df()
    finally:
        db.close()
    fire["state"] = fire.state.astype(str).str.zfill(2)
    fire["date"] = pd.to_datetime(fire.date).dt.normalize()
    base_panel = base_panel.merge(fire, on=["state", "date"], validate="one_to_one")
    if len(base_panel) != EXPECTED_STATE_DAYS or base_panel.fire_pm25.isna().any():
        raise ValueError("Existing 2018-2022 state-day climate panel is incomplete")
    daily_base = pd.read_csv(
        OLD_DAILY, dtype={"state": str}, keep_default_na=False,
    )
    if not set(STRICT_DAILY_COLUMNS).issubset(daily_base.columns):
        raise ValueError("Corrected strict procedure daily counts lack required endpoints")
    daily_base["date"] = pd.to_datetime(daily_base.date).dt.normalize()
    for column in STRICT_DAILY_COLUMNS:
        daily_base[column] = pd.to_numeric(daily_base[column], errors="raise").astype("int64")
    base_panel_for_merge = base_panel[[
        "state", "date", "broad_N20_N23", "strict_N20", "stratum",
    ]].copy()
    boundary_daily, daily_audit = make_boundary_daily(
        accepted, procedure_flags, daily_base, base_panel_for_merge
    )
    log("DAILY_COUNTS_COMPLETE",
        base_strict=daily_audit["base_counts_2018_2022_source_years"]["strict_N20"],
        added_strict=daily_audit["accepted_source_2023_additions"]["strict_N20"],
        final_strict=daily_audit["augmented_counts_source_years_2018_2023"]["strict_N20"],
        added_broad=daily_audit["accepted_source_2023_additions"]["broad_N20_N23"])

    holiday = pd.read_csv(HOLIDAY_DATES, dtype={"country": str})
    mexico_days = set(pd.to_datetime(
        holiday.loc[holiday.country.eq("MX"), "date"]
    ).dt.normalize())
    if len(mexico_days) != 38:
        raise ValueError("Expected the predeclared 38 Mexico public-date controls")
    boundary_daily["date"] = pd.to_datetime(boundary_daily.date).dt.normalize()
    boundary_daily["state"] = boundary_daily.state.astype(str).str.zfill(2)
    modeled = base_panel.merge(
        boundary_daily[[
            "state", "date", "broad_N20_N23", "admissions", *GROUPS, "TIPO_Q",
            "no_TIPO_Q",
        ]],
        on=["state", "date"], how="left", validate="one_to_one",
        suffixes=("_base", ""),
    )
    if modeled[["broad_N20_N23", "admissions", *GROUPS, "TIPO_Q", "no_TIPO_Q"]].isna().any().any():
        raise ValueError("Boundary daily counts did not cover every model state-day")
    modeled["broad_N20_N23"] = modeled.broad_N20_N23.astype("int64")
    modeled["strict_N20"] = modeled.admissions.astype("int64")
    modeled["TIPO_Q"] = modeled.TIPO_Q.astype("int64")
    modeled["no_TIPO_Q"] = modeled.no_TIPO_Q.astype("int64")
    modeled["holiday"] = modeled.date.isin(mexico_days).astype("int64")
    modeled["date_id"] = modeled.date.dt.strftime("%Y%m%d")
    modeled["fire10"] = modeled.fire_pm25 / 10.0
    if not modeled.stratum.equals(
        modeled.state + "_" + modeled.date.dt.strftime("%Y%m") + "_"
        + modeled.date.dt.dayofweek.astype(str)
    ):
        raise ValueError("Boundary model strata differ from state x year-month x weekday")
    lagged = add_lags(modeled, ["temperature", "fire10"])
    if lagged[["temperature", "fire10", "holiday"]].isna().any().any():
        raise ValueError("Boundary model exposure or holiday indicator has missing values")
    calendar = lagged.sort_values(["state", "date"]).groupby("state", observed=True).date.diff().dropna()
    if not calendar.dt.days.eq(1).all():
        raise ValueError("Boundary model calendar is not consecutive by state")

    models = {}
    results = []
    fit_two_way_lag_model(
        lagged, "broad_N20_N23", "broad_N20_N23_lag07_holiday",
        daily_audit["augmented_counts_source_years_2018_2023"]["broad_N20_N23"],
        models, results,
    )
    fit_two_way_lag_model(
        lagged, "strict_N20", "strict_N20_lag07_holiday",
        daily_audit["augmented_counts_source_years_2018_2023"]["strict_N20"],
        models, results,
    )
    fit_procedure_lag_models(lagged, models, results)
    fit_same_day_models(boundary_daily, models, results)

    audit_source = {
        "current_full_acquisition_evidence": str(ARCHIVE_EVIDENCE),
        "current_full_acquisition_evidence_sha256": sha256_file(ARCHIVE_EVIDENCE),
        "historical_2018_2022": source_provenance,
        "analysis_inputs": {
            "corrected_strict_daily_counts": {
                "path": str(OLD_DAILY),
                "sha256": sha256_file(OLD_DAILY),
            },
            "original_broad_state_day_panel": {
                "path": str(country_inputs.MEXICO_PANEL_PATH),
                "sha256": sha256_file(country_inputs.MEXICO_PANEL_PATH),
            },
            "original_same_day_model_panel": {
                "path": str(ORIGINAL_PANEL),
                "sha256": sha256_file(ORIGINAL_PANEL),
            },
            "era5_land_state_day_source": {
                "path": str(country_inputs.MEXICO_ERA5_PATH),
                "sha256": sha256_file(country_inputs.MEXICO_ERA5_PATH),
            },
            "holiday_calendar_dates": {
                "path": str(HOLIDAY_DATES),
                "sha256": sha256_file(HOLIDAY_DATES),
                "Mexico_dates": len(mexico_days),
            },
        },
        "procedure_code_groups": {
            group: {
                "label": feasibility.CODE_GROUPS[group]["label"],
                "codes": sorted(GROUP_CODES[group]),
            }
            for group in GROUPS
        },
    }
    write_outputs(
        audit_source, archive_provenance, base_cache_audit, egresos_audit,
        duplicate_audit, procedure_audit, daily_audit, boundary_daily,
        results, models,
    )
    log("BOUNDARY_COMPLETE",
        models=len(models),
        contrast_rows=len(results),
        final_broad=daily_audit["augmented_counts_source_years_2018_2023"]["broad_N20_N23"],
        final_strict=daily_audit["augmented_counts_source_years_2018_2023"]["strict_N20"],
        output=str(OUT))


if __name__ == "__main__":
    main()
