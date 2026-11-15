from __future__ import annotations

import importlib.util
import json
import os
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORK = ROOT / "work" / "revision_20261113" / "descriptive"
OUT = ROOT / "outputs" / "revision_20261113" / "descriptive"
STATUS = WORK / "mexico_demographic_reclassification_status.json"
AGGREGATE_STATUS = WORK / "demographic_source_status.json"
AUDIT = OUT / "demographic_source_audit.json"
AGGREGATOR = WORK / "aggregate_case_demographics.py"


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: object, *, atomic: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    target = path.with_suffix(path.suffix + ".tmp") if atomic else path
    target.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    if atomic:
        target.replace(path)


def main() -> None:
    started = stamp()
    write_json(STATUS, {
        "status": "running", "pid": os.getpid(), "started_at_utc": started,
        "outputs_are_aggregate_only": True, "person_ids_read_or_written": False,
    })
    try:
        spec = importlib.util.spec_from_file_location("aggregate_case_demographics", AGGREGATOR)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not load demographic aggregator at {AGGREGATOR}")
        aggregator = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = aggregator
        spec.loader.exec_module(aggregator)

        frame, mexico_meta = aggregator.run_mexico()
        if int(mexico_meta["cache_reconciliation_max_abs_difference"]) != 0:
            raise ValueError("Mexico Egresos target counts do not reconcile to the 2018–2022 valid-state cache")
        if int(frame.admissions.sum()) != int(mexico_meta["baseline_admissions"]):
            raise ValueError("Mexico aggregate row total does not match its scan metadata")
        if not AUDIT.exists():
            raise FileNotFoundError(f"Expected existing aggregate-only audit to preserve Brazil metadata: {AUDIT}")

        audit = json.loads(AUDIT.read_text(encoding="utf-8"))
        if not isinstance(audit.get("Brazil"), dict):
            raise ValueError("Existing demographic audit has no Brazil aggregate metadata to preserve")
        mexico_meta["age_codebook_sources"] = {
            "current_variable_record": "https://riisp.insp.mx/nada/index.php/catalog/18/variable/V11",
            "historical_2014_validation_rules_v2_2": "https://riisp.insp.mx/nada/files_nesstar/SAEH_SSA/Referencias/NOTAS_TECNICAS/GII-B002-SAEH-02-02.pdf",
        }
        mexico_meta["sex_codebook_source"] = "https://riisp.insp.mx/nada/index.php/catalog/18/variable/V14"
        mexico_meta["reclassification_pid"] = os.getpid()
        mexico_meta["reclassification_completed_at_utc"] = stamp()
        mexico_meta["cveedad_unmapped_counts"] = {
            str(code): int(count)
            for code, count in frame.groupby("cveedad_code", dropna=False).admissions.sum().items()
            if str(code) in {"4", "5"} or str(code) not in {"0", "1", "2", "3", "9", "missing_or_blank"}
        }
        audit["Mexico"] = mexico_meta
        audit["mexico_admissions_total"] = int(frame.admissions.sum())
        audit["status"] = "complete"
        audit["completed_at_utc"] = stamp()
        audit["aggregate_only"] = True
        audit["person_ids_read_or_written"] = False
        write_json(AUDIT, audit, atomic=True)

        if AGGREGATE_STATUS.exists():
            aggregate_status = json.loads(AGGREGATE_STATUS.read_text(encoding="utf-8"))
            aggregate_status.update({
                "status": "complete", "started_at_utc": started, "completed_at_utc": stamp(),
                "pid": os.getpid(), "Mexico": mexico_meta,
                "mexico_admissions_total": int(frame.admissions.sum()),
            })
            aggregate_status.pop("elapsed_seconds", None)
            write_json(AGGREGATE_STATUS, aggregate_status, atomic=True)

        result = {
            "status": "complete", "pid": os.getpid(), "started_at_utc": started,
            "completed_at_utc": stamp(), "source_rows_scanned": int(mexico_meta["source_rows_scanned"]),
            "bytes_scanned": int(mexico_meta["bytes"]),
            "valid_state_target_admissions": int(mexico_meta["valid_state_target_admissions"]),
            "cache_reconciliation_max_abs_difference": int(mexico_meta["cache_reconciliation_max_abs_difference"]),
            "CVEEDAD_codes": mexico_meta["cveedad_codes"],
            "age_categories": mexico_meta["age_categories"],
            "unmapped_CVEEDAD_counts": mexico_meta["cveedad_unmapped_counts"],
            "person_ids_read_or_written": False,
        }
        write_json(STATUS, result, atomic=True)
        print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    except Exception as exc:
        failure = {
            "status": "failed", "pid": os.getpid(), "started_at_utc": started,
            "failed_at_utc": stamp(), "error": repr(exc),
            "traceback": traceback.format_exc(), "person_ids_read_or_written": False,
        }
        write_json(STATUS, failure, atomic=True)
        raise


if __name__ == "__main__":
    main()
