"""Prepare or explicitly run archived model stages in an isolated workspace.

The archived model scripts are copied byte-for-byte except for narrowly scoped
path substitutions. Their model and estimation code is not rewritten here.
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE_WORK = ROOT / "analysis" / "work"
MODEL_WORKSPACE = ROOT / "artifacts" / "model_workspace"
PROCESSED = ROOT / "data" / "processed_inputs"
AGGREGATES = ROOT / "data" / "aggregate_results"

OLD_PYLIBS = "F:/analysis/climate_biliary_20261102/work/pylibs"


def _path_expr(path: Path) -> str:
    """Return a portable Python Path expression for a resolved repository path."""
    return f"Path({str(path.resolve())!r})"


def _replace_assignment(source: str, name: str, expression: str) -> str:
    """Replace exactly one top-level assignment to a known path constant."""
    tree = ast.parse(source)
    matches: list[tuple[int, int]] = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == name for target in node.targets
        ):
            matches.append((node.lineno - 1, node.end_lineno or node.lineno))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == name:
            matches.append((node.lineno - 1, node.end_lineno or node.lineno))
    if len(matches) != 1:
        raise ValueError(f"Expected one top-level {name} assignment, found {len(matches)}")

    lines = source.splitlines(keepends=True)
    start, end = matches[0]
    newline = "\r\n" if any(line.endswith("\r\n") for line in lines) else "\n"
    lines[start:end] = [f"{name} = {expression}{newline}"]
    return "".join(lines)


def _redirect_pylibs(source: str, path: Path) -> str:
    """Redirect the archive's exact historical pylibs path to the local workspace."""
    pattern = re.compile(
        r"(sys\.path\.insert\(\s*0\s*,\s*)r?['\"]"
        + re.escape(OLD_PYLIBS)
        + r"['\"](\s*\))"
    )
    def replacement(match: re.Match[str]) -> str:
        return f"{match.group(1)}str(Path({str(path.resolve())!r})){match.group(2)}"

    result, count = pattern.subn(replacement, source)
    if count != 1:
        raise ValueError(f"Expected one historical pylibs path insertion, found {count}")
    return result


def _input_specs(stage: str) -> list[dict[str, Any]]:
    """Canonical package source keys and expected materialized paths."""
    common = [
        {
            "key": "panel_composite_humidity",
            "path": "artifacts/processed_inputs/panel_composite_humidity.parquet",
            "required": True,
            "note": "Reassembled from data/processed_inputs/panel_composite_humidity.parquet.parts.",
        },
        {
            "key": "clinical_code_day_aggregate",
            "path": "data/processed_inputs/clinical_code_day_aggregate.parquet",
            "required": True,
        },
    ]
    if stage == "service":
        return common + [
            {"key": "ipca_monthly", "path": "data/aggregate_results/sensitivity_results/ipca_monthly.csv", "required": True},
            {"key": "sp_linked_daily_outcomes", "path": "data/processed_inputs/sp_linked_daily_outcomes.csv", "required": True},
            {"key": "holiday_calendar_dates", "path": "data/aggregate_results/sensitivity_results/holiday_calendar_dates.csv", "required": True},
            {"key": "acute_repair_model_results", "path": "data/aggregate_results/sensitivity_results/acute_repair_model_results.csv", "required": True},
        ]
    if stage == "nonlinear":
        return common + [
            {"key": "mexico_composite_state_day_panel", "path": "data/processed_inputs/mexico_composite_state_day_panel.parquet", "required": True},
            {"key": "mex_state_day_2015_2023", "path": "data/processed_inputs/mex_state_day_2015_2023.csv", "required": True},
            {"key": "mexico_boundary_daily_counts", "path": "data/processed_inputs/mexico_boundary_daily_counts.csv", "required": True},
            {"key": "holiday_calendar_dates", "path": "data/aggregate_results/sensitivity_results/holiday_calendar_dates.csv", "required": True},
            {"key": "n20_composite_continuous_lag", "path": "data/aggregate_results/alternative_model_results/n20_composite_continuous_lag.csv", "required": True},
            {"key": "n20_composite_groups_corrected", "path": "data/aggregate_results/alternative_model_results/n20_composite_groups_corrected.csv", "required": True},
            {
                "key": "mexico_n20_proced",
                "path": "data/processed_inputs/mexico_n20_proced.parquet",
                "required": False,
                "note": "The imported helper defines this for its Mexico strict-phenotype path; run_m3_nonlinear.py currently calls only its Brazil loader.",
            },
            {
                "key": "mexico_state_day_temperature",
                "path": "data/processed_inputs/mexico_state_day_temperature.parquet",
                "required": False,
                "note": "Defined by the imported helper but not read by run_m3_nonlinear.py's current execution path.",
            },
        ]
    raise ValueError(f"Unknown model stage: {stage}")


def _script_sources(stage: str) -> list[tuple[Path, Path]]:
    """Map archived sources to their preserved in-workspace relative structure."""
    if stage == "service":
        relative = Path("revision_20261113/service/fit_service_models.py")
        return [
            (ARCHIVE_WORK / relative, Path("work") / relative),
            (
                ROOT / "analysis/specifications/analysis_specification.json",
                Path("work/revision_20261113/service/analysis_specification.json"),
            ),
        ]
    if stage == "nonlinear":
        return [
            (
                ARCHIVE_WORK / "revision_20261113/nonlinear/run_m3_nonlinear.py",
                Path("work/revision_20261113/nonlinear/run_m3_nonlinear.py"),
            ),
            (
                ROOT / "analysis/specifications/prespecified_specification.json",
                Path("work/revision_20261113/nonlinear/prespecified_specification.json"),
            ),
            (
                ARCHIVE_WORK / "upgrade_20261112/mexico_strict_phenotype.py",
                Path("work/upgrade_20261112/mexico_strict_phenotype.py"),
            ),
        ]
    raise ValueError(f"Unknown model stage: {stage}")


def _materialize_nonlinear_clinical(workspace: Path) -> dict[str, Any]:
    """Make the strict helper's historical ROOT-relative clinical path available."""
    source = ROOT / "data/processed_inputs/clinical_code_day_aggregate.parquet"
    target = workspace / "work/upgrade_20261112/clinical_code_day_aggregate.parquet"
    record: dict[str, Any] = {
        "key": "nonlinear_clinical_code_day_aggregate",
        "source": source.relative_to(ROOT).as_posix(),
        "workspace_path": target.relative_to(workspace).as_posix(),
        "required": True,
        "exists": False,
        "method": None,
    }
    if not source.is_file():
        return record

    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file() and os.path.samefile(source, target):
        record.update(exists=True, method="hardlink")
        return record

    temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
    try:
        temporary.unlink(missing_ok=True)
        try:
            os.link(source, temporary)
            method = "hardlink"
        except OSError:
            shutil.copy2(source, temporary)
            method = "copy"
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    record.update(exists=target.is_file(), method=method if target.is_file() else None)
    return record


def _rewrite_script(stage: str, relative_path: Path, source: str, workspace: Path) -> str:
    inputs = {item["key"]: ROOT / item["path"] for item in _input_specs(stage)}
    if stage == "service" and relative_path.name == "fit_service_models.py":
        source = _redirect_pylibs(source, ROOT / "vendor" / "pylibs")
        replacements = {
            "OUT": _path_expr(workspace / "outputs" / "service"),
            "DATA": _path_expr(PROCESSED),
            "EXPOSURE_PATH": _path_expr(ROOT / "artifacts/processed_inputs/panel_composite_humidity.parquet"),
            "CLINICAL_PATH": _path_expr(inputs["clinical_code_day_aggregate"]),
            "PRICE_PATH": _path_expr(inputs["ipca_monthly"]),
            "SP_LINK_PATH": _path_expr(inputs["sp_linked_daily_outcomes"]),
            "HOLIDAY_PATH": _path_expr(inputs["holiday_calendar_dates"]),
            "LEGACY_RESULTS_PATH": _path_expr(inputs["acute_repair_model_results"]),
        }
        for name, expression in replacements.items():
            source = _replace_assignment(source, name, expression)
    elif stage == "nonlinear" and relative_path.name == "run_m3_nonlinear.py":
        source = _redirect_pylibs(source, ROOT / "vendor" / "pylibs")
        replacements = {
            "OUT": _path_expr(workspace / "outputs" / "nonlinear"),
            "BRAZIL_PANEL": _path_expr(ROOT / "artifacts/processed_inputs/panel_composite_humidity.parquet"),
            "MEXICO_PANEL": _path_expr(inputs["mexico_composite_state_day_panel"]),
            "MEXICO_ERA5": _path_expr(inputs["mex_state_day_2015_2023"]),
            "MEXICO_DAILY": _path_expr(inputs["mexico_boundary_daily_counts"]),
            "HOLIDAY_CSV": _path_expr(inputs["holiday_calendar_dates"]),
            "OLD_LINEAR": _path_expr(inputs["n20_composite_continuous_lag"]),
            "OLD_COMPOUND": _path_expr(inputs["n20_composite_groups_corrected"]),
        }
        for name, expression in replacements.items():
            source = _replace_assignment(source, name, expression)
    elif stage == "nonlinear" and relative_path.name == "mexico_strict_phenotype.py":
        source = _redirect_pylibs(source, ROOT / "vendor" / "pylibs")
        replacements = {
            "OUT": _path_expr(workspace / "outputs" / "strict_helper"),
            "SOURCE_WORK": _path_expr(PROCESSED),
            "BRAZIL_PANEL_PATH": _path_expr(ROOT / "artifacts/processed_inputs/panel_composite_humidity.parquet"),
            "BRAZIL_CLINICAL_PATH": _path_expr(inputs["clinical_code_day_aggregate"]),
            "MEXICO_PANEL_PATH": _path_expr(inputs["mexico_composite_state_day_panel"]),
            "MEXICO_CLINICAL_PATH": _path_expr(inputs["mexico_n20_proced"]),
            "MEXICO_NASA_PATH": _path_expr(inputs["mexico_state_day_temperature"]),
            "MEXICO_ERA5_PATH": _path_expr(inputs["mex_state_day_2015_2023"]),
        }
        for name, expression in replacements.items():
            source = _replace_assignment(source, name, expression)
    ast.parse(source)
    return source


def _prepare(stage: str) -> tuple[Path, dict[str, Any]]:
    workspace = MODEL_WORKSPACE / stage
    workspace.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    missing_sources: list[str] = []
    for source_path, relative_path in _script_sources(stage):
        target_path = workspace / relative_path
        if not source_path.is_file():
            missing_sources.append(source_path.relative_to(ROOT).as_posix())
            continue
        target_path.parent.mkdir(parents=True, exist_ok=True)
        if source_path.suffix == ".py":
            original = source_path.read_text(encoding="utf-8")
            rewritten = _rewrite_script(stage, relative_path, original, workspace)
            target_path.write_text(rewritten, encoding="utf-8")
        else:
            shutil.copy2(source_path, target_path)
        copied.append(relative_path.as_posix())

    inputs = []
    missing_required = []
    for item in _input_specs(stage):
        path = ROOT / item["path"]
        present = path.is_file()
        record = {**item, "exists": present}
        inputs.append(record)
        if item["required"] and not present:
            missing_required.append(item["path"])

    materialized_inputs = []
    if stage == "nonlinear":
        materialized = _materialize_nonlinear_clinical(workspace)
        materialized_inputs.append(materialized)
        if materialized["required"] and not materialized["exists"]:
            missing_required.append(materialized["workspace_path"])

    expected_reassembly = ROOT / "artifacts/processed_inputs/panel_composite_humidity.parquet"
    input_inventory = PROCESSED / "inventory.json"
    reassembly_helper = ROOT / "src/data/reassemble.py"
    if not expected_reassembly.is_file():
        reassembly_note = (
            "Reassembled panel is not present. If the parts and chunk index are available, "
            "run `python src/data/reassemble.py` before model execution."
        )
    else:
        reassembly_note = None

    manifest: dict[str, Any] = {
        "source_release": "SRR-20261114-C3",
        "stage": stage,
        "mode": "prepare-only",
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "workspace": workspace.relative_to(ROOT).as_posix(),
        "archive_work": ARCHIVE_WORK.relative_to(ROOT).as_posix(),
        "copied_stage_files": copied,
        "missing_stage_sources": missing_sources,
        "inputs": inputs,
        "materialized_inputs": materialized_inputs,
        "required_inputs_missing": missing_required,
        "reassembly_helper_exists": reassembly_helper.is_file(),
        "input_inventory_exists": input_inventory.is_file(),
        "reassembly_note": reassembly_note,
        "scientific_code_modified": False,
    }
    manifest_path = workspace / "preparation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return workspace, manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare archived service/nonlinear models; execution requires --run."
    )
    parser.add_argument("stage", choices=("service", "nonlinear"), help="model stage to prepare")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--run", action="store_true", help="execute the copied model after complete preflight")
    mode.add_argument("--prepare-only", action="store_true", help="copy, path-map, and preflight without fitting (default)")
    return parser


def main() -> int:
    args = _parser().parse_args()
    workspace, manifest = _prepare(args.stage)
    manifest_path = workspace / "preparation_manifest.json"
    print(f"Prepared {args.stage} workspace: {workspace}")
    print(f"Preparation manifest: {manifest_path}")

    if manifest["missing_stage_sources"]:
        print("Missing archived stage sources: " + ", ".join(manifest["missing_stage_sources"]), file=sys.stderr)
        return 2
    if not args.run:
        if manifest["required_inputs_missing"]:
            print("Required inputs not yet present: " + ", ".join(manifest["required_inputs_missing"]))
        if manifest["reassembly_note"]:
            print(manifest["reassembly_note"])
        return 0
    if manifest["required_inputs_missing"]:
        print("Cannot run; required inputs are missing: " + ", ".join(manifest["required_inputs_missing"]), file=sys.stderr)
        if manifest["reassembly_note"]:
            print(manifest["reassembly_note"], file=sys.stderr)
        return 2

    if args.stage == "service":
        script = workspace / "work/revision_20261113/service/fit_service_models.py"
    else:
        script = workspace / "work/revision_20261113/nonlinear/run_m3_nonlinear.py"
    print(f"Running copied model: {script}", flush=True)
    return subprocess.run([sys.executable, str(script)], cwd=workspace, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
