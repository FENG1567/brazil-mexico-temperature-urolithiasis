#!/usr/bin/env python3
"""Rebuild aggregate-result figures and tables without fitting statistical models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath
from typing import Any


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "manifest.csv"
REQUIRED_MANIFEST_COLUMNS = {"path", "sha256", "bytes", "role"}
GENERATED_ROLE_MARKERS = {"artifact", "build", "derived", "generated", "output", "render"}
GENERATED_PATH_COMPONENTS = {".git", "__pycache__", "artifacts", "build", "dist", "generated", "outputs"}


class VerificationError(RuntimeError):
    """Raised when the payload manifest or one of its files is invalid."""


def _relative_parts(raw_path: str) -> tuple[str, ...]:
    normalized = raw_path.strip().replace("\\", "/")
    candidate = PurePosixPath(normalized)
    if not normalized or candidate.is_absolute() or (candidate.parts and ":" in candidate.parts[0]):
        raise VerificationError(f"Manifest path must be repository-relative: {raw_path!r}")
    parts = tuple(part for part in candidate.parts if part not in {"", "."})
    if not parts or ".." in parts:
        raise VerificationError(f"Manifest path escapes the repository: {raw_path!r}")
    return parts


def _is_generated(path_parts: tuple[str, ...], role: str) -> bool:
    lowered = tuple(part.casefold() for part in path_parts)
    if PurePosixPath(*lowered).as_posix() == "manifest.csv":
        return True
    if GENERATED_PATH_COMPONENTS.intersection(lowered):
        return True
    role_words = set(re.split(r"[^a-z0-9]+", role.casefold()))
    return bool(GENERATED_ROLE_MARKERS.intersection(role_words))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_payloads() -> dict[str, Any]:
    """Check manifest-declared payload bytes and hashes, skipping generated files."""
    if not MANIFEST.is_file():
        raise VerificationError(f"Payload manifest not found: {MANIFEST}")

    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or ())
        missing_columns = REQUIRED_MANIFEST_COLUMNS - fields
        if missing_columns:
            raise VerificationError(f"manifest.csv is missing columns: {sorted(missing_columns)}")

        payloads: list[tuple[Path, str, int]] = []
        seen: set[str] = set()
        skipped_generated = 0
        for line_number, row in enumerate(reader, start=2):
            raw_path = row.get("path", "") or ""
            parts = _relative_parts(raw_path)
            if _is_generated(parts, row.get("role", "") or ""):
                skipped_generated += 1
                continue

            relative = PurePosixPath(*parts).as_posix()
            if relative in seen:
                raise VerificationError(f"Duplicate payload path on manifest line {line_number}: {relative}")
            seen.add(relative)

            expected_hash = (row.get("sha256", "") or "").strip().lower()
            if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
                raise VerificationError(f"Invalid SHA-256 on manifest line {line_number}: {raw_path}")
            try:
                expected_bytes = int((row.get("bytes", "") or "").strip())
            except ValueError as exc:
                raise VerificationError(f"Invalid byte count on manifest line {line_number}: {raw_path}") from exc
            if expected_bytes < 0:
                raise VerificationError(f"Negative byte count on manifest line {line_number}: {raw_path}")

            candidate = ROOT.joinpath(*parts).resolve()
            try:
                candidate.relative_to(ROOT.resolve())
            except ValueError as exc:
                raise VerificationError(f"Manifest path resolves outside the repository: {raw_path!r}") from exc
            if not candidate.is_file():
                raise VerificationError(f"Payload file is missing: {relative}")
            payloads.append((candidate, expected_hash, expected_bytes))

    if not payloads:
        raise VerificationError("manifest.csv contains no payload files to verify")

    failures: list[str] = []
    for path, expected_hash, expected_bytes in payloads:
        relative = path.relative_to(ROOT).as_posix()
        observed_bytes = path.stat().st_size
        if observed_bytes != expected_bytes:
            failures.append(f"{relative}: expected {expected_bytes} bytes, found {observed_bytes}")
        observed_hash = _sha256(path)
        if observed_hash != expected_hash:
            failures.append(f"{relative}: SHA-256 mismatch")
    if failures:
        raise VerificationError("Payload verification failed:\n- " + "\n- ".join(failures))

    return {"verified_payload_files": len(payloads), "skipped_generated_rows": skipped_generated}


def run_builders() -> None:
    commands = [
        [sys.executable, str(ROOT / "src" / "figures" / "build_figures.py")],
        [sys.executable, str(ROOT / "src" / "reporting" / "build_tables.py")],
    ]
    for command in commands:
        subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild figures and tables from frozen aggregate results; no models are fitted."
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Verify manifest-declared payload files and hashes, then exit without generating artifacts.",
    )
    args = parser.parse_args()

    try:
        verification = verify_payloads()
        if not args.verify_only:
            run_builders()
    except (OSError, VerificationError, subprocess.CalledProcessError) as exc:
        print(f"reproduction failed: {exc}", file=sys.stderr)
        return 1

    result = {"status": "verified" if args.verify_only else "built", **verification}
    if not args.verify_only:
        result["outputs"] = ["artifacts/figures/", "artifacts/tables/"]
        result["scope"] = "Aggregate-result reproduction only; no raw-data processing or model fitting."
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
