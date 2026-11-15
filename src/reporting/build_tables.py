#!/usr/bin/env python3
"""Export the released main-table CSVs and editable DOCX tables.

For local final DOCX QA, use WPS Writer first: open each final DOCX, export or
render it, and inspect the rendered pages. Record the WPS limitation and the
fallback runtime if another application is required.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "data" / "aggregate_results" / "main_tables"
OUTPUT_DIR = ROOT / "artifacts" / "tables"
METADATA_PATH = ROOT / "data" / "table_metadata.json"
TABLE_NAMES = ("table1", "table2", "table3")


def read_and_validate_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a UTF-8 CSV and check that every record matches its own header."""
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream)
            header = next(reader, None)
            if header is None:
                raise ValueError(f"Source table is empty: {path}")
            if not header or any(not cell.strip() for cell in header):
                raise ValueError(f"Table headers must be non-empty: {path}")
            if len(set(header)) != len(header):
                raise ValueError(f"Table headers must be unique: {path}")
            rows: list[list[str]] = []
            for row in reader:
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError(
                        f"Malformed CSV row in {path} at line {reader.line_num}: "
                        f"expected {len(header)} fields, found {len(row)}"
                    )
                rows.append(row)
    except UnicodeDecodeError as exc:
        raise ValueError(f"Source table must be UTF-8 encoded: {path}") from exc
    return header, rows


def write_docx(title: str, headers: list[str], rows: list[list[str]], destination: Path, note: str = "") -> None:
    """Write one editable Word table, importing python-docx only when needed."""
    try:
        from docx import Document
        from docx.shared import Inches, Pt
        from docx.enum.section import WD_ORIENT
        from docx.oxml import OxmlElement
    except ImportError as exc:
        raise RuntimeError(
            "Creating editable table DOCX files requires the python-docx dependency."
        ) from exc

    document = Document()
    section = document.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = Inches(11.7), Inches(8.3)
    section.top_margin = section.bottom_margin = Inches(0.5)
    section.left_margin = section.right_margin = Inches(0.5)
    document.styles['Normal'].font.name = 'Times New Roman'
    document.styles['Normal'].font.size = Pt(9)
    document.add_heading(title, level=1)
    table = document.add_table(rows=1, cols=len(headers))
    repeating_header = OxmlElement('w:tblHeader')
    table.rows[0]._tr.get_or_add_trPr().append(repeating_header)
    try:
        table.style = "Table Grid"
    except KeyError:
        # The Word default table style still produces an editable table.
        pass
    for cell, value in zip(table.rows[0].cells, headers):
        cell.text = value
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.bold = True
    for values in rows:
        cells = table.add_row().cells
        for cell, value in zip(cells, values):
            cell.text = value
    if note:
        paragraph = document.add_paragraph(note)
        for run in paragraph.runs:
            run.font.size = Pt(8)
    document.save(destination)


def build_tables() -> dict[str, Any]:
    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8"))
    sources: dict[str, tuple[Path, list[str], list[list[str]]]] = {}
    for name in TABLE_NAMES:
        source = INPUT_DIR / f"{name}.csv"
        if not source.is_file():
            raise FileNotFoundError(f"Required aggregate table not found: {source}")
        headers, rows = read_and_validate_table(source)
        sources[name] = (source, headers, rows)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs: list[str] = []
    for name, (source, headers, rows) in sources.items():
        csv_destination = OUTPUT_DIR / f"{name}.csv"
        docx_destination = OUTPUT_DIR / f"{name}.docx"
        # Byte-for-byte copy preserves the released CSV, including encoding and line endings.
        shutil.copyfile(source, csv_destination)
        entry = metadata[name]
        write_docx(entry["title"], headers, rows, docx_destination, entry["note"])
        outputs.extend(
            [
                csv_destination.relative_to(ROOT).as_posix(),
                docx_destination.relative_to(ROOT).as_posix(),
            ]
        )
    return {"tables": list(TABLE_NAMES), "outputs": outputs}


def main() -> int:
    parser = argparse.ArgumentParser(description="Export the three released aggregate-result tables.")
    parser.parse_args()
    try:
        result = build_tables()
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"table export failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
