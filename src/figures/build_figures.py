#!/usr/bin/env python3
"""Build the five revision figures from released, auditable result rows."""

from __future__ import annotations

import csv
import argparse
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib as mpl
from matplotlib import font_manager

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


_SCRIPT_PARENTS = Path(__file__).resolve().parents
ROOT = _SCRIPT_PARENTS[2] if len(_SCRIPT_PARENTS) > 2 else Path.cwd()
TABLES = ROOT / "data" / "figure_inputs" / "frozen" / "tables"
UPGRADE = ROOT / "data" / "figure_inputs" / "upgrade"
REVISION = ROOT / "data" / "figure_inputs" / "revision"
OUTDIR = ROOT / "artifacts" / "figures"
SOURCES_OUT = ROOT / "artifacts" / "provenance" / "figure_sources.json"
SVG_DIR = OUTDIR / "svg"

INK = "#263746"
BLUE = "#345C78"
TEAL = "#2E7D78"
RUST = "#B45E49"
SLATE = "#687C88"
PALE = "#F4F7F8"
RULE = "#CCD6DA"
GRID = "#E5EBEE"
WHITE = "#FFFFFF"
GOLD = "#A77C2E"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans", "sans-serif"],
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "font.size": 9.3,
        "axes.labelsize": 9.3,
        "axes.titlesize": 10.0,
        "xtick.labelsize": 9.0,
        "ytick.labelsize": 9.0,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": WHITE,
        "axes.facecolor": WHITE,
        "savefig.facecolor": WHITE,
        "savefig.dpi": 300,
    }
)


def register_optional_fonts() -> None:
    """Load locally staged fonts when a caller supplies FIGURE_FONT_DIR."""
    font_dir = os.environ.get("FIGURE_FONT_DIR")
    if not font_dir:
        return
    directory = Path(font_dir).expanduser()
    if not directory.is_dir():
        raise FileNotFoundError(f"FIGURE_FONT_DIR does not exist: {directory}")
    for font_path in sorted([*directory.glob("*.ttf"), *directory.glob("*.otf")]):
        font_manager.fontManager.addfont(str(font_path))


register_optional_fonts()

TABLE_CACHE: dict[str, list[dict[str, str]]] = {}
JSON_CACHE: dict[str, Any] = {}
FILE_META: dict[str, dict[str, Any]] = {}
SOURCE_RECORDS: list[dict[str, Any]] = []


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def configure_paths(repo_root: Path, outdir: Path) -> None:
    """Bind inputs to a portable package root and exports to a caller path."""
    global ROOT, TABLES, UPGRADE, REVISION, OUTDIR, SOURCES_OUT, SVG_DIR
    ROOT = repo_root.expanduser().resolve()
    TABLES = ROOT / "data" / "figure_inputs" / "frozen" / "tables"
    UPGRADE = ROOT / "data" / "figure_inputs" / "upgrade"
    REVISION = ROOT / "data" / "figure_inputs" / "revision"
    OUTDIR = outdir.expanduser().resolve()
    SOURCES_OUT = (ROOT / "artifacts/provenance/figure_sources.json" if OUTDIR == ROOT / "artifacts/figures" else OUTDIR / "figure_sources.json")
    SVG_DIR = OUTDIR / "svg"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_meta(path: Path) -> dict[str, Any]:
    key = rel(path)
    if key not in FILE_META:
        FILE_META[key] = {
            "path": key,
            "sha256": sha256(path),
            "modified_at_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"),
            "bytes": path.stat().st_size,
        }
    return FILE_META[key]


def read_csv(path: Path) -> list[dict[str, str]]:
    key = rel(path)
    if key not in TABLE_CACHE:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            TABLE_CACHE[key] = list(csv.DictReader(stream))
    return TABLE_CACHE[key]


def read_json(path: Path) -> Any:
    key = rel(path)
    if key not in JSON_CACHE:
        JSON_CACHE[key] = json.loads(path.read_text(encoding="utf-8-sig"))
    return JSON_CACHE[key]


def register_row(
    figure: str,
    panel: str,
    label: str,
    path: Path,
    selector: dict[str, str],
    columns: list[str],
) -> dict[str, str]:
    matches = [
        (line_no, row)
        for line_no, row in enumerate(read_csv(path), start=2)
        if all(row.get(key) == value for key, value in selector.items())
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one source row in {path.name} for {selector}; found {len(matches)}")
    line_no, row = matches[0]
    missing = [column for column in columns if column not in row]
    if missing:
        raise ValueError(f"Missing source columns in {path.name}: {missing}")
    SOURCE_RECORDS.append(
        {
            "figure": figure,
            "panel": panel,
            "plotted_item": label,
            "source_file": source_meta(path),
            "csv_line": line_no,
            "row_selector": selector,
            "source_columns": columns,
            "source_values": {column: row[column] for column in columns},
        }
    )
    return row


def register_series(
    figure: str,
    panel: str,
    label: str,
    path: Path,
    selector: dict[str, str],
    columns: list[str],
    *,
    order_by: str | None = None,
) -> list[dict[str, str]]:
    matches = [
        (line_no, row)
        for line_no, row in enumerate(read_csv(path), start=2)
        if all(row.get(key) == value for key, value in selector.items())
    ]
    if not matches:
        raise ValueError(f"No source rows in {path.name} for {selector}")
    if order_by:
        matches.sort(key=lambda item: float(item[1][order_by]))
    missing = [column for column in columns if any(column not in row for _, row in matches)]
    if missing:
        raise ValueError(f"Missing series columns in {path.name}: {missing}")
    statuses = sorted({row.get("status", "") for _, row in matches if row.get("status")})
    SOURCE_RECORDS.append(
        {
            "figure": figure,
            "panel": panel,
            "plotted_item": label,
            "source_file": source_meta(path),
            "row_selector": selector,
            "csv_lines": {"first": matches[0][0], "last": matches[-1][0]},
            "source_rows": len(matches),
            "series_order": order_by,
            "source_columns": columns,
            "observed_status_values": statuses,
        }
    )
    return [row for _, row in matches]


def nested_value(value: Any, path_parts: Iterable[Any]) -> Any:
    current = value
    for part in path_parts:
        current = current[part]
    return current


def register_json_value(
    figure: str,
    panel: str,
    label: str,
    path: Path,
    pointer: list[Any],
) -> Any:
    value = nested_value(read_json(path), pointer)
    SOURCE_RECORDS.append(
        {
            "figure": figure,
            "panel": panel,
            "plotted_item": label,
            "source_file": source_meta(path),
            "json_pointer_parts": pointer,
            "source_value": value,
        }
    )
    return value


def register_derived_value(
    figure: str,
    panel: str,
    label: str,
    value: int | float,
    expression: str,
    inputs: list[dict[str, Any]],
) -> int | float:
    SOURCE_RECORDS.append(
        {
            "figure": figure,
            "panel": panel,
            "plotted_item": label,
            "derived_value": value,
            "expression": expression,
            "input_values": inputs,
        }
    )
    return value


def markdown_table_value(
    figure: str,
    panel: str,
    label: str,
    path: Path,
    header: str,
    key_column: str,
    key_value: str,
    value_column: str,
) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == header)
    except StopIteration as exc:
        raise ValueError(f"Markdown table header not found: {header}") from exc
    headers = [cell.strip() for cell in header.strip().strip("|").split("|")]
    matches: list[tuple[int, dict[str, str]]] = []
    for index in range(start + 2, len(lines)):
        raw = lines[index].strip()
        if not raw.startswith("|"):
            if matches:
                break
            continue
        cells = [cell.strip() for cell in raw.strip("|").split("|")]
        if len(cells) == len(headers) and not all(re.fullmatch(r":?-{3,}:?", cell.replace(" ", "")) for cell in cells):
            matches.append((index + 1, dict(zip(headers, cells))))
    hits = [(line, row) for line, row in matches if row.get(key_column) == key_value]
    if len(hits) != 1:
        raise ValueError(f"Expected one markdown table row {key_value!r} in {path.name}; found {len(hits)}")
    line, row = hits[0]
    if value_column not in row:
        raise ValueError(f"Markdown field {value_column!r} missing from {path.name}")
    SOURCE_RECORDS.append(
        {
            "figure": figure,
            "panel": panel,
            "plotted_item": label,
            "source_file": source_meta(path),
            "source_line": line,
            "source_markdown_table": header,
            "row_selector": {key_column: key_value},
            "source_columns": headers,
            "source_values": row,
            "selected_value_column": value_column,
        }
    )
    return row[value_column]


def num(row: dict[str, str], field: str) -> float:
    return float(row[field])


def fmt_n(raw: Any) -> str:
    return f"{int(float(str(raw).replace(',', ''))):,}"


def fmt_est(row: dict[str, str], estimate: str = "RR", low: str = "CI_low", high: str = "CI_high", digits: int = 3) -> str:
    return f"{num(row, estimate):.{digits}f} ({num(row, low):.{digits}f}–{num(row, high):.{digits}f})"


def add_title(fig: plt.Figure, title: str, subtitle: str) -> None:
    height_pt = fig.get_figheight() * 72
    fig.text(0.05, 1 - 8 / height_pt, title, ha="left", va="top", fontsize=13.0, fontweight="bold", color=INK)
    fig.text(0.05, 1 - 28 / height_pt, subtitle, ha="left", va="top", fontsize=9.3, fontweight="bold", color=INK)


def heading(fig: plt.Figure, letter: str, title: str, x: float, y: float, size: float = 9.5) -> None:
    label = fig.text(x, y, letter, ha="left", va="center", fontsize=18.0, fontweight="bold", color=INK)
    label.set_gid("panel-letter")
    fig.text(x + 0.075, y, title, ha="left", va="center", fontsize=size, fontweight="bold", color=INK)


def box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    detail: str,
    color: str,
    *,
    title_size: float = 10.0,
    detail_size: float = 9.2,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.007,rounding_size=0.008",
        linewidth=0.75, edgecolor=color, facecolor=PALE,
    )
    ax.add_patch(patch)
    ax.text(x + 0.018, y + height * 0.68, title, ha="left", va="center", fontsize=title_size, fontweight="bold", color=color, linespacing=1.30, multialignment="left")
    ax.text(x + 0.018, y + height * 0.29, detail, ha="left", va="center", fontsize=detail_size, fontweight="bold", color=INK, linespacing=1.15, multialignment="left")


def draw_forest(
    fig: plt.Figure,
    rect: tuple[float, float, float, float],
    rows: list[dict[str, Any]],
    *,
    xlim: tuple[float, float],
    ticks: list[float],
    label_x: float,
    value_x: float,
    label_align: str = "right",
    label_size: float = 9.2,
    value_size: float = 9.2,
    value_digits: int = 3,
    scale: str = "linear",
    markers: list[str] | None = None,
) -> plt.Axes:
    ax = fig.add_axes(rect)
    if scale == "log":
        # Every estimate and interval endpoint must be strictly positive before using a logarithmic axis.
        if any(num(row, row.get("estimate_field", "RR")) <= 0 or num(row, row.get("low_field", "CI_low")) <= 0 or num(row, row.get("high_field", "CI_high")) <= 0 for row in rows):
            raise ValueError("Log-axis RRs and CI limits must be strictly positive")
        ax.set_xscale("log")
        ax.minorticks_off()
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda value, _position: f"{value:g}"))
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.set_yticks([])
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", pad=3.0, length=3.0, colors=INK, labelsize=9.0)
    for tick_label in ax.get_xticklabels():
        tick_label.set_fontweight("bold")
    ax.grid(axis="x", color=GRID, linewidth=0.45, zorder=0)
    ax.axvline(1.0, color=INK, linewidth=0.7, linestyle=(0, (2.0, 2.0)), zorder=1)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for i, row in enumerate(rows):
        y = len(rows) - 1 - i
        est_field = row.get("estimate_field", "RR")
        low_field = row.get("low_field", "CI_low")
        high_field = row.get("high_field", "CI_high")
        est = num(row, est_field)
        low = num(row, low_field)
        high = num(row, high_field)
        if low < xlim[0] or high > xlim[1]:
            raise ValueError(
                f"CI [{low:g}, {high:g}] falls outside plotted x limits {xlim}; "
                "expand the axis so the full interval is visible"
            )
        color = row.get("plot_color", BLUE)
        marker = markers[i] if markers else row.get("plot_marker", "o")
        ax.plot([low, high], [y, y], color=color, linewidth=1.45, solid_capstyle="round", zorder=2)
        ax.plot([low, low], [y - 0.09, y + 0.09], color=color, linewidth=0.9, zorder=2)
        ax.plot([high, high], [y - 0.09, y + 0.09], color=color, linewidth=0.9, zorder=2)
        ax.plot(est, y, marker=marker, markersize=5.0, color=color, markeredgecolor=WHITE, markeredgewidth=0.45, zorder=3)
        y_fig = rect[1] + rect[3] * ((y + 0.5) / len(rows))
        fig.text(label_x, y_fig, row["plot_label"], ha=label_align, va="center", fontsize=label_size, fontweight="bold", color=INK)
        fig.text(value_x, y_fig, fmt_est(row, est_field, low_field, high_field, value_digits), ha="left", va="center", fontsize=value_size, fontweight="bold", color=INK)
    return ax


def register_fig1_data() -> dict[str, Any]:
    descriptive_flow = REVISION / "descriptive" / "sample_flow_summary.csv"
    brazil = register_row(
        "fig1", "A", "Brazil primary outcome sample frame", descriptive_flow,
        {"country": "Brazil", "stage": "primary_panel_full_window"},
        ["period", "unit", "source_or_input_n", "valid_geography_n", "exposure_complete_n", "events_broad_or_named", "note"],
    )
    service = register_row(
        "fig1", "A", "Brazil strict-N20 service sample", descriptive_flow,
        {"country": "Brazil", "stage": "service_full_window_area_day"},
        ["period", "unit", "source_or_input_n", "valid_geography_n", "exposure_complete_n", "lag0_7_complete_n", "events_broad_or_named", "events_urgent", "events_after_fit", "rows_after_fit", "note"],
    )
    mexico_frame = register_row(
        "fig1", "B", "Mexico baseline raw-to-valid and municipality exposure frame", descriptive_flow,
        {"country": "Mexico", "stage": "baseline_raw_cache"},
        ["period", "unit", "source_or_input_n", "valid_geography_n", "exposure_complete_n", "events_broad_or_named", "note"],
    )
    municipality_frame = register_row(
        "fig1", "C", "Mexico municipality exposure-complete branch", descriptive_flow,
        {"country": "Mexico", "stage": "municipality_exposure_branch"},
        ["period", "unit", "source_or_input_n", "valid_geography_n", "exposure_complete_n", "events_broad_or_named", "note"],
    )
    mexico_addition_summary = register_row(
        "fig1", "B", "Mexico 2023 source-boundary additions summary", descriptive_flow,
        {"country": "Mexico", "stage": "source_2023_boundary_additions"},
        ["period", "unit", "source_or_input_n", "valid_geography_n", "events_broad_or_named", "note"],
    )

    admission_audit = UPGRADE / "mexico_admission_date_audit.json"
    boundary_audit = UPGRADE / "mexico_boundary_audit.json"
    mexico_raw = int(mexico_frame["source_or_input_n"])
    mexico_valid = int(mexico_frame["valid_geography_n"])
    SOURCE_RECORDS.append({"figure": "fig1", "panel": "B", "plotted_item": "Mexico admission-date audit verification", "source_file": source_meta(admission_audit), "json_pointer_parts": ["checks", 0, "records"], "source_value": int(nested_value(read_json(admission_audit), ["checks", 0, "records"]))})
    SOURCE_RECORDS.append({"figure": "fig1", "panel": "B", "plotted_item": "Mexico boundary audit passed flag", "source_file": source_meta(boundary_audit), "json_pointer_parts": ["passed"], "source_value": nested_value(read_json(boundary_audit), ["passed"])})
    mexico_additions = int(register_json_value("fig1", "B", "Mexico new source admissions", boundary_audit, ["cross_source_duplicate_audit", "new_unique_admissions_retained"]))
    mexico_added_strict = int(register_json_value("fig1", "B", "Mexico new strict-N20 source admissions", boundary_audit, ["egresos_2023_audit", "accepted_strict_N20_admissions"]))
    reconciliation = read_json(boundary_audit)["daily_count_reconciliation"]["augmented_counts_source_years_2018_2023"]
    mexico_expanded = int(register_json_value("fig1", "B", "Mexico source-expanded N20–N23 admissions", boundary_audit, ["daily_count_reconciliation", "augmented_counts_source_years_2018_2023", "broad_N20_N23"]))
    mexico_strict = int(register_json_value("fig1", "B", "Mexico source-expanded strict-N20 admissions", boundary_audit, ["daily_count_reconciliation", "augmented_counts_source_years_2018_2023", "strict_N20"]))
    eligible_difference = register_derived_value(
        "fig1", "B", "Mexico raw-to-valid qualification difference", mexico_raw - mexico_valid,
        "40,907 raw cached admissions − 40,415 valid state/date admissions",
        [{"path": rel(descriptive_flow), "row_selector": {"country": "Mexico", "stage": "baseline_raw_cache"}, "field": "source_or_input_n", "value": mexico_raw}, {"path": rel(descriptive_flow), "row_selector": {"country": "Mexico", "stage": "baseline_raw_cache"}, "field": "valid_geography_n", "value": mexico_valid}],
    )
    municipality_missing = register_derived_value(
        "fig1", "C", "Valid-state admissions without municipality exposure", int(float(municipality_frame["source_or_input_n"])) - int(float(municipality_frame["exposure_complete_n"])),
        "40,415 valid-state admissions − 39,809 complete-municipality-exposure admissions",
        [{"path": rel(descriptive_flow), "row_selector": {"country": "Mexico", "stage": "municipality_exposure_branch"}, "field": "source_or_input_n", "value": municipality_frame["source_or_input_n"]}, {"path": rel(descriptive_flow), "row_selector": {"country": "Mexico", "stage": "municipality_exposure_branch"}, "field": "exposure_complete_n", "value": municipality_frame["exposure_complete_n"]}],
    )
    municipality_match = re.search(r"across\s+([\d,]+)\s+municipalities", municipality_frame["note"])
    if not municipality_match:
        raise ValueError("Could not resolve municipality count from the descriptive sample-flow note")
    municipality_n = int(municipality_match.group(1).replace(",", ""))
    SOURCE_RECORDS.append({"figure": "fig1", "panel": "C", "plotted_item": "Mexico municipalities with complete exposure", "source_file": source_meta(descriptive_flow), "row_selector": {"country": "Mexico", "stage": "municipality_exposure_branch"}, "source_column": "note", "source_text": municipality_frame["note"], "parsed_value": municipality_n})
    service_days = int(service["source_or_input_n"]) // int(service["valid_geography_n"])
    register_derived_value(
        "fig1", "A", "Calendar days in the Brazil service panel", service_days,
        "821,700 health-region-days ÷ 450 health regions",
        [{"path": rel(descriptive_flow), "row_selector": {"country": "Brazil", "stage": "service_full_window_area_day"}, "field": "source_or_input_n", "value": service["source_or_input_n"]}, {"path": rel(descriptive_flow), "row_selector": {"country": "Brazil", "stage": "service_full_window_area_day"}, "field": "valid_geography_n", "value": service["valid_geography_n"]}],
    )
    strict_base = int(register_json_value("fig1", "B", "Mexico baseline strict N20 admissions", boundary_audit, ["daily_count_reconciliation", "base_counts_2018_2022_source_years", "strict_N20"]))
    register_derived_value(
        "fig1", "B", "Mexico strict additions cross-check", reconciliation["strict_N20"] - strict_base,
        f"{mexico_strict:,} expanded strict N20 − {strict_base:,} baseline strict N20",
        [{"path": rel(boundary_audit), "field": "daily_count_reconciliation.base_counts_2018_2022_source_years.strict_N20", "value": strict_base}, {"path": rel(boundary_audit), "field": "daily_count_reconciliation.augmented_counts_source_years_2018_2023.strict_N20", "value": mexico_strict}],
    )
    return {
        "brazil": brazil,
        "brazil_service": service,
        "mexico_frame": mexico_frame,
        "municipality_frame": municipality_frame,
        "mexico_addition_summary": mexico_addition_summary,
        "mexico_raw": mexico_raw,
        "mexico_valid": mexico_valid,
        "mexico_additions": mexico_additions,
        "mexico_added_strict": mexico_added_strict,
        "mexico_expanded": mexico_expanded,
        "mexico_strict": mexico_strict,
        "eligibility_difference": eligible_difference,
        "municipality_missing": municipality_missing,
        "service_days": service_days,
        "municipality_complete": municipality_frame["exposure_complete_n"],
        "municipality_count": municipality_n,
    }


def make_fig1() -> plt.Figure:
    d = register_fig1_data()
    fig = plt.figure(figsize=(7.2, 4.8))
    add_title(fig, "Study samples and analysis frames", "Brazil and Mexico urolithiasis admissions")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    heading(fig, "A", "Brazil • separate analysis frames", 0.05, 0.86)
    heading(fig, "B", "Mexico • state flow with parallel branches", 0.51, 0.86)

    b = d["brazil"]
    box(
        ax, (0.05, 0.58), 0.42, 0.19,
        "Primary cohort • N20–N23\n2008–2022",
        f"{fmt_n(b['valid_geography_n'])} regions • {fmt_n(b['source_or_input_n'])} region-days\n{fmt_n(b['events_broad_or_named'])} aggregate admissions",
        BLUE,
    )
    s = d["brazil_service"]
    ax.add_patch(FancyArrowPatch((0.26, 0.57), (0.26, 0.49), arrowstyle="-|>", mutation_scale=11, color=SLATE, lw=1.2))
    box(
        ax, (0.05, 0.27), 0.42, 0.19,
        "Strict N20 service frame\n2018–2022",
        f"{fmt_n(s['valid_geography_n'])} regions × {fmt_n(d['service_days'])} days\n{fmt_n(s['source_or_input_n'])} region-days • {fmt_n(s['events_urgent'])} urgent",
        TEAL,
    )

    raw = d["mexico_raw"]
    valid = d["mexico_valid"]
    expanded = d["mexico_expanded"]
    strict = d["mexico_strict"]
    box(ax, (0.51, 0.62), 0.44, 0.16, "Raw admission records • 2018–2022", f"{fmt_n(raw)} broad N20–N23 admissions", BLUE)
    ax.add_patch(FancyArrowPatch((0.73, 0.61), (0.73, 0.55), arrowstyle="-|>", mutation_scale=11, color=SLATE, lw=1.2))
    box(ax, (0.51, 0.38), 0.44, 0.16, "Valid state/date frame", f"{fmt_n(valid)} admissions • {fmt_n(d['eligibility_difference'])} excluded", TEAL)
    # The valid state/date frame forks into two distinct analysis branches.
    ax.add_patch(FancyArrowPatch((0.65, 0.37), (0.625, 0.31), arrowstyle="-|>", mutation_scale=11, color=SLATE, lw=1.2))
    ax.add_patch(FancyArrowPatch((0.81, 0.37), (0.845, 0.31), arrowstyle="-|>", mutation_scale=11, color=SLATE, lw=1.2))
    box(
        ax, (0.51, 0.07), 0.205, 0.24,
        "Source-expanded",
        f"Broad {fmt_n(expanded)};\n{fmt_n(d['mexico_additions'])} added\nStrict N20 {fmt_n(strict)};\n{fmt_n(d['mexico_added_strict'])} added",
        RUST,
    )
    box(
        ax, (0.745, 0.07), 0.205, 0.24,
        "Municipality\nframe",
        f"{fmt_n(d['municipality_complete'])} complete\n{fmt_n(d['municipality_missing'])} without exposure\n{fmt_n(d['municipality_count'])} municipalities",
        GOLD,
    )
    return fig


def make_fig2() -> plt.Figure:
    same = TABLES / "same_caliber_continuous.csv"
    era = TABLES / "mexico_era5_results.csv"
    strict = REVISION / "nonlinear" / "temperature_and_category_contrasts.csv"
    specs = [
        ("fig2", "A", "Brazil", same, {"country": "Brazil", "model": "continuous_lag07"}, ["country", "model", "RR", "CI_low", "CI_high", "p", "n"], BLUE),
        ("fig2", "A", "Mexico • ERA5-Land", era, {"exposure": "ERA5_Land_state"}, ["exposure", "RR", "CI_low", "CI_high", "p", "n"], TEAL),
        ("fig2", "A", "Mexico • NASA POWER", same, {"country": "Mexico", "model": "continuous_lag07"}, ["country", "model", "RR", "CI_low", "CI_high", "p", "n"], SLATE),
    ]
    broad: list[dict[str, Any]] = []
    for _, panel, label, path, selector, columns, color in specs:
        row = register_row("fig2", panel, label, path, selector, columns)
        row.update(plot_label=label, plot_color=color)
        broad.append(row)
    strict_rows: list[dict[str, Any]] = []
    strict_flow = REVISION / "nonlinear" / "sample_flow.csv"
    for model, label, color in [
        ("brazil_strict_N20_linear", "Brazil", BLUE),
        ("mexico_strict_N20_linear", "Mexico", TEAL),
    ]:
        row = register_row(
            "fig2", "B", label, strict,
            {"model": model, "estimand": "cumulative_temperature_lag0_7_per_C"},
            ["model", "estimand", "RR", "CI_low", "CI_high", "p_t", "t_df", "status"],
        )
        row.update(plot_label=label, plot_color=color)
        strict_rows.append(row)
        register_row(
            "fig2", "B sample/inference", f"{label} sample and cluster counts", strict_flow,
            {"model": model},
            ["family", "outcome", "source_rows_before_lag", "analysis_rows_after_complete_lag_frame", "source_events_full_window", "model_rows_retained", "events_model_retained", "G_cluster_1_retained", "G_cluster_2_retained", "cluster_1", "cluster_2", "t_df"],
        )

    fig = plt.figure(figsize=(7.2, 3.5))
    add_title(fig, "Continuous-temperature associations", "Rate ratio per 1°C • cumulative lag 0–7 • 95% confidence intervals")
    heading(fig, "A", "Broad N20–N23 • common period", 0.05, 0.80)
    draw_forest(fig, (0.35, 0.51, 0.38, 0.20), broad, xlim=(0.99, 1.06), ticks=[1.00, 1.02, 1.04, 1.06], label_x=0.32, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=4)
    heading(fig, "B", "Strict N20 • holiday-adjusted", 0.05, 0.365)
    draw_forest(fig, (0.35, 0.12, 0.38, 0.17), strict_rows, xlim=(0.99, 1.06), ticks=[1.00, 1.02, 1.04, 1.06], label_x=0.32, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=4)
    return fig


def make_fig3() -> plt.Figure:
    results = REVISION / "service" / "results.csv"
    service_rows = read_csv(results)
    def find_service(selector: dict[str, str], panel_name: str) -> dict[str, str]:
        hits = [(line, row) for line, row in enumerate(service_rows, start=2) if all(row.get(k) == v for k, v in selector.items())]
        if len(hits) != 1:
            raise ValueError(f"Expected one service result row for {selector}; found {len(hits)}")
        line, row = hits[0]
        cols = ["model_id", "model_family", "estimand", "RR_or_ratio", "CI_low", "CI_high", "p", "df", "inference", "contrast_terms", "input_rows", "model_rows", "input_events", "model_events", "formula", "cluster"]
        SOURCE_RECORDS.append({
            "figure": "fig3", "panel": panel_name, "plotted_item": selector.get("model_id", selector.get("outcome", "service result")),
            "source_file": source_meta(results), "csv_line": line, "row_selector": selector,
            "source_columns": [c for c in cols if c in row], "source_values": {c: row[c] for c in cols if c in row},
        })
        return row

    panel_a: list[dict[str, Any]] = []
    for model, label, color, marker in [
        ("broad_urgent_elective_holiday_two_way", "Broad • full", BLUE, "o"),
        ("broad_urgent_elective_holiday_two_way_excluding_2020_2021", "Broad • COVID-excluded", BLUE, "D"),
        ("strict_urgent_elective_holiday_two_way", "Strict N20 • full", TEAL, "o"),
        ("strict_urgent_elective_holiday_two_way_excluding_2020_2021", "Strict N20 • COVID-excluded", TEAL, "D"),
    ]:
        row = find_service({"model_id": model, "estimand": "heat_only_urgent_RR_over_elective_RR"}, "B")
        row.update(plot_label=label, plot_color=color, plot_marker=marker)
        row["RR"] = row["RR_or_ratio"]
        panel_a.append(row)

    panel_b: list[dict[str, Any]] = []
    for model, label, color, marker in [
        ("drainage_continuous_holiday_two_way", "Holiday-adjusted • 2018–2022", BLUE, "o"),
        ("drainage_continuous_holiday_excluding_2020_2021", "Holiday-adjusted • COVID-excluded", TEAL, "D"),
    ]:
        row = find_service({"model_id": model, "estimand": "cumulative_temperature_lag0_7_per_1C"}, "C")
        row.update(plot_label=label, plot_color=color, plot_marker=marker)
        row["RR"] = row["RR_or_ratio"]
        panel_b.append(row)
    for row, short_label in zip(panel_b, ["Full calendar", "COVID-excluded"]):
        row["plot_label"] = short_label

    panel_c: list[dict[str, Any]] = []
    outcome_labels = {
        "urgent_N20": "Urgent admissions",
        "admission_attributed_beddays": "Admission bed-days",
        "reimbursement_2022BRL": "Reimbursement • 2022 BRL",
    }
    for outcome, label, color in [
        ("urgent_N20", "Urgent N20 admissions", BLUE),
        ("admission_attributed_beddays", "Admission-attributed bed-days", TEAL),
        ("reimbursement_2022BRL", "Reimbursement • 2022 BRL", RUST),
    ]:
        for spec, model in [("full", f"service_holiday_two_way_{outcome}"), ("COVID-excluded", f"service_holiday_excluding_2020_2021_{outcome}")]:
            row = find_service({"model_id": model, "estimand": "heat_only_vs_neither"}, "A")
            row.update(plot_label=f"{label} • {spec}", plot_color=color, plot_marker="o" if spec == "full" else "D")
            row["RR"] = row["RR_or_ratio"]
            panel_c.append(row)

    burden = UPGRADE / "burden_occupancy_model_results.csv"
    occupancy: list[dict[str, Any]] = []
    for model, label, color in [
        ("occupancy_boundary_inclusive", "Admission + discharge date counted", BLUE),
        ("occupancy_boundary_exclusive", "Discharge date excluded", TEAL),
    ]:
        row = register_row("fig3", "D", label, burden, {"model": model, "outcome": "occupied", "term": "tmean"}, ["model", "outcome", "term", "outcome_sum", "input_rows", "model_rows", "formula", "cluster"])
        row.update(plot_label=label, plot_color=color)
        occupancy.append(row)

    fig = plt.figure(figsize=(7.2, 6.4))
    add_title(fig, "Urgency and service outcomes for urolithiasis", "Brazil rate ratios and resident-associated hospital-day totals")

    heading(fig, "A", "Brazil service outcomes • heat-only vs neither", 0.05, 0.89)
    service_display = [dict(row, plot_label={"full": "All years", "COVID-excluded": "Excluding 2020–2021"}[row["plot_label"].rsplit(" • ", 1)[-1]]) for row in panel_c]
    draw_forest(fig, (0.49, 0.70, 0.24, 0.16), service_display, xlim=(0.94, 1.12), ticks=[0.95, 1.00, 1.05, 1.10], label_x=0.46, value_x=0.76, label_size=9.0, value_size=9.0, value_digits=3)
    for pair_index, domain in enumerate(["Urgent N20\nadmissions", "Admission-attributed\nhospital-days", "Reimbursement\n2022 BRL"]):
        pair_y = 0.70 + 0.16 * ((5 - 2 * pair_index) / 6)
        fig.text(0.05, pair_y, domain, fontsize=9.0, fontweight="bold", color=INK, va="center", linespacing=1.30)

    heading(fig, "B", "Urgent-to-elective RR ratios", 0.05, 0.62)
    fig.text(0.05, 0.585, "N20–N23 and N20 • all years / excluding 2020–2021", fontsize=9.2, fontweight="bold", color=INK)
    urgency_display = [dict(row, plot_label=row["plot_label"].replace("COVID-excluded", "Excluding 2020–2021").replace("full", "All years")) for row in panel_a]
    draw_forest(fig, (0.38, 0.41, 0.35, 0.13), urgency_display, xlim=(0.90, 1.18), ticks=[0.9, 1.0, 1.1], label_x=0.35, value_x=0.76, label_size=9.0, value_size=9.0, value_digits=3)

    heading(fig, "C", "Recorded drainage", 0.05, 0.34)
    fig.text(0.05, 0.305, "Temperature RR per 1°C • lag 0–7", fontsize=9.2, fontweight="bold", color=INK)
    drainage_display = [dict(row, plot_label=row["plot_label"].replace("Full calendar", "All years").replace("COVID-excluded", "Excluding 2020–2021")) for row in panel_b]
    draw_forest(fig, (0.38, 0.21, 0.35, 0.075), drainage_display, xlim=(0.99, 1.03), ticks=[1.00, 1.02], label_x=0.35, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3)

    heading(fig, "D", "Resident-associated hospital-days", 0.05, 0.16)
    fig.text(0.05, 0.125, "Study-calendar totals under two admission/discharge date conventions", fontsize=9.2, fontweight="bold", color=INK)
    axd = fig.add_axes([0.36, 0.045, 0.58, 0.06])
    totals = [num(row, "outcome_sum") for row in occupancy]
    colors = [row["plot_color"] for row in occupancy]
    y = [1, 0]
    axd.barh(y, totals, color=colors, height=0.44, alpha=0.9)
    axd.set_yticks(y, ["Admission + discharge date counted", "Discharge date excluded"])
    axd.tick_params(axis="y", labelsize=9.0, labelcolor=INK, pad=4, length=0)
    axd.tick_params(axis="x", labelsize=9.0, labelcolor=INK, pad=2, length=2)
    for tick_label in [*axd.get_xticklabels(), *axd.get_yticklabels()]:
        tick_label.set_fontweight("bold")
    axd.set_xlim(0, 1_350_000)
    axd.set_xticks([0, 500_000, 1_000_000], ["0", "500,000", "1,000,000"])
    axd.grid(axis="x", color=GRID, linewidth=0.4, zorder=0)
    axd.set_axisbelow(True)
    axd.spines["left"].set_visible(False)
    axd.spines["bottom"].set_color(RULE)
    axd.spines["bottom"].set_linewidth(0.7)
    for yi, total in zip(y, totals):
        axd.text(total + 18_000, yi, fmt_n(total), va="center", fontsize=9.0, fontweight="bold", color=INK)
    return fig


def make_fig4() -> plt.Figure:
    acute = UPGRADE / "acute_repair_model_results.csv"
    lag_rows: list[dict[str, Any]] = []
    for phenotype, prefix, color in [
        ("Broad N20–N23 urgent", "broad_N20_N23_urgent", BLUE),
        ("Strict N20 urgent", "strict_N20_urgent", TEAL),
    ]:
        for lag in (7, 14, 21):
            contrast = f"cumulative_temperature_lag_0_{lag}_per_C"
            model = f"{prefix}_temperature_lag_0_{lag}"
            row = register_row(
                "fig4", "A", f"{phenotype} • lag 0–{lag}", acute,
                {"model": model, "contrast": contrast},
                ["model", "phenotype", "contrast", "RR", "CI_low", "CI_high", "p", "input_n", "model_n", "input_events", "model_events", "max_lag_days", "cluster"],
            )
            row.update(plot_label=f"{phenotype} • lag 0–{lag}", plot_color=color)
            lag_rows.append(row)

    curve_file = REVISION / "nonlinear" / "temperature_curves_psd.csv"
    curve_specs = [
        ("brazil_strict_N20_spline3df", "Brazil • strict N20", BLUE, "Brazil"),
        ("brazil_primary_spline3df_temp_fire_two_way", "Brazil • primary N20–N23", BLUE, "Brazil_primary_2008_2022"),
        ("mexico_strict_N20_spline3df", "Mexico • strict N20", TEAL, "Mexico"),
    ]
    percentiles_file = ROOT / "data" / "aggregate_results" / "nonlinear" / "knots_and_percentiles.json"
    percentile_data = read_json(percentiles_file)
    curve_data: list[tuple[str, str, str, list[dict[str, str]], float, str]] = []
    for model, label, color, percentile_key in curve_specs:
        rows = register_series(
            "fig4", "B", label, curve_file,
            {"model": model, "estimand": "temperature_curve", "reference": "P50", "status": "ok"},
            ["model", "reference", "exposure_value_C", "RR", "CI_low", "CI_high", "p_t", "t_df", "curve_index", "country_or_family", "covariance_version"],
            order_by="exposure_value_C",
        )
        if len(rows) != 100:
            raise ValueError(f"Expected 100 grid points for {model}; found {len(rows)}")
        p50 = float(percentile_data[percentile_key]["P50_C"])
        curve_data.append((model, label, color, rows, p50, percentile_key))
    tests_file = REVISION / "nonlinear" / "joint_nonlinearity_tests_psd.csv"
    tests: dict[str, dict[str, str]] = {}
    for model, label, _, _ in curve_specs:
        key = model
        row = register_row("fig4", "B test", f"{label} global nonlinearity test", tests_file, {"model": key}, ["model", "family", "restrictions", "rank", "F", "df_num", "df_denom", "p_F", "restriction_covariance_min_eigenvalue", "status", "covariance_version", "raw_test_status", "raw_p_F", "raw_restriction_covariance_min_eigenvalue", "raw_covariance_preserved"])
        tests[model] = row
    raw_curves = REVISION / "nonlinear" / "temperature_curves.csv"
    ref_values = {
        model: float(percentile_data[key]["P50_C"])
        for model, _, _, key in curve_specs
    }
    SOURCE_RECORDS.append({
        "figure": "fig4", "panel": "B", "plotted_item": "Raw covariance sensitivity archive and plotted-model P50 references",
        "source_file": source_meta(raw_curves), "row_selectors": [spec[0] for spec in curve_specs],
        "role": "raw curves retained for supplementary sensitivity; P50 values label main-figure curve references",
        "reference_source_file": source_meta(percentiles_file), "reference_temperatures_C": ref_values,
        "reference_json_pointers": {model: [key, "P50_C"] for model, _, _, key in curve_specs},
    })

    fig = plt.figure(figsize=(7.2, 5.6))
    add_title(fig, "Cumulative temperature effects and nonlinear responses", "Rate ratios • 95% confidence intervals")
    heading(fig, "A", "Brazil urgent admissions • cumulative RR per 1°C", 0.05, 0.87)
    draw_forest(fig, (0.38, 0.63, 0.35, 0.20), lag_rows, xlim=(0.99, 1.05), ticks=[1.00, 1.02, 1.04], label_x=0.35, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3)

    fig.text(0.05, 0.54, "Temperature-spline sensitivities • RR vs model P50", fontsize=10.0, fontweight="bold", color=INK)
    all_low = [num(row, "CI_low") for _, _, _, rows, _, _ in curve_data for row in rows]
    all_high = [num(row, "CI_high") for _, _, _, rows, _, _ in curve_data for row in rows]
    curve_ymin, curve_ymax = min(all_low), max(all_high)
    pad = max((curve_ymax - curve_ymin) * 0.04, 0.02)
    curve_ylim = (min(0.65, curve_ymin - pad), max(1.30, curve_ymax + pad))
    for j, (model, label, color, rows, p50, _) in enumerate(curve_data):
        rect = (0.08 + j * 0.32, 0.10, 0.245, 0.27)
        test = tests[model]
        p_text = f"{float(test['p_F']):.2g}" if test.get("p_F") else "not estimable"
        heading(fig, "BCD"[j], label, rect[0] - 0.01, 0.46, 9.0)
        fig.text(rect[0], 0.422, f"P50 reference: {p50:.2f}°C", fontsize=9.0, fontweight="bold", color=INK, ha="left")
        fig.text(rect[0], 0.393, f"Nonlinearity p={p_text}", fontsize=9.0, fontweight="bold", color=INK, ha="left")
        ax = fig.add_axes(rect)
        x = [num(row, "exposure_value_C") for row in rows]
        rr = [num(row, "RR") for row in rows]
        low = [num(row, "CI_low") for row in rows]
        high = [num(row, "CI_high") for row in rows]
        ax.fill_between(x, low, high, color=color, alpha=0.14, linewidth=0)
        ax.plot(x, rr, color=color, linewidth=1.25, linestyle="--" if model.startswith("brazil_primary") else "-")
        ax.axhline(1, color=INK, linewidth=0.65, linestyle=(0, (2, 2)))
        ax.axvline(p50, color=SLATE, linewidth=0.6, linestyle=(0, (2, 2)))
        ax.set_ylim(*curve_ylim)
        ax.set_xlim(min(x), max(x))
        ax.set_xticks([v for v in mpl.ticker.MaxNLocator(nbins=4).tick_values(min(x), max(x)) if min(x) <= v <= max(x)])
        ax.set_yticks([v for v in mpl.ticker.MaxNLocator(nbins=4).tick_values(*curve_ylim) if curve_ylim[0] <= v <= curve_ylim[1]])
        ax.set_xlabel("Mean temperature (°C)", fontsize=9.0, fontweight="bold", labelpad=3)
        ax.tick_params(axis="both", labelsize=9.0, labelcolor=INK, pad=3, length=3)
        for tick_label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
            tick_label.set_fontweight("bold")
        ax.spines["left"].set_color(RULE)
        ax.spines["bottom"].set_color(RULE)
        ax.spines["left"].set_linewidth(0.55)
        ax.spines["bottom"].set_linewidth(0.55)
        ax.grid(axis="y", color=GRID, linewidth=0.4)
    return fig


def make_fig5() -> plt.Figure:
    nonlinear_contrasts = REVISION / "nonlinear" / "temperature_and_category_contrasts.csv"
    brazil_compound = register_row(
        "fig5", "A", "Brazil primary both-vs-neither compound-day contrast", nonlinear_contrasts,
        {"model": "brazil_primary_same_day_compound_two_way", "estimand": "compound_vs_neither"},
        ["model", "estimand", "RR", "CI_low", "CI_high", "p_t", "t_df", "status", "country_or_family"],
    )
    interaction_file = REVISION / "covariance" / "brazil_compound_interaction_two_way.csv"
    brazil_interaction = register_row(
        "fig5", "A interaction", "Brazil multiplicative heat-smoke interaction contrast", interaction_file,
        {"model": "brazil_primary_same_day_compound_two_way", "estimand": "multiplicative_heat_smoke_interaction"},
        ["model", "estimand", "RR", "CI_low", "CI_high", "p", "df", "contrast", "covariance", "model_rows"],
    )
    mexico_compound_file = TABLES / "mexico_composite_results.csv"
    mexico_compound = register_row(
        "fig5", "A", "Mexico state both-vs-neither compound-day contrast", mexico_compound_file,
        {"model": "three_group", "term": "C(group)[T.composite]"},
        ["model", "term", "RR", "CI_low", "CI_high", "p", "n"],
    )
    brazil_compound.update(plot_label="Brazil • N20–N23 primary, two-way", plot_color=BLUE)
    brazil_compound["plot_marker"] = "o"
    mexico_compound.update(plot_label="Mexico • state-level source", plot_color=TEAL)
    mexico_compound["plot_marker"] = "D"
    panel_a = [brazil_compound, mexico_compound]

    municipality_file = TABLES / "mexico_municipality_strict_results.csv"
    municipal_fire = register_row(
        "fig5", "B", "Mexico municipality same-day fire PM2.5", municipality_file,
        {"model": "municipality_ym_dow", "term": "fire10"},
        ["model", "term", "RR", "CI_low", "CI_high", "p", "n"],
    )
    municipal_fire.update(plot_label="Mexico municipality • fire PM2.5 per 10 µg/m³", plot_color=RUST)

    cov_file = REVISION / "covariance" / "mexico_procedure_covariance_contrasts.csv"
    proc_names = {
        "nephrostomy": "Nephrostomy",
        "ureteral_catheterization": "Ureteral catheterization",
        "urinary_stone_fragmentation_lithotripsy": "Stone fragmentation / lithotripsy",
    }
    proc_rows: dict[str, list[dict[str, Any]]] = {"temperature": [], "fire_pm25": []}
    endpoint_ci_differences: list[float] = []
    for endpoint, label in proc_names.items():
        for exposure, unit_label, color, unit_value in [
            ("temperature", "1 degree C", BLUE, "1°C"),
            ("fire_pm25", "10 ug/m3", RUST, "10 µg/m³"),
        ]:
            row = register_row(
                "fig5", "C", f"Mexico {label} • {exposure} PSD-corrected interval", cov_file,
                {"endpoint": endpoint, "exposure": exposure, "unit": unit_label, "covariance_method": "eigenvalue_PSD_CGM"},
                ["model", "endpoint", "exposure", "unit", "covariance_method", "RR", "CI_low", "CI_high", "p", "t_df", "variance_ratio_PSD_to_original", "q_BH_six_contrasts"],
            )
            original = register_row(
                "fig5", "C original-CI audit", f"Mexico {label} • {exposure} original covariance interval (audit only)", cov_file,
                {"endpoint": endpoint, "exposure": exposure, "unit": unit_label, "covariance_method": "original_CGM"},
                ["RR", "CI_low", "CI_high", "p", "t_df", "variance_ratio_PSD_to_original"],
            )
            endpoint_ci_differences.extend([abs(float(row["CI_low"]) - float(original["CI_low"])), abs(float(row["CI_high"]) - float(original["CI_high"]))])
            row.update(plot_label=label, plot_color=color, exposure_unit=unit_value)
            proc_rows[exposure].append(row)
    max_endpoint_change = max(endpoint_ci_differences)
    register_derived_value(
        "fig5", "C original-CI audit", "Maximum PSD-to-original CI endpoint difference", max_endpoint_change,
        "Maximum absolute difference in lower/upper RR limits across six Mexico procedure contrasts",
        [{"path": rel(cov_file), "method_pair": ["eigenvalue_PSD_CGM", "original_CGM"], "contrasts": 6}],
    )
    diag_file = REVISION / "covariance" / "mexico_procedure_covariance_diagnostics.csv"
    for row in read_csv(diag_file):
        SOURCE_RECORDS.append({
            "figure": "fig5", "panel": "C diagnostics", "plotted_item": f"PSD covariance diagnostic • {row['endpoint']}",
            "source_file": source_meta(diag_file), "row_selector": {"endpoint": row["endpoint"]},
            "source_values": row,
        })

    boundary = UPGRADE / "mexico_boundary_model_results.csv"
    subgroup: list[dict[str, Any]] = []
    for model, label, color in [
        ("strict_N20_TIPO_Q_same_day", "Recorded TIPO=Q", BLUE),
        ("strict_N20_no_TIPO_Q_same_day", "No recorded TIPO=Q", TEAL),
    ]:
        row = register_row(
            "fig5", "D", f"Mexico {label} same-day temperature", boundary,
            {"model": model, "contrast": "same_day_temperature_per_C"},
            ["model", "outcome", "contrast", "RR", "CI_low", "CI_high", "p", "events_full_window", "model_rows", "cluster", "inference"],
        )
        row.update(plot_label=label, plot_color=color)
        subgroup.append(row)

    fig = plt.figure(figsize=(7.2, 7.2))
    py = lambda points: points / (fig.get_figheight() * 72)
    add_title(fig, "Exploratory compound and procedure estimates", "Urolithiasis outcomes • rate ratios with 95% confidence intervals")

    heading(fig, "A", "Compound days • both exposures vs neither", 0.05, 0.88)
    brazil_compound["plot_label"] = "Brazil primary"
    mexico_compound["plot_label"] = "Mexico state"
    draw_forest(fig, (0.35, py(402), 0.38, py(32)), [brazil_compound, mexico_compound], xlim=(0.90, 1.12), ticks=[0.90, 1.00, 1.10], label_x=0.32, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3)
    fig.text(0.05, py(390), "Brazil heat–smoke interaction", fontsize=9.2, fontweight="bold", color=INK, ha="left")
    fig.text(0.05, py(375), f"Multiplicative contrast • p={float(brazil_interaction['p']):.3f}", fontsize=9.2, fontweight="bold", color=INK, ha="left")
    interaction_row = dict(brazil_interaction, plot_label="", plot_color=RUST)
    draw_forest(fig, (0.35, py(346), 0.38, py(12)), [interaction_row], xlim=(0.98, 1.08), ticks=[1.00, 1.04, 1.08], label_x=0.32, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3)

    heading(fig, "B", "Municipal fire PM2.5", 0.05, py(311))
    fig.text(0.05, py(293), "Mexico • broad N20–N23 • same-day RR per 10 µg/m³", fontsize=9.2, fontweight="bold", color=INK)
    municipal_fire["plot_label"] = "Municipality model"
    draw_forest(fig, (0.35, py(253), 0.38, py(12)), [municipal_fire], xlim=(0.70, 1.10), ticks=[0.8, 1.0], label_x=0.32, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3)
    fig.text(0.05, py(281), f"{fmt_n(municipal_fire['n'])} observations analyzed", fontsize=9.0, fontweight="bold", color=INK, ha="left", va="center")

    heading(fig, "C", "Mexico strict N20 procedures • lag 0–7", 0.05, py(218))
    fig.text(0.05, py(200), "Temperature RR per 1°C", fontsize=9.2, fontweight="bold", color=INK)
    temp_rows = [dict(r) for r in proc_rows["temperature"]]
    fire_rows = [dict(r) for r in proc_rows["fire_pm25"]]
    draw_forest(fig, (0.40, py(154), 0.33, py(33)), temp_rows, xlim=(0.90, 1.10), ticks=[0.9, 1.0, 1.1], label_x=0.37, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3)
    fig.text(0.05, py(134), "Fire PM2.5 RR per 10 µg/m³", fontsize=9.2, fontweight="bold", color=INK)
    draw_forest(fig, (0.40, py(89), 0.33, py(33)), fire_rows, xlim=(0.20, 2.0), ticks=[0.2, 0.5, 1.0, 2.0], label_x=0.37, value_x=0.76, label_size=9.2, value_size=9.0, value_digits=3, scale="log")

    heading(fig, "D", "Mexico strict N20 • same-day temperature by operation record", 0.05, py(54))
    for row in subgroup:
        operation_label = {"Recorded TIPO=Q": "Recorded operation", "No recorded TIPO=Q": "No recorded operation"}[row["plot_label"]]
        row["plot_label"] = f"{operation_label} • {fmt_n(row['events_full_window'])} admissions"
    draw_forest(fig, (0.45, py(18), 0.28, py(24)), subgroup, xlim=(0.98, 1.04), ticks=[0.98, 1.00, 1.02, 1.04], label_x=0.43, value_x=0.76, label_size=9.0, value_size=9.0, value_digits=3)
    return fig


def export_figure(fig: plt.Figure, stem: str) -> dict[str, Any]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    def layout_issue(message: str) -> None:
        if os.environ.get("FIGURE_LAYOUT_DIAGNOSTIC") == "1":
            print("LAYOUT_DIAGNOSTIC " + message, flush=True)
        else:
            raise ValueError(message)
    for artist in fig.findobj(mpl.text.Text):
        if artist.get_visible() and artist.get_text().strip():
            artist.set_fontweight("bold")
            if artist.get_gid() == "panel-letter":
                artist.set_fontsize(18.0)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    width,height=fig.canvas.get_width_height()
    text_boxes: list[tuple[str, Any, Any]] = []
    for artist in fig.findobj(mpl.text.Text):
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        bbox=artist.get_window_extent(renderer)
        if bbox.x0 < -1 or bbox.y0 < -1 or bbox.x1 > width+1 or bbox.y1 > height+1:
            layout_issue(f"Text outside {stem} canvas: {artist.get_text()!r}")
        if artist.get_gid() == "panel-letter" and artist.get_fontsize() != 18.0:
            raise ValueError(f"Panel letter in {stem} is not exactly 18 pt: {artist.get_text()!r}")
        text_boxes.append((artist.get_text().strip(), artist, bbox))
    for i, (text_a, artist_a, bbox_a) in enumerate(text_boxes):
        for text_b, artist_b, bbox_b in text_boxes[i + 1:]:
            x_overlap = min(bbox_a.x1, bbox_b.x1) - max(bbox_a.x0, bbox_b.x0)
            y_overlap = min(bbox_a.y1, bbox_b.y1) - max(bbox_a.y0, bbox_b.y0)
            if x_overlap <= 1.0 or y_overlap <= 1.0:
                continue
            # Matplotlib can create coincident copies of a tick label for paired axes.
            if text_a == text_b and abs(bbox_a.x0 - bbox_b.x0) < 1 and abs(bbox_a.y0 - bbox_b.y0) < 1:
                continue
            layout_issue(f"Overlapping text in {stem}: {text_a!r} / {text_b!r}")
    paths = {
        "pdf": OUTDIR / f"{stem}.pdf",
        "png": OUTDIR / f"{stem}.png",
        "tiff": OUTDIR / f"{stem}.tiff",
        "editable_svg": OUTDIR / f"{stem}.svg",
    }
    fig.savefig(paths["editable_svg"], format="svg", bbox_inches=None, metadata={"Title": stem, "Creator": "Matplotlib; editable working source"})
    fig.savefig(paths["pdf"], format="pdf", bbox_inches=None, metadata={"Title": stem, "Creator": "Matplotlib; source-mapped revision figures"})
    fig.savefig(paths["png"], format="png", dpi=300, bbox_inches=None)
    fig.savefig(paths["tiff"], format="tiff", dpi=600, bbox_inches=None, pil_kwargs={"compression": "tiff_lzw"})
    return {
        key: {"path": rel(path), "sha256": sha256(path), "bytes": path.stat().st_size}
        for key, path in paths.items()
    }


def build_bundle(figure_ids: list[str]) -> dict[str, Any]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    factories = {"fig1": make_fig1, "fig2": make_fig2, "fig3": make_fig3, "fig4": make_fig4, "fig5": make_fig5}
    outputs: dict[str, Any] = {}
    for figure_id in figure_ids:
        fig = factories[figure_id]()
        outputs[figure_id] = export_figure(fig, figure_id)
        plt.close(fig)
    return outputs


def write_source_map(outputs: dict[str, Any]) -> None:
    SOURCES_OUT.parent.mkdir(parents=True, exist_ok=True)
    for record in SOURCE_RECORDS:
        figure,panel=record['figure'],record['panel']
        if figure=='fig1' and panel=='C':record['panel']='B'
        elif figure=='fig3' and panel in ['A','B','C']:
            record['panel']={'A':'B','B':'C','C':'A'}[panel]
        elif figure=='fig4' and panel!='A':
            label=record.get('plotted_item','')
            if label.startswith('Brazil • strict N20'):record['panel']='B'
            elif label.startswith('Brazil • primary N20'):record['panel']='C'
            elif label.startswith('Mexico • strict N20'):record['panel']='D'
            else:record['panel']='B–D'
    SOURCES_OUT.write_text(
        json.dumps(
            {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "backend": f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro} / Matplotlib {mpl.__version__} Agg",
                "scope": "Exactly five main figures; no model fitting or effect-estimate recalculation.",
                "selected_figures": list(outputs),
                "reuse_level": "Style-only inheritance from an earlier project figure script; all panels and data mappings rebuilt for this revision.",
                "export_contract": {"width_mm": 182.88, "pdf_text": "TrueType font embedding (pdf.fonttype=42)", "png_dpi": 300, "tiff_dpi": 600, "tiff_compression": "LZW", "editable_svg_source": "fig1.svg–fig5.svg adjacent to the raster/vector exports"},
                "source_records": SOURCE_RECORDS,
                "outputs": outputs,
            },
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build source-mapped revision figures.")
    parser.add_argument("--figures", default="fig1,fig2,fig3,fig4,fig5", help="Comma-separated figure IDs; default renders all five.")
    parser.add_argument("--repo-root", type=Path, default=ROOT, help="Package root containing data/figure_inputs.")
    parser.add_argument("--outdir", type=Path, help="Export directory; defaults to <repo-root>/artifacts/figures.")
    args = parser.parse_args()
    selected = [part.strip() for part in args.figures.split(",") if part.strip()]
    allowed = {"fig1", "fig2", "fig3", "fig4", "fig5"}
    if not selected or len(set(selected)) != len(selected) or set(selected) - allowed:
        raise SystemExit(f"Invalid figure selection: {selected}")
    repo_root = args.repo_root.expanduser().resolve()
    outdir = (args.outdir or (repo_root / "artifacts" / "figures")).expanduser().resolve()
    if not (repo_root / "data" / "figure_inputs").is_dir():
        raise SystemExit(f"No data/figure_inputs directory under --repo-root {repo_root}")
    configure_paths(repo_root, outdir)
    for cache in (TABLE_CACHE, JSON_CACHE, FILE_META):
        cache.clear()
    SOURCE_RECORDS.clear()
    outputs = build_bundle(selected)
    if selected == ["fig1", "fig2", "fig3", "fig4", "fig5"] and len(SOURCE_RECORDS) != 73:
        raise RuntimeError(f"Expected 73 source provenance records; found {len(SOURCE_RECORDS)}")
    write_source_map(outputs)
    print(json.dumps({"figures": list(outputs), "output_dir": rel(OUTDIR), "source_map": rel(SOURCES_OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
