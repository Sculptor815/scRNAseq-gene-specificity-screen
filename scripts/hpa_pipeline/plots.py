from __future__ import annotations

import csv
import math
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Sequence


LABEL_LEFT_BASE = 42
PLOT_LEFT_BASE = 660
PLOT_RIGHT_BASE = 1715
VALUE_X_BASE = 1750
BASE_WIDTH = 1950
HEADER_HEIGHT = 165
FOOTER_HEIGHT = 58

COLORS = {
    "background": "#FFFFFF",
    "text": "#18212B",
    "muted": "#6B7280",
    "grid": "#DCE3E8",
    "axis": "#AAB4BF",
    "other": "#EF5A5A",
    "megakaryocytic": "#167C9C",
    "hematopoietic": "#8B949E",
    "bone_non_hema": "#E1C45A",
    "bone": "#167C9C",
    "megakaryocytic_bg": "#EAF6FA",
    "hematopoietic_bg": "#F1F3F5",
    "bone_non_hema_bg": "#FFF7D6",
    "bone_bg": "#EAF6FA",
    "na_bg": "#F7F8FA",
}


def _pil_imports():
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise RuntimeError("Plot output requires Pillow; install requirements.txt") from exc
    return Image, ImageDraw, ImageFont


def _load_font(size: int, bold: bool = False, italic: bool = False):
    _, _, ImageFont = _pil_imports()
    windows = r"C:\Windows\Fonts"
    if bold:
        candidates = [
            f"{windows}\\arialbd.ttf",
            f"{windows}\\calibrib.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        ]
    elif italic:
        candidates = [
            f"{windows}\\ariali.ttf",
            f"{windows}\\calibrii.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
        ]
    else:
        candidates = [
            f"{windows}\\arial.ttf",
            f"{windows}\\calibri.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return ImageFont.truetype(candidate, size=size)
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def _text_width(draw, text: str, font) -> int:
    box = draw.textbbox((0, 0), text, font=font)
    return box[2] - box[0]


def _ellipsize(draw, text: str, font, max_width: int) -> str:
    if _text_width(draw, text, font) <= max_width:
        return text
    suffix = "..."
    low, high = 0, len(text)
    while low < high:
        middle = (low + high + 1) // 2
        candidate = text[:middle].rstrip() + suffix
        if _text_width(draw, candidate, font) <= max_width:
            low = middle
        else:
            high = middle - 1
    return text[:low].rstrip() + suffix


def _nice_axis_max(value: float) -> float:
    if not math.isfinite(value) or value <= 0:
        return 1.0
    rough = value * 1.05
    exponent = math.floor(math.log10(rough))
    fraction = rough / (10**exponent)
    nice = 1 if fraction <= 1 else 2 if fraction <= 2 else 2.5 if fraction <= 2.5 else 5 if fraction <= 5 else 10
    return nice * (10**exponent)


def _format_tick(value: float) -> str:
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    if abs(value) >= 10 or float(value).is_integer():
        return f"{value:.0f}"
    return f"{value:.1f}"


def _format_value(value: float) -> str:
    return f"{value:,.1f}" if abs(value) >= 1000 else f"{value:.1f}"


def _row_layout(row_count: int) -> tuple[int, int]:
    if row_count <= 65:
        return 30, 17
    if row_count <= 110:
        return 27, 15
    return 24, 14


def _role_colors(role: str | None) -> tuple[str | None, str]:
    if role == "megakaryocytic":
        return COLORS["megakaryocytic_bg"], COLORS["megakaryocytic"]
    if role == "hematopoietic":
        return COLORS["hematopoietic_bg"], COLORS["hematopoietic"]
    if role == "bone_non_hema":
        return COLORS["bone_non_hema_bg"], COLORS["bone_non_hema"]
    if role == "bone":
        return COLORS["bone_bg"], COLORS["bone"]
    return None, COLORS["other"]


def _draw_legend(draw, items: Sequence[tuple[str, str]], y: int, font) -> None:
    x = LABEL_LEFT_BASE
    for label, color in items:
        draw.rounded_rectangle((x, y + 2, x + 18, y + 16), radius=3, fill=color)
        x += 26
        draw.text((x, y), label, font=font, fill=COLORS["muted"])
        x += _text_width(draw, label, font) + 30


def _render_chart(
    output_path: Path,
    title: str,
    subtitle: str,
    categories: Sequence[str],
    values: Sequence[float | None],
    roles: Sequence[str | None],
    unit: str,
    legend: Sequence[tuple[str, str]],
    source_note: str,
    dpi: int,
    canvas_width: int,
) -> None:
    Image, ImageDraw, _ = _pil_imports()
    scale = canvas_width / BASE_WIDTH
    label_left = int(LABEL_LEFT_BASE * scale)
    plot_left = int(PLOT_LEFT_BASE * scale)
    plot_right = int(PLOT_RIGHT_BASE * scale)
    value_x = int(VALUE_X_BASE * scale)
    row_count = len(categories)
    row_height, label_size = _row_layout(row_count)
    height = HEADER_HEIGHT + row_count * row_height + FOOTER_HEIGHT
    image = Image.new("RGB", (canvas_width, height), COLORS["background"])
    draw = ImageDraw.Draw(image)
    title_font = _load_font(30, bold=True)
    subtitle_font = _load_font(18)
    legend_font = _load_font(16)
    axis_font = _load_font(15)
    axis_bold_font = _load_font(16, bold=True)
    label_font = _load_font(label_size)
    value_font = _load_font(max(13, label_size - 1))
    value_bold_font = _load_font(max(13, label_size - 1), bold=True)
    italic_font = _load_font(max(13, label_size - 1), italic=True)
    footer_font = _load_font(14)

    draw.text((label_left, 22), title, font=title_font, fill=COLORS["text"])
    draw.text((label_left, 61), subtitle, font=subtitle_font, fill=COLORS["muted"])
    _draw_legend(draw, legend, 92, legend_font)
    numeric = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    axis_max = _nice_axis_max(max(numeric, default=0.0))
    ticks = [(axis_max * index / 5, _format_tick(axis_max * index / 5)) for index in range(6)]
    axis_y = 146
    draw.text((plot_left, 112), f"RNA expression ({unit})", font=axis_bold_font, fill=COLORS["text"])
    for tick, label in ticks:
        x = plot_left + int((plot_right - plot_left) * tick / axis_max)
        draw.line((x, axis_y - 5, x, axis_y), fill=COLORS["axis"], width=1)
        draw.text((x - _text_width(draw, label, axis_font) // 2, 123), label, font=axis_font, fill=COLORS["muted"])
    draw.text((value_x, 112), "Value", font=axis_bold_font, fill=COLORS["text"])
    draw.line((plot_left, axis_y, plot_right, axis_y), fill=COLORS["axis"], width=1)

    label_width = plot_left - label_left - 40
    for index, (category, value, role) in enumerate(zip(categories, values, roles)):
        y0 = HEADER_HEIGHT + index * row_height
        y1 = y0 + row_height
        background, bar_color = _role_colors(role)
        if background:
            draw.rectangle((18, y0, canvas_width - 18, y1), fill=background)
        x = label_left
        while x < plot_right:
            draw.line((x, y1 - 1, min(x + 5, plot_right), y1 - 1), fill=COLORS["grid"], width=1)
            x += 10
        display = _ellipsize(draw, str(category), label_font, label_width)
        box = draw.textbbox((0, 0), display, font=label_font)
        label_y = y0 + (row_height - (box[3] - box[1])) // 2 - box[1]
        draw.text((label_left, label_y), display, font=label_font, fill=COLORS["text"])
        if value is None or not math.isfinite(float(value)):
            draw.text((plot_left + 7, label_y), "N/A", font=italic_font, fill=COLORS["muted"])
            draw.text((value_x, label_y), "N/A", font=italic_font, fill=COLORS["muted"])
            continue
        numeric_value = float(value)
        fraction = max(0.0, min(1.0, numeric_value / axis_max))
        bar_end = plot_left + int((plot_right - plot_left) * fraction)
        bar_top = y0 + max(5, row_height // 5)
        bar_bottom = y1 - max(5, row_height // 5)
        if numeric_value == 0:
            middle = (bar_top + bar_bottom) // 2
            draw.ellipse((plot_left - 4, middle - 4, plot_left + 4, middle + 4), fill=bar_color)
        else:
            draw.rounded_rectangle(
                (plot_left, bar_top, max(plot_left + 3, bar_end), bar_bottom),
                radius=3,
                fill=bar_color,
            )
        draw.text(
            (value_x, label_y),
            _format_value(numeric_value),
            font=value_bold_font if role in {"megakaryocytic", "bone"} else value_font,
            fill=COLORS["text"],
        )
    draw.text((label_left, height - 38), source_note, font=footer_font, fill=COLORS["muted"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(output_path, format="PNG", optimize=False, compress_level=6, dpi=(dpi, dpi))


def _safe_name(value: str) -> str:
    return re.sub(r'[<>:"/\\|?*]+', "_", value).strip(" .") or "unnamed"


def _matrix_map(matrix: dict[str, Any]) -> dict[str, list[Any]]:
    return {
        str(row[0]).upper(): list(row)
        for row in matrix["rows"]
        if row and row[0] not in (None, "")
    }


def _prepare_jobs(
    expression: dict[str, Any],
    gene_rows: Sequence[dict[str, Any]],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    sc_headers = [str(value) for value in expression["single_cell"]["headers"][2:]]
    tissue_headers = [str(value) for value in expression["tissues"]["headers"][2:]]
    sc_map = _matrix_map(expression["single_cell"])
    tissue_map = _matrix_map(expression["tissues"])
    megakaryocytic = set(config["megakaryocytic_cell_types"])
    hematopoietic = set(config["hematopoietic_cell_types"])
    bone_non_hema = set(config["bone_marrow_non_hematopoietic_cell_types"])
    jobs: list[dict[str, Any]] = []
    for row in gene_rows:
        gene = str(row["Gene"])
        sc_record = sc_map.get(gene.upper())
        tissue_record = tissue_map.get(gene.upper())
        ensembl = (sc_record[1] if sc_record else None) or (tissue_record[1] if tissue_record else None)
        sc_values = list(sc_record[2:]) if sc_record else [None] * len(sc_headers)
        tissue_values = list(tissue_record[2:]) if tissue_record else [None] * len(tissue_headers)
        sc_roles = [
            "megakaryocytic"
            if name in megakaryocytic
            else "hematopoietic"
            if name in hematopoietic
            else "bone_non_hema"
            if name in bone_non_hema
            else None
            for name in sc_headers
        ]
        jobs.append(
            {
                "Gene": gene,
                "Ensembl": ensembl,
                "GEP": row.get("GEP") or "UNGROUPED",
                "GEP_Zscore": row.get("GEP_Zscore"),
                "GEP_Loading": row.get("GEP_Loading"),
                "GEP_Zscore_Rank": row.get("GEP_Zscore_Rank"),
                "sc_headers": sc_headers,
                "sc_values": sc_values,
                "sc_roles": sc_roles,
                "tissue_headers": tissue_headers,
                "tissue_values": tissue_values,
                "tissue_roles": ["bone" if name.lower() == "bone marrow" else None for name in tissue_headers],
                "dpi": int(config.get("plot_dpi", 150)),
                "width": int(config.get("plot_width_pixels", BASE_WIDTH)),
            }
        )
    return jobs


def _render_gene(job: dict[str, Any], output_root_text: str) -> dict[str, Any]:
    output_root = Path(output_root_text)
    gene = str(job["Gene"])
    gep = str(job["GEP"])
    ensembl = job.get("Ensembl")
    gene_dir = output_root / _safe_name(gep) / _safe_name(gene)
    prefix = _safe_name(gene)
    cell_path = gene_dir / f"{prefix}_01_RNA_every_cell_type.png"
    tissue_path = gene_dir / f"{prefix}_02_RNA_every_tissue.png"
    subtitle = f"{gene} | {gep} | {ensembl or 'No exact HPA match'} | Linear scale"
    source_note = f"Source: Human Protein Atlas public API | Gene: {gene} ({ensembl or 'unmatched'})"
    sc_values = [float(value) if isinstance(value, (int, float)) else None for value in job["sc_values"]]
    tissue_values = [float(value) if isinstance(value, (int, float)) else None for value in job["tissue_values"]]
    _render_chart(
        cell_path,
        f"{gene} - RNA across all HPA cell types",
        subtitle,
        job["sc_headers"],
        sc_values,
        job["sc_roles"],
        "nCPM",
        [
            ("Megakaryocyte / progenitor / platelet", COLORS["megakaryocytic"]),
            ("Other hematopoietic cells", COLORS["hematopoietic"]),
            ("Bone-marrow non-hematopoietic", COLORS["bone_non_hema"]),
            ("Other cells", COLORS["other"]),
        ],
        source_note + " | RNA unit: nCPM",
        job["dpi"],
        job["width"],
    )
    _render_chart(
        tissue_path,
        f"{gene} - RNA across all HPA tissues",
        subtitle,
        job["tissue_headers"],
        tissue_values,
        job["tissue_roles"],
        "nTPM",
        [("Bone marrow", COLORS["bone"]), ("Other tissues", COLORS["other"])],
        source_note + " | RNA unit: nTPM",
        job["dpi"],
        job["width"],
    )
    return {
        "GEP": gep,
        "Gene": gene,
        "Ensembl": ensembl,
        "GEP_Zscore": job.get("GEP_Zscore"),
        "GEP_Loading": job.get("GEP_Loading"),
        "GEP_Zscore_Rank": job.get("GEP_Zscore_Rank"),
        "RNA_cell_type_count": len(job["sc_headers"]),
        "RNA_tissue_count": len(job["tissue_headers"]),
        "RNA_every_cell_type_path": str(cell_path.relative_to(output_root)).replace("\\", "/"),
        "RNA_every_tissue_path": str(tissue_path.relative_to(output_root)).replace("\\", "/"),
    }


def _write_manifest(output_root: Path, rows: Sequence[dict[str, Any]]) -> Path:
    fields = [
        "GEP",
        "Gene",
        "Ensembl",
        "GEP_Zscore",
        "GEP_Loading",
        "GEP_Zscore_Rank",
        "RNA_cell_type_count",
        "RNA_tissue_count",
        "RNA_every_cell_type_path",
        "RNA_every_tissue_path",
    ]
    path = output_root / "HPA_RNA_plot_manifest.csv"
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (str(row["GEP"]), str(row["Gene"]))))
    return path


def _write_readme(output_root: Path, count: int) -> None:
    content = f"""# HPA RNA plots

- Genes: {count}
- Images per gene: 2
- Cell-type RNA unit: nCPM
- Tissue RNA unit: nTPM
- Scale: linear

## Cell-type colors

- Blue: megakaryocytes, megakaryocyte progenitors, MEPs, and platelets.
- Gray: other hematopoietic cells.
- Light yellow: configured bone-marrow non-hematopoietic/niche cell types.
- Red: all remaining cell types.

The light-yellow list is editable in `config.json`. HPA's body-wide cell-type matrix can aggregate a named cell type across tissues, so the color is a plotting annotation rather than a new expression metric.
"""
    (output_root / "README_plot_colors.md").write_text(content, encoding="utf-8")


def generate_plots(
    expression: dict[str, Any],
    gene_rows: Sequence[dict[str, Any]],
    output_root: Path,
    config: dict[str, Any],
    workers: int = 1,
) -> Path:
    started = time.time()
    jobs = _prepare_jobs(expression, gene_rows, config)
    output_root.mkdir(parents=True, exist_ok=True)
    print(f"Rendering {len(jobs)} genes / {len(jobs) * 2} RNA PNGs with {workers} worker(s)", flush=True)
    results: list[dict[str, Any]] = []
    if workers <= 1:
        for index, job in enumerate(jobs, start=1):
            results.append(_render_gene(job, str(output_root)))
            if index == 1 or index % 25 == 0 or index == len(jobs):
                print(f"Plot progress: {index}/{len(jobs)} genes", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_render_gene, job, str(output_root)) for job in jobs]
            for index, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                if index == 1 or index % 25 == 0 or index == len(jobs):
                    print(f"Plot progress: {index}/{len(jobs)} genes", flush=True)
    manifest = _write_manifest(output_root, results)
    _write_readme(output_root, len(jobs))
    print(f"PLOTS_DONE {manifest} ({time.time() - started:.1f}s)", flush=True)
    return manifest
