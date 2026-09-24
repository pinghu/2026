#!/usr/bin/env python3
"""
Convert every sheet in an Excel workbook into a nicely formatted table on its
own PowerPoint slide. Each slide includes a footnote indicating the input
file, output file, and sheet name.

Usage:
    python3 excel_to_pptx.py input.xlsx [-o output.pptx] [--max-rows N]

    input.xlsx     Excel workbook to convert
    -o, --output   Output .pptx file path (default: <input_stem>.pptx)
    --max-rows     Max data rows per slide before splitting into multiple
                    slides (default: 20)
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn

HEADER_FILL = RGBColor(0x1F, 0x4E, 0x79)
HEADER_FONT = RGBColor(0xFF, 0xFF, 0xFF)
ROW_FILL_ALT = RGBColor(0xEA, 0xF1, 0xF9)
ROW_FILL = RGBColor(0xFF, 0xFF, 0xFF)
BORDER_COLOR = "9CB7D9"


def set_cell_border(cell):
    """Add thin borders to a table cell (python-pptx has no direct API)."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    for tag in ("a:lnL", "a:lnR", "a:lnT", "a:lnB"):
        ln = tcPr.makeelement(qn(tag), {"w": "6350", "cap": "flat"})
        fill = ln.makeelement(qn("a:solidFill"), {})
        clr = fill.makeelement(qn("a:srgbClr"), {"val": BORDER_COLOR})
        fill.append(clr)
        ln.append(fill)
        tcPr.append(ln)


def format_value(val):
    if pd.isna(val):
        return ""
    if isinstance(val, float):
        if val == int(val) and abs(val) < 1e15:
            return str(int(val))
        return f"{val:.4g}"
    return str(val)


def add_table_slide(prs, df, sheet_name, input_path, output_path, part_label=""):
    blank_layout = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank_layout)
    slide_w, slide_h = prs.slide_width, prs.slide_height

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.3), Inches(0.2), slide_w - Inches(0.6), Inches(0.5))
    tf = title_box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = f"{sheet_name}{part_label}"
    p.font.size = Pt(20)
    p.font.bold = True
    p.font.color.rgb = RGBColor(0x1F, 0x4E, 0x79)

    n_rows = len(df) + 1  # + header
    n_cols = len(df.columns)

    table_top = Inches(0.85)
    table_left = Inches(0.3)
    table_width = slide_w - Inches(0.6)
    table_height = slide_h - table_top - Inches(0.7)

    graphic_frame = slide.shapes.add_table(n_rows, n_cols, table_left, table_top, table_width, table_height)
    table = graphic_frame.table

    # Column widths: distribute evenly
    col_width = Emu(int(table_width / n_cols))
    for c in range(n_cols):
        table.columns[c].width = col_width

    # Header row
    for c, col_name in enumerate(df.columns):
        cell = table.cell(0, c)
        cell.text = str(col_name)
        cell.fill.solid()
        cell.fill.fore_color.rgb = HEADER_FILL
        cell.vertical_anchor = 1  # middle
        for para in cell.text_frame.paragraphs:
            para.alignment = PP_ALIGN.CENTER
            for run in para.runs:
                run.font.size = Pt(11)
                run.font.bold = True
                run.font.color.rgb = HEADER_FONT
        set_cell_border(cell)

    # Data rows
    font_size = Pt(10) if n_cols <= 10 else Pt(9)
    for r in range(len(df)):
        row_fill = ROW_FILL_ALT if r % 2 == 0 else ROW_FILL
        for c in range(n_cols):
            cell = table.cell(r + 1, c)
            cell.text = format_value(df.iat[r, c])
            cell.fill.solid()
            cell.fill.fore_color.rgb = row_fill
            for para in cell.text_frame.paragraphs:
                para.alignment = PP_ALIGN.CENTER
                for run in para.runs:
                    run.font.size = font_size
                    run.font.color.rgb = RGBColor(0x20, 0x20, 0x20)
            set_cell_border(cell)

    table.first_row = True
    table.horz_banding = False

    # Footnote
    footnote_box = slide.shapes.add_textbox(
        Inches(0.25), slide_h - Inches(0.45), slide_w - Inches(0.5), Inches(0.35)
    )
    ftf = footnote_box.text_frame
    ftf.word_wrap = True
    fp = ftf.paragraphs[0]
    fp.text = f"Input: {input_path.name} | Output: {output_path.name} | Sheet: {sheet_name}"
    fp.font.size = Pt(8)
    fp.font.color.rgb = RGBColor(128, 128, 128)
    fp.alignment = PP_ALIGN.LEFT


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Input .xlsx file")
    parser.add_argument("-o", "--output", help="Output .pptx file path")
    parser.add_argument("--max-rows", type=int, default=20, help="Max data rows per slide (default: 20)")
    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    if not input_path.is_file():
        print(f"Error: '{input_path}' is not a file", file=sys.stderr)
        sys.exit(1)

    output_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else input_path.with_suffix(".pptx")
    )

    xl = pd.ExcelFile(input_path)
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    total_slides = 0
    for sheet_name in xl.sheet_names:
        df = xl.parse(sheet_name)
        if df.empty:
            print(f"Skipping empty sheet: {sheet_name}")
            continue

        n_chunks = max(1, -(-len(df) // args.max_rows))  # ceil div
        for i in range(n_chunks):
            chunk = df.iloc[i * args.max_rows : (i + 1) * args.max_rows]
            part_label = f" (part {i + 1}/{n_chunks})" if n_chunks > 1 else ""
            add_table_slide(prs, chunk, sheet_name, input_path, output_path, part_label)
            total_slides += 1
        print(f"Added: {sheet_name} ({len(df)} rows x {len(df.columns)} cols, {n_chunks} slide(s))")

    prs.save(output_path)
    print(f"\nSaved {total_slides} slide(s) to: {output_path}")


if __name__ == "__main__":
    main()
