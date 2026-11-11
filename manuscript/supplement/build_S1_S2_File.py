"""Combine the supporting-information parts into the files uploaded to PLOS ONE.

S1 File (PDF): Supplementary Methods (parts/S1_Methods.pdf, from build_S1_File_methods.py; added 10-29), then Fig A
(earlier WGCNA soft-threshold PDFs, stacked here into parts/S1_Fig.pdf) and Fig B (parts/S2_Fig.pdf).
S2 File (xlsx): Table A-M (parts/S1_Table.xlsx ... S13_Table.xlsx), one sheet per table part.
S1 Data is built separately by build_S1_Data.py.

Formatting only: values, styles, merged cells, widths and freeze panes are copied unchanged; the only edits are the
renaming "Sn Table" -> "Table <letter>" and "S1/S2 Fig" -> "Fig A/B in S1 File" inside the copied text.
PLOS rule (journals.plos.org/plosone/s/supporting-information): combine parts into one file and cite elements as
"Table A in S1 File".

Run from the project root:
  _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S1_S2_File.py
"""
import copy
import re
import subprocess
from pathlib import Path

import openpyxl
from openpyxl.styles import Alignment, Font

R = Path(__file__).resolve().parents[3]
PARTS = R / "06_resubmission_20261008/code/supplement/parts"
OUT = R / "06_resubmission_20261008/supplementary"
LETTER = {n: chr(64 + n) for n in range(1, 14)}  # S1 Table -> Table A ... S13 Table -> Table M
FIG = {"1": "A", "2": "B"}


def rename(s):
    s = re.sub(r"S(\d+) Table", lambda m: "Table " + LETTER[int(m.group(1))], s)
    return re.sub(r"S([12]) Fig", lambda m: "Fig " + FIG[m.group(1)] + " in S1 File", s)


# ---- S2 File: all tables in one workbook ----------------------------------------------------------------------------
wb = openpyxl.Workbook()
toc = wb.active
toc.title = "Contents"
toc["A1"] = "S2 File. Supplementary tables (Table A–M)."
toc["A1"].font = Font(bold=True, size=12)
toc["A2"] = "Each sheet holds one table or table part; the title is in row 1, the legend in row 2 and the column headers in row 4."
toc.append([])
toc.append(["Sheet", "Title"])
for c in toc[4]:
    c.font = Font(bold=True)
toc.column_dimensions["A"].width = 14
toc.column_dimensions["B"].width = 120

for n in range(1, 14):
    src = openpyxl.load_workbook(PARTS / f"S{n}_Table.xlsx")
    for ws in src.worksheets:
        part = ws.title[len(f"S{n}"):]  # "a", "b", ...
        name = f"Table {LETTER[n]}" + (f" ({part})" if part != "a" else "")
        dst = wb.create_sheet(name)
        for row in ws.iter_rows():
            for c in row:
                d = dst.cell(row=c.row, column=c.column, value=rename(c.value) if isinstance(c.value, str) else c.value)
                if c.has_style:
                    d.font, d.fill, d.border = copy.copy(c.font), copy.copy(c.fill), copy.copy(c.border)
                    d.alignment, d.number_format, d.protection = copy.copy(c.alignment), c.number_format, copy.copy(c.protection)
        for k, dim in ws.column_dimensions.items():
            dst.column_dimensions[k].width = dim.width
        for k, dim in ws.row_dimensions.items():
            if dim.height:
                dst.row_dimensions[k].height = dim.height
        for m in ws.merged_cells.ranges:
            dst.merge_cells(str(m))
        dst.freeze_panes = ws.freeze_panes
        toc.append([name, dst["A1"].value])
        toc.cell(row=toc.max_row, column=1).hyperlink = f"#'{name}'!A1"
        toc.cell(row=toc.max_row, column=1).style = "Hyperlink"
for row in toc.iter_rows(min_row=5):
    row[1].alignment = Alignment(wrap_text=False)
wb.save(OUT / "S2_File.xlsx")

# ---- S1 File: Fig A = the earlier soft-threshold figures (vector, unchanged) stacked with a/b labels; Fig B = S2_Fig ----
SOFT = [(R / "WGCNA_figures/01_COPD_soft_threshold_connectivity.pdf", "a", "COPD lung network (GSE76925)"),
        (R / "WGCNA_figures/02_AD_soft_threshold_connectivity.pdf", "b", "AD DLPFC network (GSE33000)")]
panels = ",".join(f'("{p}", "{t}", "{n}")' for p, t, n in SOFT)
swift = f"""
import PDFKit
import CoreText
let tmp = URL(fileURLWithPath: "{PARTS}/S1_Fig.pdf")
let parts: [(String, String, String)] = [{panels}]
let pages = parts.map {{ CGPDFDocument(URL(fileURLWithPath: $0.0) as CFURL)!.page(at: 1)! }}
let w = pages[0].getBoxRect(.mediaBox).width, head: CGFloat = 26
let h = pages.reduce(0) {{ $0 + $1.getBoxRect(.mediaBox).height + head }}
var box = CGRect(x: 0, y: 0, width: w, height: h)
let ctx = CGContext(tmp as CFURL, mediaBox: &box, nil)!
ctx.beginPDFPage(nil)
func text(_ s: String, _ size: CGFloat, _ x: CGFloat, _ y: CGFloat) {{
    let font = CTFontCreateWithName("Arial-BoldMT" as CFString, size, nil)
    let line = CTLineCreateWithAttributedString(NSAttributedString(string: s, attributes: [.init(kCTFontAttributeName as String): font]))
    ctx.textPosition = CGPoint(x: x, y: y); CTLineDraw(line, ctx)
}}
var top = h
for (i, pg) in pages.enumerated() {{
    let r = pg.getBoxRect(.mediaBox)
    text(parts[i].1, 14, 12, top - 18); text(parts[i].2, 12, 30, top - 18)
    top -= head + r.height
    ctx.saveGState(); ctx.translateBy(x: 0, y: top); ctx.drawPDFPage(pg); ctx.restoreGState()
}}
ctx.endPDFPage(); ctx.closePDF()
let out = PDFDocument()
for p in ["{PARTS}/S1_Methods.pdf", tmp.path, "{PARTS}/S2_Fig.pdf"] {{
    let d = PDFDocument(url: URL(fileURLWithPath: p))!
    for i in 0..<d.pageCount {{ out.insert(d.page(at: i)!, at: out.pageCount) }}
}}
precondition(out.write(to: URL(fileURLWithPath: "{OUT}/S1_File.pdf")))
print("pages", out.pageCount)
"""
print(subprocess.run(["swift", "-"], input=swift, text=True, capture_output=True, check=True).stdout.strip())
print("S2 File sheets", len(wb.sheetnames) - 1)
