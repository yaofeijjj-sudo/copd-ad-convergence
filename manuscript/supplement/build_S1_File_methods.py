"""Render the Supplementary Methods of S1 File (parts/S1_File_methods_en.md) to parts/S1_Methods.pdf.

The markdown source holds the complete Methods text moved out of the main manuscript on 2026-10-29; its Chinese master
is 06_resubmission_20261008/补充方法_简体中文.md. Only headings, paragraphs, *italics*, `code` and [links](url) occur, so
a minimal converter is used; the PDF is printed by headless Google Chrome (Times New Roman, 12 pt, US Letter, 1-inch margins, as in the main text since 2026-11-09).
build_S1_S2_File.py then places these pages before Fig A and Fig B.

Run from the project root:
  _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S1_File_methods.py
"""
import html
import re
import subprocess
import tempfile
from pathlib import Path

R = Path(__file__).resolve().parents[3]
PARTS = R / "06_resubmission_20261008/code/supplement/parts"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = s.replace("₁", "<sub>1</sub>")  # Times New Roman has no U+2081 glyph; use a real subscript
    return re.sub(r"\*([^*\s][^*]*)\*", r"<i>\1</i>", s)


body = []
for block in re.split(r"\n\s*\n", (PARTS / "S1_File_methods_en.md").read_text().strip()):
    block = " ".join(block.split("\n")).strip()
    if block.startswith("# "):
        body.append(f"<h1>S1 File. {inline(block[2:])}</h1>")
    elif block.startswith("### "):
        body.append(f"<h2>{inline(block[4:])}</h2>")
    else:
        body.append(f"<p>{inline(block)}</p>")

page = f"""<!doctype html><html><head><meta charset="utf-8"><style>
@page {{ size: Letter; margin: 2.54cm 2.54cm; }}
body {{ font-family: 'Times New Roman', Times, serif; font-size: 12pt; line-height: 1.5; color: #000; }}
h1 {{ font-size: 14pt; margin: 0 0 12pt; }}
h2 {{ font-size: 12pt; margin: 14pt 0 6pt; page-break-after: avoid; }}
p {{ margin: 0 0 7pt; text-align: justify; }}
code {{ font-family: 'Times New Roman', Times, serif; }}
a {{ color: #000; text-decoration: none; }}
sub {{ font-size: 0.7em; line-height: 0; }}
</style></head><body>{''.join(body)}</body></html>"""

with tempfile.TemporaryDirectory() as tmp:
    src = Path(tmp) / "S1_Methods.html"
    src.write_text(page)
    out = PARTS / "S1_Methods.pdf"
    subprocess.run([CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer", f"--user-data-dir={tmp}/profile",
                    f"--print-to-pdf={out}", src.as_uri()], check=True, capture_output=True)
print("wrote", out, out.stat().st_size, "bytes")
