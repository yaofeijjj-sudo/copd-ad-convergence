"""Set every cell of S2 File and S1 Data to Times New Roman, the font of the main text (2026-11-09).
Size, bold, italic and colour are kept. Run after build_S1_S2_File.py and build_S1_Data.py (see run_all.sh).
"""
import copy
from pathlib import Path
import openpyxl
SUP = Path.cwd() / '06_resubmission_20261008/supplementary'
for name in ('S2_File.xlsx', 'S1_Data.xlsx'):
    f = SUP / name
    wb = openpyxl.load_workbook(f)
    n = 0
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row:
                if c.value is None and not c.has_style: continue
                ft = copy.copy(c.font); ft.name = 'Times New Roman'; c.font = ft; n += 1
    wb.save(f); print(name, 'cells set:', n)
