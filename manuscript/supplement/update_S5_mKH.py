"""Insert the modified Hartung-Knapp (mKH) row into S5 Table (a) after 'REML + Hartung-Knapp' (2026-10-28).

Values from 03_candidate_analysis/AxisA_meta_mKH_20261028/results/pooled_mKH.tsv. Edits the workbook in place (idempotent).
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_S5_mKH.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font

R = Path.cwd()
X = R / "06_resubmission_20261008/code/supplement/parts/S5_Table.xlsx"
with open(R / "03_candidate_analysis/AxisA_meta_mKH_20261028/results/pooled_mKH.tsv", newline="", encoding="utf-8") as f:
    res = {r["analysis"]: r for r in csv.DictReader(f, delimiter="\t")}
m = res["REML + modified Hartung-Knapp (mKH)"]
LABEL = "REML + modified Hartung-Knapp (mKH)"
ADD = (" REML + modified Hartung-Knapp (mKH): the Hartung-Knapp variance scaling factor truncated at 1 (Röver et al. 2015); "
       f"here the factor was {float(m['KH_scaling_q']):.3f}, so the interval equals the unmodified Hartung-Knapp interval.")

wb = openpyxl.load_workbook(X)
ws = wb["S5a"]
labels = [ws.cell(i, 1).value for i in range(1, ws.max_row + 1)]
if LABEL not in labels:
    i = labels.index("REML + Hartung-Knapp") + 1  # 1-based row of the HK line
    ws.insert_rows(i + 1)
    for j, v in enumerate([LABEL, m["effect_SD"], m["ci_low"], m["ci_high"], m["p"], m["I2"]], 1):
        ws.cell(i + 1, j, v if j == 1 else float(v)).font = Font(name="Arial", size=10)
if ADD.strip() not in ws["A2"].value:
    ws["A2"].value = ws["A2"].value + ADD
wb.save(X)
print("updated", X.name)
