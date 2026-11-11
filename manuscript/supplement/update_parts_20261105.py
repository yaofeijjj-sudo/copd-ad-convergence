"""Add Table D (e) to S4_Table.xlsx: prespecified secondary comparisons within GSE47460 (2026-11-05).

The main text now cites these comparisons in the disease scope paragraph. Values are read unchanged from the frozen
03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/results/axisA_specificity_results.tsv.
Runs after update_parts_20261030.py (run_all.sh gives the order). Idempotent.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_parts_20261105.py
"""
from copy import copy
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Alignment, Font

R = Path.cwd()
PATH = R / "06_resubmission_20261008/code/supplement/parts/S4_Table.xlsx"
SRC = R / "03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/results/axisA_specificity_results.tsv"
IDS = ["GSE47460-GPL14550 IPF", "GSE47460-GPL14550 otherILD", "GSE47460-GPL14550 COPD-IPF",
       "GSE47460-GPL6480 IPF", "GSE47460-GPL6480 otherILD", "GSE47460-GPL6480 COPD-IPF"]
HDR = ["Platform", "Contrast", "Cases (n)", "Comparison group (n)", "Axis A genes measured (of 137)", "Covariates",
       "Effect (SD)", "95% CI low", "95% CI high", "P"]
LABEL = {"IPF/UIP vs control": "IPF/UIP vs control", "other ILD vs control": "Other interstitial lung disease vs control",
         "COPD vs IPF/UIP": "COPD vs IPF/UIP"}

with open(SRC, newline="", encoding="utf-8") as f:
    by = {r["id"]: r for r in csv.DictReader(f, delimiter="\t")}
wb = openpyxl.load_workbook(PATH)
if "S4e" in wb.sheetnames:
    del wb["S4e"]
ref = wb["S4d"]
ws = wb.create_sheet("S4e")
ws["A1"] = "S4 Table (e). Secondary comparisons within GSE47460 (IPF/UIP, other interstitial lung disease, and COPD vs IPF/UIP)."
ws["A1"].font = copy(ref["A1"].font)
ws["A2"] = ("Prespecified descriptive comparisons on both GSE47460 platforms with the disease-association model of Axis A "
            "(age, sex and smoking). Effects are standardized differences in the Axis A score. P values are unadjusted, "
            "and these comparisons are outside the BH-corrected disease scope families.")
ws["A2"].font = copy(ref["A2"].font)
ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=len(HDR))
ws.row_dimensions[2].height = 45
for j, h in enumerate(HDR, 1):
    c = ws.cell(4, j, h)
    c.font = copy(ref["A4"].font); c.fill = copy(ref["A4"].fill); c.border = copy(ref["A4"].border)
    c.alignment = copy(ref["A4"].alignment)
    ws.column_dimensions[c.column_letter].width = 18
ws.column_dimensions["B"].width = 38
for i, k in enumerate(IDS, 5):
    r = by[k]
    vals = [r["dataset"].split("-")[1], LABEL[r["contrast"]], int(r["n_case"]), int(r["n_control"]), int(r["axisA_covered"]),
            "Age, sex, smoking", float(r["effect_SD"]), float(r["ci_low"]), float(r["ci_high"]), float(r["p"])]
    for j, v in enumerate(vals, 1):
        c = ws.cell(i, j, v); c.font = Font(name="Arial", size=10)
ws.freeze_panes = "A5"
wb.save(PATH)
print("S4e written", len(IDS))
