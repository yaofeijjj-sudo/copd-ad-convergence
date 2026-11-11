"""Update S4 Table (a) after removing chronic hypersensitivity pneumonitis (CHP) from the lung BH family (2026-10-28).

Values come from 03_candidate_analysis/AxisA_disease_scope_lung_family_20261028/results/lung_family_without_CHP.tsv.
Edits supplementary/S4_Table.xlsx in place (idempotent): CHP FDR cell emptied, IPF and smoking FDR replaced, legend sentence added.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_S4_lung_family.py
"""
from pathlib import Path
import csv
import openpyxl

R = Path.cwd()
X = R / "06_resubmission_20261008/code/supplement/parts/S4_Table.xlsx"
src = R / "03_candidate_analysis/AxisA_disease_scope_lung_family_20261028/results/lung_family_without_CHP.tsv"
with open(src, newline="", encoding="utf-8") as f:
    new = {r["id"]: r for r in csv.DictReader(f, delimiter="\t")}
CONTRAST = {"Chronic hypersensitivity pneumonitis vs control": "GSE150910 CHP", "IPF vs control": "GSE150910 IPF",
            "Ever- vs never-smokers without lung disease": "GSE150910 smoking"}
OLD = "BH FDR is computed separately within the lung and the frontal-cortex disease-scope families."
NEW = ("BH FDR is computed separately within the lung (IPF and ever-smoking) and the frontal-cortex disease-scope families. "
       "Chronic hypersensitivity pneumonitis is shown descriptively without FDR: in GSE150910 all its samples were processed on dates "
       "without controls, so disease and processing batch cannot be separated.")

wb = openpyxl.load_workbook(X)
ws = wb["S4a"]
hdr = [c.value for c in ws[4]]
ci, pi, fi = hdr.index("Contrast"), hdr.index("P"), hdr.index("BH FDR (family)")
done = 0
for row in ws.iter_rows(min_row=5):
    k = CONTRAST.get(row[ci].value)
    if not k:
        continue
    assert abs(float(row[pi].value) - float(new[k]["p"])) < 1e-12, k
    v = new[k]["fdr_without_CHP"]
    row[fi].value = None if v == "NA" else float(v)
    done += 1
assert done == 3
leg = ws["A2"].value
if NEW not in leg:
    assert OLD in leg
    ws["A2"].value = leg.replace(OLD, NEW)
wb.save(X)
print("updated", X.name)
