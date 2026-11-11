"""Apply the 2026-10-30 supplement corrections found when the main text was checked against S2 File.

Runs after update_lung_donor_merge_20261029.py (run_all.sh gives the order). Edits are idempotent.
  S3 (a), (b)  GSE48350 co-expression recomputed without GSM1176257, which the replication analyses exclude (S4 Table (b));
               values from 03_candidate_analysis/AxisA_sigQC_GSE48350_dedup_20261030; the legend states the exclusion
  S4 (c)       legend: the pointer "(main text, Methods)" is dropped, as the main text no longer reports these two values
  S12 (a)      legend: Spearman correlation of shared decreases with total out-degree, the value quoted in the main text
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_parts_20261030.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font

R = Path.cwd()
PARTS = R / "06_resubmission_20261008/code/supplement/parts"
DEDUP = R / "03_candidate_analysis/AxisA_sigQC_GSE48350_dedup_20261030/results"
F10 = Font(name="Arial", size=10)


def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(x):
    return None if x in ("", "NA") else float(x)


def legend_height(ws):
    ws.row_dimensions[2].height = max(ws.row_dimensions[2].height or 45, 15 * (len(ws["A2"].value) // 140 + 1))


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r, i = [0.0] * len(v), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    a, b = ranks(x), ranks(y)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    num = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    return num / (sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b)) ** 0.5


# ---- S3: GSE48350 without GSM1176257 ----
wb = openpyxl.load_workbook(PARTS / "S3_Table.xlsx")
ws = wb["S3a"]
new = {r["dataset"]: r for r in rows(DEDUP / "set_level_coherence.tsv")}["GSE48350"]
hdr = [c.value for c in ws[4]]
col = {h: i + 1 for i, h in enumerate(hdr)}
row = next(r for r in range(5, ws.max_row + 1) if ws.cell(r, 1).value == "GSE48350")
for h, k, cast in [("Samples (n)", "n_samples", int), ("Axis A members measured (n)", "axisA_measured", int),
                   ("Coherence (observed)", "autocor_median", float), ("Null mean", "null_mean", float),
                   ("Null 95th percentile", "null_q95", float), ("Null maximum", "null_max", float),
                   ("Empirical P (one-sided)", "empirical_p_one_sided", float)]:
    ws.cell(row, col[h], cast(float(new[k]))).font = F10
assert ws.cell(row, col["Samples (n)"]).value == 40
old_s = "Preprocessing is identical to the replication analyses."
new_s = "Preprocessing is identical to the replication analyses, which exclude GSM1176257 from GSE48350 (S4 Table (b))."
if new_s not in ws["A2"].value:
    assert old_s in ws["A2"].value
    ws["A2"] = ws["A2"].value.replace(old_s, new_s)
    legend_height(ws)

ws = wb["S3b"]
hdr = [c.value for c in ws[4]]
gl = rows(DEDUP / "gene_level_coherence.tsv")
assert len(gl) == 137 and ws.max_row == 4 + 137
before = {ws.cell(r, 1).value for r in range(5, ws.max_row + 1)}
key = {"Entrez ID": ("entrez_id", int), "Symbol": ("symbol", str), "Datasets measured (n)": ("n_measured", int),
       "Datasets with median > 0 (n)": ("n_positive", int), "Median across datasets": ("median_across_datasets", float),
       "Class": ("class", str)}
for i, g in enumerate(gl, 5):
    for j, h in enumerate(hdr, 1):
        k, cast = key.get(h, (h, float))
        v = g[k]
        ws.cell(i, j, None if v in ("", "NA") else cast(v) if cast is not float else fnum(v)).font = F10
assert {ws.cell(r, 1).value for r in range(5, ws.max_row + 1)} == before
wb.save(PARTS / "S3_Table.xlsx")
print("S3: GSE48350 row and member table updated")

# ---- S4 (c): stale pointer to the main text ----
wb = openpyxl.load_workbook(PARTS / "S4_Table.xlsx")
ws = wb["S4c"]
ws["A2"] = ws["A2"].value.replace(" (main text, Methods)", "")
assert "main text" not in ws["A2"].value
wb.save(PARTS / "S4_Table.xlsx")
print("S4 (c): legend pointer removed")

# ---- S12 (a): Spearman of shared decreases with total out-degree ----
wb = openpyxl.load_workbook(PARTS / "S12_Table.xlsx")
ws = wb["S12a"]
hdr = [c.value for c in ws[4]]
c_dec, c_ol, c_ob = (hdr.index(h) for h in ("Decreased in both (n)", "Out-degree, lung", "Out-degree, brain"))
data = [[c.value for c in r] for r in ws.iter_rows(min_row=5) if r[0].value is not None]
assert len(data) == 357
rho = spearman([r[c_dec] for r in data], [r[c_ol] + r[c_ob] for r in data])
assert round(rho, 2) == 0.68, rho
add = f" Across the 357 TFs, the number of shared decreases increases with the total out-degree (lung plus brain; Spearman ρ = {rho:.2f})."
if "total out-degree" not in ws["A2"].value:
    ws["A2"] = ws["A2"].value + add
    legend_height(ws)
wb.save(PARTS / "S12_Table.xlsx")
print(f"S12 (a): Spearman rho = {rho:.4f} added to the legend")
