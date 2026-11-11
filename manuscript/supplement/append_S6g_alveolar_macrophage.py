"""Append S6 Table (g): alveolar macrophage donor-level tests (2026-10-28), and update the S6a title.

S6 Table d-f were added to supplementary/S6_Table.xlsx on 2026-10-28 outside build_S1_S13_tables.py, so rerunning that builder
drops d-g; this script edits the existing workbook in place and is idempotent (an existing S6g is replaced).
Values are copied from 03_candidate_analysis/AxisA_AM_pseudobulk_20261028/results; layout follows Book.add in the builder.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/append_S6g_alveolar_macrophage.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font, Alignment

R = Path.cwd()
X = R / "06_resubmission_20261008/code/supplement/parts/S6_Table.xlsx"
A = R / "03_candidate_analysis/AxisA_AM_pseudobulk_20261028/results"


def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        v = float(x)
        return int(v) if x.strip().lstrip("-").isdigit() else v
    except (TypeError, ValueError):
        return x


pb = {r["analysis"]: r for r in rows(A / "pseudobulk_fry_camera.tsv")}
uc = {r["analysis"]: r for r in rows(A / "donor_ucell_wilcoxon.tsv")}
share = {r["disease_group"]: r for r in rows(A / "AM_share_of_myeloid_summary.tsv")}
assert list(pb) == list(uc)
LABEL = {"primary: AM": "Alveolar macrophages (primary)",
         "sensitivity a: AM, >=50 cells": "Alveolar macrophages, ≥50 cells per donor",
         "sensitivity b: AM, flagged libraries removed": "Alveolar macrophages, 2 flagged libraries removed",
         "sensitivity c: AM family": "Alveolar macrophage family (4 labels)",
         "reference: Myeloid (S6 Table)": "All myeloid cells (reference; S6 Table a, b)"}
header = ["Analysis", "Donors (n)", "COPD donors (n)", "Control donors (n)", "Cells (n)", "Axis A members after filtering (n)",
          "fry direction", "fry P", "CAMERA direction", "CAMERA P", "CAMERA inter-gene correlation",
          "Median donor UCell, control", "Median donor UCell, COPD", "Difference in medians (COPD − control)",
          "Hodges–Lehmann shift (COPD − control)", "95% CI low", "95% CI high", "Wilcoxon P"]
data = []
for k in pb:
    p, u = pb[k], uc[k]
    assert p["n_donors"] == u["n_donors"], k
    data.append([LABEL[k], p["n_donors"], p["n_COPD"], p["n_control"], p["n_cells"], p["axis_after_filter"],
                 p["fry_direction"], p["fry_p"], p["camera_direction"], p["camera_p"], p["camera_corr"],
                 u["median_control"], u["median_COPD"], u["diff_medians"], u["HL_shift"], u["HL_low"], u["HL_high"], u["wilcoxon_p"]])
c, k = share["COPD"], share["Control"]
legend = ("Alveolar macrophages: cells with the CellTypist (Human_Lung_Atlas, best match) label 'Alveolar macrophages' (GSE302339). "
          "Pseudobulk: raw UMI counts summed per donor; edgeR filterByExpr, TMM and dispersion estimation, then fry (self-contained) and, after voom, "
          "CAMERA (competitive, inter-gene correlation estimated); model ~ disease + age + sex + smoking, as in S6 Table (b). "
          "UCell: per-library median score, then median across a donor's libraries; two-sided Wilcoxon rank-sum test with the Hodges–Lehmann "
          "location shift and its 95% CI. Donors with at least 20 cells (50 in the stated sensitivity analysis). Flagged libraries: S6 Table (d). "
          "Alveolar macrophage family: 'Alveolar macrophages', 'Alveolar Mph CCL3+', 'Alveolar Mph proliferating' and 'Alveolar Mph MT-positive'. "
          f"Alveolar macrophages as a share of myeloid cells per donor, median (IQR): COPD {float(c['median_AM_share']):.2f} "
          f"({float(c['q25']):.2f}–{float(c['q75']):.2f}), control {float(k['median_AM_share']):.2f} ({float(k['q25']):.2f}–{float(k['q75']):.2f}).")

wb = openpyxl.load_workbook(X)
assert wb.sheetnames[:6] == ["S6a", "S6b", "S6c", "S6d", "S6e", "S6f"], wb.sheetnames
if "S6g" in wb.sheetnames:
    del wb["S6g"]
wb["S6a"]["A1"] = ("S6 Table. Donor-level tests of Axis A in single-cell data (UCell, pseudobulk, alveolar macrophages, emphysema subtypes "
                   "and lung function) and quality control and annotation of the lung single-cell data.")
ws = wb.create_sheet("S6g")
ws["A1"] = "S6 Table (g). Donor-level tests of Axis A within lung alveolar macrophages."
ws["A1"].font = Font(name="Arial", bold=True, size=11)
ws["A2"] = legend; ws["A2"].font = Font(name="Arial", size=10)
ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=len(header))
ws.row_dimensions[2].height = max(45, 15 * (len(legend) // 140 + 1))
for j, h in enumerate(header, 1):
    ws.cell(4, j, h).font = Font(name="Arial", bold=True, size=10)
for i, r in enumerate(data, 5):
    for j, v in enumerate(r, 1):
        ws.cell(i, j, num(v)).font = Font(name="Arial", size=10)
for j in range(1, len(header) + 1):
    ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = 16
ws.column_dimensions["A"].width = 44
ws.freeze_panes = "A5"
wb.save(X)
print("written", X.name, wb.sheetnames)
