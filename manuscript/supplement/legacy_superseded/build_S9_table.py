"""Build S9 Table (peripheral-blood Axis A fgsea adjusted for EPIC leukocyte fractions).

Formatting only; values copied from the frozen result files. Layout as S1-S8.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S9_table.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font, Alignment

R = Path.cwd()
OUT = R / "06_resubmission_20261008/code/supplement/parts/S9_Table.xlsx"
D = R / "03_candidate_analysis/AxisA_blood_composition_adjustment_20261016/results"
SAMPLE = {"GSE54837": "whole blood", "GSE248493": "PBMC", "E-MTAB-5278": "whole blood", "GSE63060": "whole blood",
          "GSE63061": "whole blood", "GSE248417": "whole blood"}
DISEASE = {"GSE54837": "COPD", "GSE248493": "COPD", "E-MTAB-5278": "COPD", "GSE63060": "AD", "GSE63061": "AD", "GSE248417": "AD"}


def rows(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return x if x not in ("", "NA") else None


def sheet(wb, name, title, legend, header, data, first=False):
    ws = wb.active if first else wb.create_sheet()
    ws.title = name
    ws["A1"] = title; ws["A1"].font = Font(name="Arial", bold=True, size=11)
    ws["A2"] = legend; ws["A2"].font = Font(name="Arial", size=10)
    ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(len(header), 6))
    ws.row_dimensions[2].height = 100
    for j, h in enumerate(header, 1):
        ws.cell(4, j, h).font = Font(name="Arial", bold=True, size=10)
    for i, r in enumerate(data, 5):
        for j, v in enumerate(r, 1):
            ws.cell(i, j, num(v)).font = Font(name="Arial", size=10)
    for j in range(1, len(header) + 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = 18


wb = openpyxl.Workbook()
fg = rows(D / "fgsea_adjusted_vs_baseline.csv")
sheet(wb, "S9a fgsea adjusted",
      "S9 Table. Peripheral-blood Axis A enrichment adjusted for estimated leukocyte fractions.",
      "(a) Leukocyte fractions were estimated with EPIC 1.1.7 (BRef reference of circulating immune cells; default settings; "
      "microarray input 2^log2 expression, RNA-seq input counts per million). Standardized neutrophil and monocyte fractions "
      "(whole blood) or the monocyte fraction (PBMC) were added to the age- and sex-adjusted limma model; Axis A was tested with "
      "fgseaMultilevel exactly as in Fig 5b, c. Baseline values reproduce Fig 5b, c. (b) EPIC cell fractions, case vs control: "
      "scale(fraction) ~ group + age + sex. (c) Per-sample EPIC cell fractions.",
      ["Dataset", "Disease", "Sample type", "Axis A genes ranked (n)", "Adjusted for", "NES, baseline", "P, baseline",
       "NES, adjusted", "P, adjusted"],
      [[r["dataset_id"], DISEASE[r["dataset_id"]], SAMPLE[r["dataset_id"]], r["axis_present"], r["adjusted_for"],
        r["NES_base"], r["P_base"], r["NES_adj"], r["P_adj"]] for r in fg], first=True)
ft = rows(D / "epic_fraction_case_control.csv")
sheet(wb, "S9b fractions case-control", "S9 Table (b). EPIC fractions in cases vs controls.",
      "Effect in SD units of the fraction, adjusted for age and sex; two-sided P.",
      ["Dataset", "Cell type", "Mean fraction, control", "Mean fraction, case", "Effect (SD)", "P"],
      [[r["dataset_id"], r["cell"], r["mean_control"], r["mean_case"], r["effect_SD"], r["P"]] for r in ft])
ps = rows(D / "epic_fractions_per_sample.csv")
cells = [c for c in ps[0].keys() if c not in ("dataset_id", "sample_id")]
sheet(wb, "S9c per-sample fractions", "S9 Table (c). Per-sample EPIC cell fractions.",
      "otherCells = fraction not attributed to the six reference cell types.",
      ["Dataset", "Sample"] + cells, [[r["dataset_id"], r["sample_id"]] + [r[c] for c in cells] for r in ps])
wb.save(OUT)
w = openpyxl.load_workbook(OUT)
print([(ws.title, ws.max_row - 4) for ws in w.worksheets])
