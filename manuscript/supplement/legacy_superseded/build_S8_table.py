"""Build S8 Table (MAGMA sensitivity analyses: APOE-region exclusion and conditioning on a general microglia set).

Formatting only; values copied from the frozen result files. Layout as S1-S7 (title row 1, legend row 2, headers row 4).
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S8_table.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font, Alignment

R = Path.cwd()
OUT = R / "06_resubmission_20261008/code/supplement/parts/S8_Table.xlsx"
APOE = R / "03_candidate_analysis/AxisA_MAGMA_APOE_exclusion_20261016/results/AD_MAGMA_APOE_exclusion.csv"
COND = R / "03_candidate_analysis/AxisA_MAGMA_conditional_microglia_20261016/results"


def rows(p, delim=","):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


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
a = rows(APOE)
sheet(wb, "S8a APOE exclusion",
      "S8 Table. MAGMA sensitivity analyses of the Axis A and Axis B enrichment in the AD GWAS.",
      "(a) APOE-region exclusion. Region: APOE (Entrez 348) gene boundaries +/- 1 Mb, NCBI build 37.3 "
      "(chr19:44,409,039-46,412,650). 'set_only': region genes removed from the gene sets, background unchanged; "
      "'whole_region': the 69 region genes removed from the whole analysis with MAGMA --settings gene-exclude "
      "(reported in the main text). One-sided competitive P; Holm across the two axes within each variant. "
      "(b) Conditional analysis: MAGMA --model condition=MicrogliaRef, where MicrogliaRef is the top 250 (primary) or all 1,000 "
      "(sensitivity) human microglia markers of BRETIGEA 1.0.3 mapped to Entrez IDs. Rows with model 'conditional' come from "
      "joint models (MODEL column) containing the listed axis and MicrogliaRef. (c) Reference-set composition and overlap.",
      ["Variant", "Axis", "Genes (n)", "Beta", "Beta (standardized)", "SE", "P (one-sided)", "Holm P (two axes)",
       "Region genes in set (n)", "Background genes (n)"],
      [[r["variant"], r["axis"], r["NGENES"], r["BETA"], r["BETA_STD"], r["SE"], r["P_one_sided"],
        r["holm_P_two_axes_within_variant"], r["n_region_genes_in_set"], r["background_genes"]] for r in a], first=True)
c = rows(COND / "conditional_magma_results.csv")
sheet(wb, "S8b conditional", "S8 Table (b). Conditioning on a general microglia gene set.",
      "Unconditional rows test each set alone; conditional rows test the variable jointly with MicrogliaRef (same MODEL id). "
      "AxisA_nonRef = Axis A members outside the reference set (descriptive).",
      ["Reference set", "Model", "Model ID", "Variable", "Genes (n)", "Beta", "Beta (standardized)", "SE", "P (one-sided)",
       "Holm P (two axes)"],
      [[r["reference"], r["model"], r.get("MODEL"), r["VARIABLE"], r["NGENES"], r["BETA"], r["BETA_STD"], r["SE"], r["P"],
        r.get("holm_P_two_axes")] for r in c])
m = rows(COND / "reference_set_manifest.csv")
sheet(wb, "S8c reference set", "S8 Table (c). Microglia reference sets.",
      "Symbols mapped to Entrez IDs with org.Hs.eg.db 3.22.0 (official symbols only, one-to-many dropped).",
      ["Reference set", "Marker symbols (n)", "Mapped Entrez IDs (n)", "In AD gene analysis (n)", "Overlap with Axis A (n)",
       "Overlap with Axis A in analysis (n)", "Overlap with Axis B (n)"],
      [[r["reference"], r["n_symbols"], r["n_entrez"], r["n_in_gene_analysis"], r["overlap_AxisA"],
        r["overlap_AxisA_analysed"], r["overlap_AxisB"]] for r in m])
wb.save(OUT)
w = openpyxl.load_workbook(OUT)
print([(ws.title, ws.max_row - 4) for ws in w.worksheets])
