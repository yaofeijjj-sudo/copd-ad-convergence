"""Build S6 Table (sigQC co-expression) and S7 Table (pseudobulk + cell-proportion adjustment).

Formatting only: values are copied from the frozen result files of the candidate analyses; no statistics
are computed here. Layout follows S1-S5: title in row 1, legend in row 2, column headers in row 4.
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S6_S7_tables.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font, Alignment

R = Path.cwd()
OUT = R / "06_resubmission_20261008/code/supplement/parts"
SIG = R / "03_candidate_analysis/AxisA_sigQC_coherence_20261014/results"
PB = R / "03_candidate_analysis/AxisA_pseudobulk_disease_20261014/results"
DC = R / "03_candidate_analysis/Deconvolution_composition_adjustment_20261002/results"


def tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return x


def sheet(wb, name, title, legend, header, rows, first=False):
    ws = wb.active if first else wb.create_sheet()
    ws.title = name
    ws["A1"] = title; ws["A1"].font = Font(name="Arial", bold=True, size=11)
    ws["A2"] = legend; ws["A2"].font = Font(name="Arial", size=10)
    ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(len(header), 6))
    ws.row_dimensions[2].height = 110
    for j, h in enumerate(header, 1):
        c = ws.cell(4, j, h); c.font = Font(name="Arial", bold=True, size=10)
    for i, r in enumerate(rows, 5):
        for j, v in enumerate(r, 1):
            ws.cell(i, j, num(v)).font = Font(name="Arial", size=10)
    for j in range(1, len(header) + 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = 18
    return ws


# ---------------- S6 Table ----------------
wb = openpyxl.Workbook()
s = tsv(SIG / "set_level_coherence.tsv")
sheet(wb, "S6a set level",
      "S6 Table. Co-expression of Axis A members in eight independent datasets (sigQC).",
      "(a) Set level. Coherence is the median of the pairwise Spearman correlation matrix of measured Axis A members "
      "(sigQC autocor_median), computed without disease labels. The null is 1,000 random gene sets matched to the "
      "measured members on mean expression (20 bins); empirical P (one-sided) = (1 + number of null values >= observed) / 1,001. "
      "Preprocessing is identical to the replication analyses. (b) Member level: each member's median Spearman "
      "correlation with the other members in each dataset. 'consistent' = median > 0 in every dataset where measured "
      "and measured in >= 6 of 8 datasets; labels are descriptive and do not redefine Axis A.",
      ["Dataset", "Tissue", "Samples (n)", "Axis A members measured (n)", "Coherence (observed)",
       "Null mean", "Null 95th percentile", "Null maximum", "Empirical P (one-sided)"],
      [[r["dataset"], r["tissue"], r["n_samples"], r["axisA_measured"], r["autocor_median"], r["null_mean"],
        r["null_q95"], r["null_max"], r["empirical_p_one_sided"]] for r in s], first=True)
g = tsv(SIG / "gene_level_coherence.tsv")
ds = ["GSE47460-GPL14550", "GSE47460-GPL6480", "GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980"]
sheet(wb, "S6b member level", "S6 Table (b). Member-level co-expression.",
      "Median Spearman correlation of each Axis A member with the other members, per dataset (blank = not measured).",
      ["Entrez ID", "Symbol"] + ds + ["Datasets measured (n)", "Datasets with median > 0 (n)", "Median across datasets", "Class"],
      [[r["entrez_id"], r["symbol"]] + [r[d] if r[d] not in ("", "NA") else None for d in ds] +
       [r["n_measured"], r["n_positive"], r["median_across_datasets"], r["class"]] for r in g])
wb.save(OUT / "S6_Table.xlsx")

# ---------------- S7 Table ----------------
wb = openpyxl.Workbook()
pb = tsv(PB / "fry_summary.tsv")
sheet(wb, "S7a pseudobulk",
      "S7 Table. Pseudobulk tests within myeloid cells and cell-proportion adjustment in the discovery cohorts.",
      "(a) Pseudobulk. Raw UMI counts were summed per donor within lung myeloid cells (GSE302339) or brain Microglia 2 "
      "(SEA-AD DFC Micro-PVM_2). edgeR v4: filterByExpr, TMM normalization, estimateDisp; Axis A tested with fry "
      "(self-contained, directional). Primary model per organ is the pre-specified test; the unadjusted lung model is secondary. "
      "Member-level log2 fold changes are descriptive. (b) Disease effect of the Axis A source modules (lung red, brain blue) and "
      "Axis B source modules before and after adding estimated cell proportions to the Fig 2d models (dtangle for lung GSE76925, "
      "BRETIGEA for brain GSE33000). Marker version B (markers overlapping Axis A removed) is primary; A (all markers) and "
      "C (markers in the whole source module removed) and multi-type adjustment are sensitivity analyses. "
      "Percent change = 1 - adjusted beta / baseline beta. (c) Estimated proportions, case vs control, same covariates (SD units).",
      ["Organ", "Model", "Coefficient", "Donors (n)", "Axis A members after filtering (n)", "fry direction", "fry P",
       "Members up (n)", "Members down (n)", "Median member log2FC"],
      [[r["organ"], r["model"], r["coefficient"], r["n_donors"], r["axis_after_filter"], r["fry_direction"],
        r["fry_p"], r["members_up"], r["members_down"], r["median_member_logFC"]] for r in pb], first=True)
am = tsv(DC / "adjusted_models.tsv")
role_order = {"primary": 0, "secondary": 1, "sensitivity": 2, "sensitivity (multi-type)": 3}
am.sort(key=lambda r: (role_order[r["role"]], r["version"], r["organ"], r["module"]))
sheet(wb, "S7b proportion adjustment", "S7 Table (b). Module-disease effects adjusted for estimated cell proportions.",
      "Standardized module eigengene regressed on disease plus the Fig 2d covariates, with and without the listed proportion(s).",
      ["Organ", "Module", "Role", "Method", "Marker version", "Adjusted for", "Baseline beta", "Baseline 95% CI low",
       "Baseline 95% CI high", "Baseline P", "Adjusted beta", "Adjusted 95% CI low", "Adjusted 95% CI high",
       "Adjusted P", "Percent change"],
      [[r["organ"], r["module"], r["role"], r["method"], r["version"], r["adjusted_for"], r["beta_base"], r["lo_base"],
        r["hi_base"], r["p_base"], r["beta_adj"], r["lo_adj"], r["hi_adj"], r["p_adj"], r["pct_change"]] for r in am])
pc = tsv(DC / "proportion_case_control.tsv")
sheet(wb, "S7c proportions", "S7 Table (c). Estimated cell proportions in cases vs controls.",
      "scale(proportion) ~ disease + covariates of Fig 2d. Lung proportions from dtangle are relative; brain values from "
      "BRETIGEA are surrogate proportions (no mean proportion reported).",
      ["Organ", "Marker version", "Cell type", "Disease effect (SD)", "P", "Mean proportion, control", "Mean proportion, case"],
      [[r["organ"], r["version"], r["cell"], r["beta_disease_SD"], r["p"], r["mean_prop_control"] or None,
        r["mean_prop_case"] or None] for r in pc])
wb.save(OUT / "S7_Table.xlsx")

for f in ("S6_Table.xlsx", "S7_Table.xlsx"):
    w = openpyxl.load_workbook(OUT / f)
    print(f, [(ws.title, ws.max_row - 4) for ws in w.worksheets])
