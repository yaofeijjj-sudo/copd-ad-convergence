"""Apply the 2026-10-27/28 supplement edits that were first made by hand in parts/S*_Table.xlsx.

Runs after build_S1_S13_tables.py, update_S4_lung_family.py and update_S5_mKH.py, and before append_S6g_alveolar_macrophage.py
(run_all.sh gives the order). Values are copied from frozen result files; only titles, legends and labels are written here.
  S4 (a)  GSE47460-GPL6480 row; two-platform sentence in the legend
  S4 (d)  GSE76925 pairs with GSE47460-GPL6480 donors added; platform labels
  S6      (d) per-library QC, (e) flagged-library sensitivity, (f) CellTypist labels; legend of (a)
  S10     (f) annotation-window sensitivity (MAGMA gene-set output); title
  S11     leukocyte-adjustment sheets (b)-(d) removed, duplicate-donor check (e) becomes (b); title and legend
  all     US spelling (standardised, summarised, normalised)
Edits are idempotent. Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_parts_20261027_28.py
"""
from pathlib import Path
import csv
import openpyxl
from openpyxl.styles import Font, Alignment

R = Path.cwd()
PARTS = R / "06_resubmission_20261008/code/supplement/parts"
C = R / "03_candidate_analysis"
QC = C / "Lung_sc_library_QC_check_20261028/results"
WIN = C / "AxisA_MAGMA_window_sensitivity_20261027/results"
F10, F11 = Font(name="Arial", size=10), Font(name="Arial", bold=True, size=10)
US = {"standardised": "standardized", "Standardised": "Standardized", "summarised": "summarized", "normalised": "normalized"}


def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        v = float(x)
        return int(v) if x.strip().lstrip("-").isdigit() else v
    except (TypeError, ValueError):
        return x


def sheet(wb, name, title, legend, header, data, widths=None):
    """(Re)create a sheet in the layout of the hand-made 10-28 sheets: title row 1, legend row 2, headers row 4."""
    if name in wb.sheetnames:
        del wb[name]
    ws = wb.create_sheet(name)
    ws["A1"] = title; ws["A1"].font = Font(name="Arial", bold=True, size=11)
    ws["A2"] = legend; ws["A2"].font = F10; ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    for j, h in enumerate(header, 1):
        ws.cell(4, j, h).font = F11
    for i, r in enumerate(data, 5):
        for j, v in enumerate(r, 1):
            ws.cell(i, j, v).font = F10
    for j, w in enumerate(widths or [], 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = w
    return ws


def write_rows(ws, start, data):
    for i, r in enumerate(data, start):
        for j, v in enumerate(r, 1):
            ws.cell(i, j, v).font = F10


# ---------------- S4 (a) and (d) ----------------
wb = openpyxl.load_workbook(PARTS / "S4_Table.xlsx")
ws = wb["S4a"]
TWO = " GSE47460 was profiled on two platforms (GPL14550 and GPL6480) in different donors; each platform was analyzed separately."
if not ws["A2"].value.endswith(TWO):
    ws["A2"] = ws["A2"].value + TWO
ds = [ws.cell(i, 2).value for i in range(1, ws.max_row + 1)]
assert ds[4] == "GSE47460-GPL14550", ds[4]
g = {r["id"]: r for r in rows(C / "AxisA_singscore_official_rerun_20261004/results/axisA_specificity_results.tsv")}["GSE47460-GPL6480 COPD"]
row6 = ["COPD replication (lung)", "GSE47460-GPL6480", "Agilent microarray (GPL6480)", "Lung tissue", "COPD vs control",
        num(g["n_case"]), num(g["n_control"]), "Unique donors, none shared with GSE47460-GPL14550; replicate arrays averaged within donor",
        "Submitter series matrix; probes mapped to current Entrez IDs, one probe per gene (highest MAD)", num(g["axisA_covered"]),
        "Age, sex, smoking", num(g["effect_SD"]), None, num(g["ci_low"]), num(g["ci_high"]), num(g["p"])]
if "GSE47460-GPL6480" not in ds:
    ws.insert_rows(6)
write_rows(ws, 6, [row6])

ws = wb["S4d"]
ws["A2"] = ("Four GSE76925–GSE47460-GPL14550 pairs share age, sex and FEV1 % predicted; their gene-centered expression correlation ranks "
            "72–231 of 236 GPL14550 donors (one pair also differs in diagnosis). Five GSE76925–GSE47460-GPL6480 pairs share age, sex and "
            "FEV1 % predicted (pre- or post-bronchodilator); their correlation ranks 19–91 of 92 GPL6480 donors. None is the same donor. "
            "doppelgangR flagged no cross-dataset expression duplicates.")
hdr = ["GSE76925 sample", "GSE47460 sample (platform)", "Age", "Sex", "FEV1 % predicted", "Expression correlation",
       "Rank among GSE47460 donors on the same platform", "Best correlation with any other donor"]
d14 = [[r["GSE76925"], f'{r["GSE47460"]} (GPL14550)', num(r["age"]), r["sex"], num(r["fev1_pp"]), num(r["cor"]),
        f'{r["rank_among_GSE47460"]} of {r["n_GSE47460"]}', num(r["best_other_cor"])]
       for r in rows(C / "AxisA_member_direction_robustness_20261018/results/duplicate_lung_fingerprint_expression.tsv")]
d64 = [[r["s76"], f'{r["donor"]} (GPL6480)', num(r["age"]), r["sex"], num(r["fev1"]), num(r["cor"]), f'{r["rank"]} of {r["n"]}',
        num(r["best_other"])]
       for r in rows(C / "AxisA_singscore_official_rerun_20261004/qa/GPL6480_vs_GSE76925_fingerprint_pairs_20261027.tsv")]
assert len(d14) == 4 and len(d64) == 5
ws.delete_rows(4, ws.max_row)
write_rows(ws, 4, [hdr]); [setattr(c, "font", F11) for c in ws[4]]
write_rows(ws, 5, d14 + d64)
wb.save(PARTS / "S4_Table.xlsx")

# ---------------- S6 (d)-(f) ----------------
wb = openpyxl.load_workbook(PARTS / "S6_Table.xlsx")
for s in ("S6d", "S6e", "S6f", "S6g"):  # S6g is re-appended by append_S6g_alveolar_macrophage.py
    if s in wb.sheetnames:
        del wb[s]
TAIL = (" (d) Per-library quality control of the lung data. (e) Lung myeloid test after removing flagged libraries. "
        "(f) CellTypist labels and broad classes.")
if not wb["S6a"]["A2"].value.endswith(TAIL):
    wb["S6a"]["A2"] = wb["S6a"]["A2"].value + TAIL
lq = rows(QC / "library_qc_with_flag.tsv")
assert len(lq) == 65 and sum(r["flagged"] == "TRUE" for r in lq) == 2
sheet(wb, "S6d", "S6 Table (d). Per-library quality control of the lung single-cell data (GSE302339).",
      "Quality control within each 10x library: scuttle quickPerCellQC removed cells more than 3 median absolute deviations below the "
      "library median for total counts or detected genes (log scale) or above it for the mitochondrial percentage; scDblFinder then "
      "removed doublets within each library. Thresholds are on the original scale. A library is flagged when its post-QC median "
      "mitochondrial percentage exceeds the median across the 65 libraries by more than 3 median absolute deviations on the log scale "
      "(cut-off 17.5%); flagged libraries were removed only in the sensitivity analysis in S6 Table (e).",
      ["GEO sample", "Donor", "Disease", "Emphysema group (original study)", "Cells in matrix (n)", "Lower threshold, total counts",
       "Lower threshold, detected genes", "Upper threshold, mitochondrial %", "Cells after QC metrics (n)", "Doublets (n)",
       "Cells retained (n)", "Median mitochondrial % after QC", "Flagged library"],
      [[r["gsm"], r["subject_id"], r["disease_group"], r["emphysema_group"], num(r["n_input"]), num(r["min_umi"]), num(r["min_genes"]),
        num(r["max_mito_pct"]), num(r["n_core_qc"]), num(r["n_doublet"]), num(r["n_pass"]), num(r["median_percent_mito_pass"]),
        {"TRUE": "yes", "FALSE": "no"}[r["flagged"]]] for r in lq],
      [16, 12, 12, 34, 21, 31, 33, 34, 28, 14, 20, 33, 17])
me = rows(QC / "myeloid_donor_test_library_check.tsv")
sheet(wb, "S6e", "S6 Table (e). Lung myeloid donor-level test after removing flagged libraries.",
      "Same aggregation and test as S6 Table (a), on the current Axis A (137 genes): per-library median UCell score of myeloid cells, "
      "median across a donor's libraries, donors with ≥20 cells, two-sided Wilcoxon rank-sum test of COPD vs control. Removed libraries: "
      "GSM9102301 (donor ua72) and GSM9102263 (one of two libraries of donor copd4).",
      ["Analysis", "Donors (n)", "Control donors (n)", "COPD donors (n)", "Myeloid cells (n)", "Median donor UCell, control",
       "Median donor UCell, COPD", "Difference in medians (COPD − control)", "P"],
      [[r["analysis"]] + [num(r[k]) for k in ("donors", "control", "copd", "cells", "median_control", "median_copd", "delta", "P")]
       for r in me])
sheet(wb, "S6f", "S6 Table (f). CellTypist labels of the lung single-cell data and their broad classes.",
      "Labels from the CellTypist Human_Lung_Atlas model (CellTypist 1.7.1, best-match mode, no majority voting) for the 211,464 "
      "QC-passed cells, and the broad class used in the analyses. Other/Progenitor (hematopoietic stem cells) was not analyzed as a class.",
      ["Broad class", "CellTypist label", "Cells (n)"],
      [[r["broad_celltype"], r["celltypist_label"], num(r["n_cells"])] for r in rows(QC / "celltypist_label_to_broad.tsv")])
wb.save(PARTS / "S6_Table.xlsx")

# ---------------- S10 (f) ----------------
wb = openpyxl.load_workbook(PARTS / "S10_Table.xlsx")
wb["S10a"]["A1"] = ("S10 Table. Sensitivity analyses of MAGMA enrichment in the AD GWAS (removal of the APOE region, conditioning on "
                    "microglial gene sets, restricted background and annotation windows).")
win = []
for f, lab in [("AD_genic_only", "Genic only (window-analysis input)"), ("AD_10kb", "±10 kb"),
               ("AD_35up10down", "35 kb upstream / 10 kb downstream")]:
    t = [l.split() for l in (WIN / f"{f}.gsa.out").read_text().splitlines() if l.strip() and not l.startswith("#")]
    h = {k: i for i, k in enumerate(t[0])}
    sets = {r[h["VARIABLE"]]: r for r in t[1:]}
    p = {a: float(sets[a][h["P"]]) for a in ("AxisA", "AxisB")}
    lo, hi = sorted(p, key=p.get)
    holm = {lo: 2 * p[lo], hi: max(p[hi], 2 * p[lo])}
    for a in ("AxisA", "AxisB"):
        r = sets[a]
        win.append([lab, f"Axis {a[-1]}", int(r[h["NGENES"]]), float(r[h["BETA"]]), float(r[h["BETA_STD"]]), float(r[h["SE"]]),
                    p[a], min(1.0, holm[a])])
sheet(wb, "S10f", "S10 Table (f). Prespecified annotation-window sensitivity analyses of the MAGMA competitive tests in the AD GWAS.",
      "Gene-set step rerun for the fixed Axis A and Axis B gene sets on the gene-level results of the same AD GWAS (GCST90704648; "
      "13,648,198 SNPs) under genic-only annotation and the two windows prespecified on 2026-07-15. One-sided competitive tests "
      "conditioned on gene size, gene density, sample size, inverse minor allele count and their logarithms; Holm correction across the "
      "two axes within each window. All three windows use the gene-level inputs of the 2026-07-15 protocol, in which the per-SNP sample "
      "size is the effective sample size; the primary analysis (S10 Table a–e; Fig 4c) uses N cases + N controls and excludes duplicated "
      "rsIDs, with identical SNP P values. The genic-only rows are therefore the within-protocol reference for the two windows and differ "
      "slightly from the primary result (Axis A P = 7.12 × 10⁻⁶ vs 7.18 × 10⁻⁶).",
      ["Annotation window", "Axis", "Genes tested (n)", "Beta", "Beta (standardized)", "SE", "P (one-sided)", "Holm P (two axes)"], win)
wb.save(PARTS / "S10_Table.xlsx")

# ---------------- S11: keep CAMERA/scores (a) and the duplicate-donor check ----------------
wb = openpyxl.load_workbook(PARTS / "S11_Table.xlsx")
if "S11e" in wb.sheetnames:
    for s in ("S11b", "S11c", "S11d"):
        del wb[s]
    wb["S11e"].title = "S11b"
ws = wb["S11b"]
assert ws["A1"].value.endswith("Duplicate-donor check between the two AddNeuroMed datasets (GSE63060, GSE63061)."), ws["A1"].value
ws["A1"] = "S11 Table (b). Duplicate-donor check between the two AddNeuroMed datasets (GSE63060, GSE63061)."
wb["S11a"]["A1"] = "S11 Table. CAMERA tests, participant-level effects and duplicate-donor checks for Axis A in peripheral blood."
a2 = wb["S11a"]["A2"].value
OLD = "fgsea (fgseaMultilevel on moderated t) is listed for comparison with earlier versions; it treats genes as independent."
NEW = ("fgsea (fgseaMultilevel on moderated t) NES is shown for direction only; fgsea treats genes as independent, so its P values are "
       "smaller than those of CAMERA, which accounts for inter-gene correlation.")
assert OLD in a2 or NEW in a2
wb["S11a"]["A2"] = a2.replace(OLD, NEW)
wb.save(PARTS / "S11_Table.xlsx")

# ---------------- US spelling in every table part ----------------
for n in range(1, 14):
    p = PARTS / f"S{n}_Table.xlsx"
    wb = openpyxl.load_workbook(p); hit = 0
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row:
                if isinstance(c.value, str) and any(k in c.value for k in US):
                    v = c.value
                    for k, u in US.items():
                        v = v.replace(k, u)
                    c.value = v; hit += 1
    if hit:
        wb.save(p)
print("updated S4, S6, S10, S11 and spelling")
