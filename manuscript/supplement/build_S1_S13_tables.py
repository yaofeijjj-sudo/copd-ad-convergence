"""Build S1-S13 Table (PLOS ONE numbering of the 10-19 Chinese manuscript, captions at the end of the Chinese manuscript).

Formatting only: every value is copied from a frozen result file (or from the previous, still-valid supplementary table);
the only computation is the DerSimonian-Laird pooled row of S4 Table, which is checked against the adopted result.
Layout: title in row 1, legend in row 2, column headers in row 4; one workbook per table, sheets a, b, c ...
Run from the project root:
    _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S1_S13_tables.py
"""
from pathlib import Path
import csv, gzip, math
import openpyxl
from openpyxl.styles import Font, Alignment

R = Path.cwd()
OUT = R / "06_resubmission_20261008/code/supplement/parts"
OLD = R / "99_deprecated_archive/20261023_supplement_old_numbering"   # previous S1-S9 (10-08/10-16 numbering)
C = R / "03_candidate_analysis"
REF = R / "05_github_release/copd-ad-convergence/validation/reference"
SRC = R / "06_resubmission_20261008/code/main_figures/source_data"
WG = R / "WGCNA_figures/reproduction_materials/analysis/02_brain_wgcna/results/T04_AD_age_outlier_fixed_module_sensitivity_20260724/tables"
BC = C / "AxisA_blood_composition_adjustment_20261016/results"


def rows(p, delim=None):
    p = Path(p)
    op = gzip.open if p.suffix == ".gz" else open
    with op(p, "rt", newline="", encoding="utf-8") as f:
        d = delim or ("\t" if ".tsv" in p.name else ",")
        return list(csv.DictReader(f, delimiter=d))


def num(x):
    if x is None or x == "" or x == "NA":
        return None
    if isinstance(x, (int, float)):
        return x
    try:
        v = float(x)
        return int(v) if x.strip().lstrip("-").isdigit() else v
    except (TypeError, ValueError):
        return x


def old_sheet(fname, idx=0):
    """Return (legend, header, data rows) of a previous table sheet."""
    ws = openpyxl.load_workbook(OLD / fname, read_only=True).worksheets[idx]
    r = [list(x) for x in ws.iter_rows(values_only=True)]
    hdr = [c for c in r[3] if c is not None]
    data = [x[:len(hdr)] for x in r[4:] if any(v is not None for v in x)]
    return r[1][0], hdr, data


# Figure / sub-table references inside legends carried over from the previous numbering
LEGEND_FIXES = [
    ("(as in Fig 1c and 1d)", "(as in Fig 1a and 1b)"),
    ("Display names follow Fig 5;", "Display names follow Fig 3;"),
    ("plus the Fig 2d covariates", "plus the Fig 1e covariates"),
    ("covariates of Fig 2d.", "covariates of Fig 1e."),
    ("(a) Leukocyte fractions were estimated", "(b) Leukocyte fractions were estimated"),
    ("Axis A was tested with fgseaMultilevel exactly as in Fig 5b, c. Baseline values reproduce Fig 5b, c. (b) EPIC cell fractions, case vs control: scale(fraction) ~ group + age + sex. (c) Per-sample EPIC cell fractions.",
     "Axis A was tested with fgseaMultilevel as in S11 Table (a); baseline values reproduce its fgsea columns. (c) EPIC cell fractions, case vs control: scale(fraction) ~ group + age + sex. (d) Per-sample EPIC cell fractions."),
]


class Book:
    def __init__(self, n):
        self.n, self.wb, self.k = n, openpyxl.Workbook(), 0

    def add(self, title, legend, header, data):
        ws = self.wb.active if self.k == 0 else self.wb.create_sheet()
        letter = "abcdefghij"[self.k]; self.k += 1
        ws.title = f"S{self.n}{letter}"
        ws["A1"] = f"S{self.n} Table{'' if letter == 'a' else ' (' + letter + ')'}. {title}" if letter == "a" else f"S{self.n} Table ({letter}). {title}"
        ws["A1"].font = Font(name="Arial", bold=True, size=11)
        for old, new in LEGEND_FIXES:
            legend = legend.replace(old, new)
        ws["A2"] = legend; ws["A2"].font = Font(name="Arial", size=10)
        ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(len(header), 8))
        ws.row_dimensions[2].height = max(45, 15 * (len(legend) // 140 + 1))
        for j, h in enumerate(header, 1):
            ws.cell(4, j, h).font = Font(name="Arial", bold=True, size=10)
        for i, r in enumerate(data, 5):
            for j, v in enumerate(r, 1):
                ws.cell(i, j, num(v)).font = Font(name="Arial", size=10)
        for j in range(1, len(header) + 1):
            ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = 16
        ws.freeze_panes = "A5"
        return self

    def save(self):
        p = OUT / f"S{self.n}_Table.xlsx"; self.wb.save(p); print("written", p.name, [w.title for w in self.wb.worksheets])


def from_dicts(d, cols):
    """cols: list of (header, key)"""
    return [h for h, _ in cols], [[r.get(k) for _, k in cols] for r in d]


# ---------------- S1: module overlap ----------------
b = Book(1)
leg, hdr, dat = old_sheet("S1_Table.xlsx")
b.add("Gene overlap of all 128 lung–brain module pairs, and sensitivity analyses with the grey-module background and with restored age outliers.", leg, hdr, dat)
d = rows(C / "Methods_checks_20261019/results/fisher_128_pairs_background_comparison.tsv")
h, dd = from_dicts(d, [("Lung module", "lung_module"), ("Brain module", "brain_module"),
                       ("Overlap (n), non-grey background", "overlap_nongrey"), ("Odds ratio, non-grey", "odds_ratio_nongrey"),
                       ("P, non-grey", "p_value_nongrey"), ("BH FDR, non-grey", "fdr_nongrey"),
                       ("Overlap (n), grey-inclusive background", "overlap_grey"), ("Lung module genes (n), grey-inclusive", "lung_genes_grey"),
                       ("Brain module genes (n), grey-inclusive", "brain_genes_grey"), ("Background (n), grey-inclusive", "background_grey"),
                       ("Odds ratio, grey-inclusive", "odds_ratio_grey"), ("P, grey-inclusive", "p_value_grey"), ("BH FDR, grey-inclusive", "fdr_grey")])
dd.sort(key=lambda r: float(r[12]))
b.add("Sensitivity of the module-pair overlap to the background gene set.",
      "One-sided Fisher's exact tests for all 128 lung–brain module pairs repeated with the 4,832 genes that passed filtering in both "
      "networks (including grey-module genes) as background, next to the primary analysis with the 3,623 shared non-grey module genes. "
      "BH FDR is computed across the 128 pairs within each background. Rows are ordered by the grey-inclusive FDR.", h, dd)
d = rows(WG / "t04_fixed_module_disease_effect_comparison_20260724.csv")
h, dd = from_dicts(d, [("Brain module", "module"), ("Beta, age-filtered (n = 459)", "beta_current"), ("95% CI low", "ci_low_current"),
                       ("95% CI high", "ci_high_current"), ("BH FDR", "fdr_current"), ("Beta, 7 age outliers restored (n = 466)", "beta_sensitivity"),
                       ("95% CI low ", "ci_low_sensitivity"), ("95% CI high ", "ci_high_sensitivity"), ("BH FDR ", "fdr_sensitivity")])
rk = {r["module"]: (r["non_grey_abs_rank_current"], r["non_grey_abs_rank_sensitivity"]) for r in rows(WG / "t04_figure_source_data_20260724.csv")}
h += ["Rank of |beta| among 16 non-grey modules, age-filtered", "Rank, outliers restored"]
dd = [r + list(rk.get(r[0], (None, None))) for r in dd]
b.add("Brain module–AD effects after restoring the seven age outliers (GSE33000).",
      "The seven donors removed by the within-group 1.5 × IQR age rule were restored (466 donors); module eigengenes were projected "
      "with the frozen loadings and module–AD effects re-estimated with age and sex as covariates. Beta is the eigengene difference "
      "between AD and control; BH FDR across the non-grey modules.", h, dd)
b.save()

# ---------------- S2: gene sets ----------------
leg, hdr, dat = old_sheet("S5_Table.xlsx")
Book(2).add("The fixed Axis A and Axis B gene sets.", leg, hdr, dat).save()

# ---------------- S3: sigQC ----------------
b = Book(3)
leg, hdr, dat = old_sheet("S6_Table.xlsx", 0); b.add("Co-expression of Axis A members in eight independent datasets (sigQC).", leg, hdr, dat)
leg, hdr, dat = old_sheet("S6_Table.xlsx", 1); b.add("Member-level co-expression.", leg, hdr, dat)
b.save()

# ---------------- S4: datasets, replication effects, donor checks ----------------
bc = {r["cohort"]: r for r in rows(C / "AxisA_replication_after_exclusion_20261019/results/by_cohort.tsv")}
ps = {r["analysis"]: r for r in rows(C / "AxisA_replication_after_exclusion_20261019/results/pooled_sensitivity.tsv")}
AD = ["GSE84422", "GSE118553", "GSE313307", "GSE140841", "GSE48350", "GSE36980"]
y = [float(bc[k]["effect_SD"]) for k in AD]; s = [float(bc[k]["se"]) for k in AD]
w = [1 / v ** 2 for v in s]; yf = sum(a * b_ for a, b_ in zip(w, y)) / sum(w)
Q = sum(a * (b_ - yf) ** 2 for a, b_ in zip(w, y)); k = len(y)
tau2 = max(0.0, (Q - (k - 1)) / (sum(w) - sum(a * a for a in w) / sum(w)))
ws_ = [1 / (v ** 2 + tau2) for v in s]; mu = sum(a * b_ for a, b_ in zip(ws_, y)) / sum(ws_); se_mu = math.sqrt(1 / sum(ws_))
I2 = max(0.0, (Q - (k - 1)) / Q) * 100
from math import erf
pQ = None  # chi-square upper tail via regularized gamma
def chi2_sf(x, df):
    # series for regularized upper incomplete gamma (df/2 small integer or half-integer is fine here)
    a, xx = df / 2.0, x / 2.0
    term, ssum, n = 1.0 / a, 1.0 / a, 0
    while term > 1e-15 * ssum:
        n += 1; term *= xx / (a + n); ssum += term
    return 1 - math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * ssum
pQ = chi2_sf(Q, k - 1)
pz = 2 * (1 - 0.5 * (1 + erf(abs(mu / se_mu) / math.sqrt(2))))
assert abs(mu - float(ps["DL (primary)"]["effect_SD"])) < 1e-9 and abs(I2 - float(ps["DL (primary)"]["I2"])) < 1e-6
leg, hdr, dat = old_sheet("S4_Table.xlsx")
leg = leg.replace("(metafor 5.0.1; S1 Fig)", "(metafor 5.0.1; Fig 2b)") + (
    " GSE48350: GSM1176257 (annotated AD) was excluded because its raw CEL file is byte-identical to control GSM350078 in GSE11882 "
    "and its sex-gene expression is male, consistent with that control (S4 Table b).")
hi = {h_: i for i, h_ in enumerate(hdr)}
for r in dat:
    if r[hi["Dataset"]] == "GSE48350":
        b_ = bc["GSE48350"]
        r[hi["Cases (n)"]] = b_["n_case"]; r[hi["Controls (n)"]] = b_["n_control"]
        r[hi["Sample selection"]] = str(r[hi["Sample selection"]]) + "; GSM1176257 excluded (duplicate of control GSM350078, GSE11882)"
        for col, key in [("Effect (SD)", "effect_SD"), ("SE", "se"), ("95% CI low", "ci_low"), ("95% CI high", "ci_high"), ("P", "p")]:
            r[hi[col]] = b_[key]
    if str(r[hi["Dataset"]]).startswith("Pooled, 6 datasets"):
        r[hi["Cases (n)"]] = sum(int(bc[k]["n_case"]) for k in AD); r[hi["Controls (n)"]] = sum(int(bc[k]["n_control"]) for k in AD)
        r[hi["Effect (SD)"]] = mu; r[hi["SE"]] = se_mu; r[hi["95% CI low"]] = mu - 1.959963984540054 * se_mu
        r[hi["95% CI high"]] = mu + 1.959963984540054 * se_mu; r[hi["P"]] = pz
        r[hi["I² (%)"]] = I2; r[hi["τ²"]] = tau2; r[hi["Heterogeneity P"]] = pQ
    if r[hi["Dataset"]] in bc and r[hi["Dataset"]] != "GSE48350":
        assert abs(float(r[hi["Effect (SD)"]]) - float(bc[r[hi["Dataset"]]]["effect_SD"])) < 1e-9, r[hi["Dataset"]]
b = Book(4).add("Datasets, sample selection, covariates, data processing and donor checks for the Axis A replication and disease scope analyses.", leg, hdr, dat)
d = rows(C / "AxisA_member_direction_robustness_20261018/results/adopted_sample_exclusions_20261019.tsv")
h, dd = from_dicts(d, [("Dataset", "cohort"), ("Sample", "sample"), ("Reason", "reason")])
b.add("Samples excluded after donor checks.", "Samples removed from all analyses after the duplicate and sex checks below.", h, dd)
d = rows(C / "AxisA_member_direction_robustness_20261018/results/sex_check_mismatches.tsv")
h, dd = from_dicts(d, [("Dataset", "cohort"), ("Sample", "sample"), ("Disease (1 = case)", "disease"), ("Annotated sex", "annotated"),
                       ("Predicted sex", "predicted"), ("Sex score (Y genes − XIST)", "sex_score"), ("Y-linked genes used (n)", "y_genes"),
                       ("XIST used", "xist_used"), ("Separation", "separation"), ("Informative dataset", "informative")])
b.add("Samples whose annotated sex disagrees with sex-gene expression.",
      "Predicted sex from the mean of Y-linked genes minus XIST (standardized within dataset). Mismatching samples were kept in the "
      "primary analyses; excluding them gave 0.46 SD in GSE47460 and 0.49 SD for the AD pool (main text, Methods). Datasets whose sex "
      "score was not bimodal are marked as not informative.", h, dd)
d = rows(C / "AxisA_member_direction_robustness_20261018/results/duplicate_lung_fingerprint_expression.tsv")
h, dd = from_dicts(d, [("GSE76925 sample", "GSE76925"), ("GSE47460 sample", "GSE47460"), ("Age", "age"), ("Sex", "sex"), ("FEV1 % predicted", "fev1_pp"),
                       ("Expression correlation", "cor"), ("Rank among 236 GSE47460 donors", "rank_among_GSE47460"), ("Best correlation with any other donor", "best_other_cor")])
b.add("Lung discovery–validation pairs with identical age, sex and FEV1.",
      "Four GSE76925–GSE47460 pairs share age, sex and FEV1 % predicted. Their gene-centered expression correlation ranks 72–231 of 236 "
      "candidates (one pair also differs in diagnosis), so none is the same donor. doppelgangR flagged no cross-dataset expression duplicates.", h, dd)
b.save()

# ---------------- S5: replication sensitivity and members ----------------
b = Book(5)
pl = rows(C / "AxisA_replication_after_exclusion_20261019/results/pooled_sensitivity.tsv")
lo = rows(C / "AxisA_replication_after_exclusion_20261019/results/leave_one_cohort_out.tsv")
dd = [[r["analysis"], r["effect_SD"], r["ci_low"], r["ci_high"], r["p"], r["I2"]] for r in pl]
dd += [["Leave out " + r["dropped"], r["effect_SD"], r["ci_low"], r["ci_high"], r["p"], r["I2"]] for r in lo]
b.add("Sensitivity analyses of Axis A replication and member-level disease effects.",
      "Pooled standardized case–control difference in Axis A score across the six AD frontal-cortex datasets. DL, DerSimonian–Laird; "
      "REML + Hartung–Knapp, restricted maximum likelihood with Hartung–Knapp confidence intervals. 'common members (108)' scores only the "
      "Axis A members measurable in all seven independent datasets; 'GSE84422 probable/definite only' removes possible AD from GSE84422; "
      "'GSE313307 non-infected only' removes donors with a terminal respiratory infection (infection is then dropped from the covariates).",
      ["Analysis", "Effect (SD)", "95% CI low", "95% CI high", "P", "I² (%)"], dd)
me = rows(SRC / "Fig3d_member_effects.tsv")
common = {r["entrez_id"] for r in rows(C / "AxisA_member_direction_robustness_20261018/results/common_measurable_members.tsv")}
cls = {"both_up": "up in both", "both_down": "down in both", "opposite": "opposite"}
dd = [[r["entrez_id"], r["symbol"], r["lung"], r["lung_low"], r["lung_high"], r["ad"], r["ad_low"], r["ad_high"], cls[r["class"]],
       r["both_CI_above_0"], "yes" if r["entrez_id"] in common else "no"] for r in me]
assert len(dd) == 137 and sum(x[8] == "up in both" for x in dd) == 110
b.add("Member-level disease effects of the 137 Axis A members.",
      "Each member's expression was standardized within dataset and regressed on disease with the covariates of the score model; the six "
      "AD effects were pooled with the DerSimonian–Laird model. Effects are in SD units. Class: direction in GSE47460 vs the AD pool.",
      ["Entrez Gene ID", "Symbol", "Effect, GSE47460 (SD)", "95% CI low", "95% CI high", "Effect, AD pooled (SD)", "95% CI low ", "95% CI high ",
       "Class", "Both 95% CIs above 0", "Measurable in all 7 datasets"], dd)
lm = rows(C / "AxisA_replication_after_exclusion_20261019/results/leave_one_member_out.tsv")
dd = [[r["entrez_id"], r["symbol"], r["GSE47460"], r["delta_GSE47460"], r["AD_pooled"], r["delta_AD"], r["AD_pooled_p"]] for r in lm]
b.add("Leave-one-member-out Axis A scores.",
      "Axis A was rescored without each member in turn (singscore), and the disease effects re-estimated in GSE47460 and pooled across the "
      "six AD datasets. Change = effect without the member minus the effect with all members.",
      ["Entrez Gene ID", "Removed member", "Effect, GSE47460 (SD)", "Change, GSE47460", "Effect, AD pooled (SD)", "Change, AD pooled", "P, AD pooled"], dd)
b.save()

# ---------------- S6: single-cell donor-level tests ----------------
b = Book(6)
leg, hdr, dat = old_sheet("S3_Table.xlsx")
b.add("Donor-level tests of Axis A in single-cell data (UCell, pseudobulk, emphysema subtypes and lung function).", leg, hdr, dat)
d = rows(C / "SingleCell_pseudobulk_supplements_20261019/results/pseudobulk_fry_camera.tsv")
h, dd = from_dicts(d, [("Organ", "organ"), ("Cell class", "class"), ("Model", "model"), ("Coefficient", "coefficient"), ("Donors (n)", "n_donors"),
                       ("Axis A members after filtering (n)", "axis_after_filter"), ("fry direction", "fry_direction"), ("fry P", "fry_p"),
                       ("CAMERA direction", "camera_direction"), ("CAMERA P", "camera_p"), ("CAMERA inter-gene correlation", "camera_corr"),
                       ("fry BH FDR (organ, primary)", "fry_BH_within_organ"), ("CAMERA BH FDR (organ, primary)", "camera_BH_within_organ")])
b.add("Donor pseudobulk set tests of Axis A within myeloid cells.",
      "Raw UMI counts summed per donor and cell class; edgeR filterByExpr, TMM and dispersion estimation, then fry (self-contained) and, "
      "after voom, CAMERA (competitive, inter-gene correlation estimated). Lung models adjust for age, sex and smoking; brain models for age, "
      "sex, PMI, RIN and study. E-COPD/NE-COPD: emphysematous and non-emphysematous COPD as defined by the original authors; CPS_Local: "
      "SEA-AD local continuous pseudo-progression score of the DLPFC; ADNC: AD neuropathologic change. BH FDR within organ for the primary models.", h, dd)
d = rows(C / "Methods_checks_20261019/results/single_cell_donor_flow.tsv")
h, dd = from_dicts(d, [("Organ", "organ"), ("Analysis", "analysis"), ("Unit set", "unit_set"), ("Cells or nuclei (n)", "nuclei_or_cells"),
                       ("Donors (n)", "donors"), ("Inclusion rule", "rule")])
b.add("Cells, nuclei and donors used in each single-cell analysis.", "Numbers of cells or nuclei and donors after the inclusion rule of each analysis.", h, dd)
b.save()

# ---------------- S7: composition adjustment ----------------
b = Book(7)
CX = C / "CIBERSORTx_composition_adjustment_20261018"
d = rows(CX / "results/adjusted_models_by_cohort.tsv")
h, dd = from_dicts(d, [("Reference version", "version"), ("Dataset", "cohort"), ("Organ", "organ"), ("Model", "model"), ("Adjusted for", "adjusted_for"),
                       ("Effect (SD)", "effect_SD"), ("SE", "se"), ("95% CI low", "ci_low"), ("95% CI high", "ci_high"), ("P", "p"), ("Samples (n)", "n"),
                       ("Base effect (SD)", "base_effect"), ("Attenuation (fraction)", "pct_change")])
b.add("Disease effects of Axis A before and after adjustment for cell proportions, and the estimated proportions.",
      "Version B (primary): signature built after removing all 137 Axis A genes; version A: all genes. 'target proportion' adds the "
      "standardized myeloid (lung) or immune (brain) proportion to the replication model; 'all proportions except reference' adds every "
      "estimated proportion except the named reference class. Attenuation = 1 − adjusted/base effect.", h, dd)
d = rows(CX / "results/adjusted_models_AD_pooled.tsv")
h, dd = from_dicts(d, [("Reference version", "version"), ("Model", "model"), ("Effect (SD)", "effect_SD"), ("SE", "se"), ("95% CI low", "ci_low"),
                       ("95% CI high", "ci_high"), ("P", "p"), ("I² (%)", "I2"), ("Datasets (n)", "k"), ("Samples (n)", "n_total"),
                       ("Base effect (SD)", "base"), ("Attenuation (fraction)", "pct_change")])
b.add("Pooled AD effects (DerSimonian–Laird) before and after proportion adjustment.", "Six AD frontal-cortex datasets; models as in (a).", h, dd)
d = rows(CX / "results/GSE76925_red_module_adjusted.tsv")
h, dd = from_dicts(d, [("Reference version", "version"), ("Model", "model"), ("Beta", "beta"), ("95% CI low", "ci_low"), ("95% CI high", "ci_high"), ("P", "p")])
b.add("Lung discovery cohort: red-module eigengene–COPD effect with proportion adjustment (GSE76925).",
      "Covariate-adjusted model of the red (Axis A source) module eigengene with age, sex and pack-years, with the estimated myeloid proportion added.", h, dd)
d = rows(CX / "results/proportion_case_control.tsv")
h, dd = from_dicts(d, [("Reference version", "version"), ("Dataset", "cohort"), ("Cell class", "cell"), ("Mean proportion, control", "mean_control"),
                       ("Mean proportion, case", "mean_case"), ("Case–control difference (SD)", "diff_SD"), ("SE", "se"), ("P", "p"),
                       ("Correlation with Axis A score (r)", "r_with_score")])
b.add("Estimated cell proportions in cases and controls and their correlation with the Axis A score.",
      "Case–control differences in standardized proportions use the covariates of the replication model; r is the Pearson correlation "
      "between the proportion and the standardized Axis A score across samples.", h, dd)
d = rows(CX / "qa/cibersortx_fit_quality.tsv")
h, dd = from_dicts(d, [("Reference version", "version"), ("Dataset", "cohort"), ("Samples (n)", "n"), ("Median correlation", "median_correlation"),
                       ("Median RMSE", "median_RMSE"), ("Samples with deconvolution P < 0.05 (n)", "n_p_lt_0.05"), ("Reference class for 'all proportions'", "reference_cell")])
b.add("CIBERSORTx goodness of fit.",
      "Correlation and RMSE between each observed mixture and its reconstruction from the estimated proportions (CIBERSORTx output), "
      "and the number of samples with deconvolution P < 0.05 (100 permutations).", h, dd)
for idx, t in [(1, "Discovery cohorts: module–disease effects adjusted for proportions estimated with dtangle (lung) and BRETIGEA (brain)."),
               (2, "Discovery cohorts: estimated proportions in cases vs controls (dtangle, BRETIGEA).")]:
    leg, hdr, dat = old_sheet("S7_Table.xlsx", idx); b.add(t, leg, hdr, dat)
b.save()

# ---------------- S8: enrichment ----------------
leg, hdr, dat = old_sheet("S2_Table.xlsx")
Book(8).add("All functional enrichment tests for Axis A and Axis B.", leg, hdr, dat).save()

# ---------------- S9: myeloid signatures ----------------
b = Book(9); MR = C / "AxisA_myeloid_reference_20261019/results"
LAB = {"BRETIGEA_microglia_top250": "Microglia, BRETIGEA top 250", "Patir2019_core_microglia": "Core microglia, Patir et al.",
       "Travaglini2020_lung_macrophage": "Lung macrophage, Travaglini et al.", "GO0006956_complement_activation": "Complement activation (GO:0006956)",
       "GO0006909_phagocytosis": "Phagocytosis (GO:0006909)"}
d = rows(MR / "signature_manifest.tsv")
b.add("Comparison of Axis A with myeloid identity and function signatures.",
      "Signatures and their overlap with the 137 Axis A members (genes mapped to current Entrez IDs).",
      ["Signature", "Genes (n)", "Shared with Axis A (n)"], [[LAB.get(r["signature"], r["signature"]), r["n_genes"], r["overlap_axisA"]] for r in d])
d = rows(MR / "signature_models_by_cohort.tsv")
dd = [[LAB.get(r["signature"], r["signature"]), r["cohort"], r["analysis"], r["effect_SD"], r["se"], r["ci_low"], r["ci_high"], r["p"],
       r["vif_ref"], r["r_with_axisA"]] for r in d]
b.add("Per-dataset signature effects and Axis A effects adjusted for each signature.",
      "'signature alone': disease effect of the signature score; 'Axis A adjusted for signature': Axis A disease effect with the signature "
      "score as an additional covariate; 'shared genes removed from both': the same after removing shared members from both scores. "
      "VIF, variance inflation factor of the signature term; r, correlation of the signature score with the Axis A score.",
      ["Signature", "Dataset", "Analysis", "Effect (SD)", "SE", "95% CI low", "95% CI high", "P", "VIF", "r with Axis A"], dd)
d = rows(MR / "signature_models_AD_pooled.tsv")
dd = [[LAB.get(r["signature"], r["signature"]), r["analysis"], r["effect_SD"], r["se"], r["ci_low"], r["ci_high"], r["p"], r["I2"], r["median_r"], r["max_vif"]] for r in d]
b.add("Pooled AD results (DerSimonian–Laird).", "Six AD frontal-cortex datasets; analyses as in (b). Median r across datasets; maximum VIF.",
      ["Signature", "Analysis", "Effect (SD)", "SE", "95% CI low", "95% CI high", "P", "I² (%)", "Median r with Axis A", "Maximum VIF"], dd)
d = rows(MR / "matched_null_empirical_p.tsv")
h, dd = from_dicts(d, [("Dataset", "cohort"), ("Axis A effect (SD)", "axisA_effect"), ("Null median (SD)", "null_median"),
                       ("Empirical P, matched on expression and myeloid specificity", "empirical_p_one_sided"),
                       ("Empirical P, matched on expression only", "expression_matched_only_p")])
b.add("Axis A disease effect compared with random gene sets matched on expression and myeloid specificity.",
      "1,000 random sets drawn within bins of mean expression (10) and organ-specific myeloid specificity (5) to match the measured Axis A "
      "members; one-sided empirical P = proportion of random effects not smaller than the Axis A effect.", h, dd)
b.save()

# ---------------- S10: MAGMA sensitivity ----------------
b = Book(10)
for idx, t in [(0, "Exclusion of the APOE region."), (1, "Conditioning on BRETIGEA microglia gene sets."), (2, "Microglia reference sets.")]:
    leg, hdr, dat = old_sheet("S8_Table.xlsx", idx); b.add(t if idx else "Sensitivity analyses of MAGMA enrichment in the AD GWAS (removal of the APOE region, conditioning on microglial gene sets and restricted background).", leg, hdr, dat)
d = rows(C / "AxisA_auxiliary_checks_20261019/results/magma_patir_conditional.tsv")
h, dd = from_dicts(d, [("Model", "model"), ("Model ID", "MODEL"), ("Variable", "VARIABLE"), ("Genes (n)", "NGENES"), ("Beta", "BETA"),
                       ("Beta (standardized)", "BETA_STD"), ("SE", "SE"), ("P (one-sided)", "P"), ("Holm P (two axes)", "holm_two_axes")])
b.add("Conditioning on the core microglia signature of Patir et al.",
      "MAGMA gene-set analysis with the Patir core microglia genes as a covariate (condition option); Axis A and Axis B tested in separate models.", h, dd)
d = rows(C / "AxisA_auxiliary_checks_20261019/results/magma_background_sensitivity.tsv")
h, dd = from_dicts(d, [("GWAS", "disease"), ("Background", "background"), ("Variable", "VARIABLE"), ("Genes (n)", "NGENES"), ("Beta", "BETA"),
                       ("Beta (standardized)", "BETA_STD"), ("SE", "SE"), ("P (one-sided)", "P"), ("Holm P (four tests)", "holm_4tests")])
b.add("Background restricted to the 4,832 network-eligible genes.",
      "The four primary tests repeated with the gene-include option limited to genes that passed filtering in both networks; Holm correction across the four tests.", h, dd)
b.save()

# ---------------- S11: peripheral blood ----------------
b = Book(11)
d = rows(C / "AxisA_auxiliary_checks_20261019/results/blood_camera_and_effects.tsv")
h, dd = from_dicts(d, [("Dataset", "dataset_id"), ("Disease", "disease"), ("Sample type", "sample_type"), ("Subjects (n)", "n"), ("Axis A genes measured (n)", "axis_present"),
                       ("CAMERA direction", "camera_direction"), ("CAMERA P", "camera_P"), ("CAMERA BH FDR (6 datasets)", "camera_BH"),
                       ("CAMERA inter-gene correlation", "camera_corr"), ("Axis A score difference (SD)", "score_diff_SD"), ("95% CI low", "ci_low"),
                       ("95% CI high", "ci_high"), ("Score P", "score_P"), ("fgsea NES", "fgsea_NES"), ("fgsea P", "fgsea_p")])
b.add("CAMERA tests, participant-level effects, leukocyte proportion adjustment and duplicate-donor checks for Axis A in peripheral blood.",
      "Each dataset analyzed separately with limma models adjusted for age and sex. CAMERA estimates the inter-gene correlation within the "
      "set; BH across the six datasets. Score difference: standardized singscore Axis A score, case minus control. fgsea (fgseaMultilevel on "
      "moderated t) is listed for comparison with earlier versions; it treats genes as independent.", h, dd)
for idx, t in [(0, "Axis A enrichment adjusted for EPIC-estimated leukocyte fractions."), (1, "EPIC leukocyte fractions in cases vs controls."),
               (2, "Per-sample EPIC leukocyte fractions.")]:
    leg, hdr, dat = old_sheet("S9_Table.xlsx", idx); b.add(t, leg, hdr, dat)
d = rows(C / "AxisA_auxiliary_checks_20261019/results/addneuromed_cross_dataset_matches.tsv")
h, dd = from_dicts(d, [("GSE63060 sample", "GSE63060"), ("Best-matching GSE63061 sample", "GSE63061"), ("Correlation, best", "r_best"),
                       ("Correlation, second best", "r_second"), ("Mutual best match", "mutual"), ("Flagged as duplicate", "flagged")])
b.add("Duplicate-donor check between the two AddNeuroMed datasets (GSE63060, GSE63061).",
      "Gene-centered expression correlation between samples of the two datasets, with age and sex compared for mutual best matches; a pair is "
      "flagged only if it is a mutual best match with matching age and sex and a clear margin over the second-best correlation.", h, dd)
b.save()

# ---------------- S12: CellOracle ranking ----------------
d = rows(C / "CellOracle_ranking_scale_20261019/results/all357_TF_ranking.tsv")
for org in ("lung", "brain"):
    srt = sorted(d, key=lambda r: -float(r[f"mean_abs_delta_{org}"]))
    for i, r in enumerate(srt, 1):
        r[f"rank_{org}"] = i
assert [r["rank_lung"] for r in d if r["TF"] in ("KLF4", "PPARG")] and next(r for r in d if r["TF"] == "MEF2C")["rank_brain"] == 2
h, dd = from_dicts(d, [("Rank (shared decreases)", "rank_shared_down"), ("TF", "TF"), ("Axis A genes evaluable (n)", "n_evaluable"),
                       ("Decreased in both (n)", "shared_down"), ("Increased in both (n)", "shared_up"), ("Opposite direction (n)", "opposite"),
                       ("Changed in lung only (n)", "lung_only_change"), ("Changed in brain only (n)", "brain_only_change"),
                       ("Mean signed delta_X, lung", "mean_signed_delta_lung"), ("Mean signed delta_X, brain", "mean_signed_delta_brain"),
                       ("Mean |delta_X|, lung", "mean_abs_delta_lung"), ("Mean |delta_X|, brain", "mean_abs_delta_brain"),
                       ("Within-organ response rank, lung", "rank_lung"), ("Within-organ response rank, brain", "rank_brain"),
                       ("Response percentile, lung", "response_percentile_lung"), ("Response percentile, brain", "response_percentile_brain"),
                       ("Out-degree, lung", "degree_lung"), ("Out-degree, brain", "degree_brain"),
                       ("Shared decreases per 100 out-degree", "shared_down_per_100_outdegree_mean")])
Book(12).add("Complete virtual knockout ranking of 357 transcription factors.",
             "Each TF was set to zero in the lung alveolar macrophage and brain Microglia 2 models (three propagation rounds); delta_X is the "
             "predicted change on CellOracle's natural-log scale for the 136 Axis A genes evaluable in both models. Within-organ response rank "
             "orders TFs by mean |delta_X| in that organ (1 = largest). Out-degree = targets with non-zero network coefficients.", h, dd).save()

# ---------------- S13: TF activity ----------------
d = rows(C / "SingleCell_pseudobulk_supplements_20261019/results/tf_activity_all.tsv.gz")
h, dd = from_dicts(d, [("Organ", "organ"), ("Cell class", "class"), ("Model", "model"), ("Regulatory prior", "network"), ("TF", "TF"),
                       ("Targets measured (n)", "n_targets_measured"), ("Activity logFC (exposure)", "logFC"), ("Mean activity", "AveExpr"),
                       ("t", "t"), ("P", "P.Value"), ("BH FDR, top-10 CellOracle TFs", "FDR_top10"), ("BH FDR, all TFs", "FDR_all"),
                       ("Correlation with donor Axis A score (r)", "r_with_axisA_score")])
Book(13).add("Transcription factor activity inferred in donor pseudobulk.",
             "Activities from decoupleR's univariate linear model on donor pseudobulk log-CPM, with the full CollecTRI prior or the prior after "
             "removing Axis A targets. Each activity was tested with the design matrix of the pseudobulk tests (S6 Table b); BH FDR within the "
             "ten top CellOracle TFs (primary family) and within all TFs. r: Pearson correlation of activity with the donor Axis A score.", h, dd).save()
print("done")
