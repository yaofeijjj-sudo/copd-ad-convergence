"""Build S1 Data (values underlying Fig 2-8 and Fig A-B in S1 File) from the frozen sources used by code/main_figures/fig2-8.R
and code/supplement/build_S1_S2_Fig.R (Fig A-B sheets added 2026-10-30).

Panels already exported by the figure scripts are read from code/main_figures/source_data/; the remaining panels
are read from the same files the figure scripts read. Nothing is recomputed. Fig 1 is a schematic. Per-cell UMAP values (Fig 4a-d) are
summarised by cell group because the brain table (1,357,145 nuclei) exceeds the spreadsheet row limit.
Run from the project root:
  _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/build_S1_Data.py
"""
import gzip, csv, io
from pathlib import Path
import openpyxl
from openpyxl.styles import Font

ROOT = Path.cwd()
REF = ROOT / '05_github_release/copd-ad-convergence/validation/reference'
CAND = ROOT / '03_candidate_analysis'
SRC = ROOT / '06_resubmission_20261008/code/main_figures/source_data'
SC = ROOT / 'cell_localization/SYN01_Entrez_downstream_20261001/source_data'
EN = ROOT / 'functional_enrichment/SYN01_Entrez_downstream_20261001/source_data'
LM = CAND / 'Lung_sc_donor_merge_20261029/results/tables'  # 10-29: GSE302339 cood8 = copd8, lung donor-level values for 55 donors
CO = LM  # CellOracle ranking with 55-donor lung means (CellOracle_ranking_scale_20261019 layout)
OUT = ROOT / '06_resubmission_20261008/supplementary/S1_Data.xlsx'


def read(path, delim=None):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f:
        text = f.read()
    if delim is None:
        delim = ',' if str(path).endswith('.csv') else '\t'
    rows = list(csv.reader(io.StringIO(text), delimiter=delim))
    return rows[0], rows[1:]


def num(x):
    try:
        v = float(x)
        return int(v) if v.is_integer() and 'e' not in x.lower() and '.' not in x else v
    except (ValueError, AttributeError):
        return x


def umap_summary(path):
    head, rows = read(path)
    g, s = head.index('group'), head.index('AxisA_Entrez_UCell')
    acc = {}
    for r in rows:
        v = float(r[s]); a = acc.setdefault(r[g], [])
        a.append(v)
    out = []
    for grp, vals in sorted(acc.items()):
        vals.sort(); n = len(vals)
        med = vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
        out.append([grp, n, sum(vals) / n, med])
    return ['group', 'cells_or_nuclei', 'mean_AxisA_UCell', 'median_AxisA_UCell'], out, len(rows)


def top10(path):
    head, rows = read(path)
    keep = [head.index(c) for c in ('TF', 'rank_shared_down', 'shared_down', 'shared_up', 'opposite')]
    rk = head.index('rank_shared_down')
    sel = sorted([r for r in rows if float(r[rk]) <= 10], key=lambda r: float(r[rk]))
    assert len(sel) == 10
    return [head[k] for k in keep], [[r[k] for k in keep] for r in sel]


sheets = []  # (sheet name, panel, description, header, rows, source)
EN_DIR = {'cell_localization': 'cell_localization', 'functional_enrichment': 'functional_enrichment',
          'WGCNA_figures': 'WGCNA_figures', 'reproduction_materials': 'reproduction_materials'}  # English labels for the submitted file


def add(name, panel, desc, src, header=None, rows=None):
    if header is None:
        header, rows = read(src)
    rel = Path(src).relative_to(ROOT).parts
    sheets.append((name, panel, desc, header, rows, '/'.join(EN_DIR.get(p, p) for p in rel)))


add('Fig2a', 'Fig 2a', 'Lung module-trait Pearson correlations and P values (GSE76925, n = 147)', REF / 'lung_module_traits.csv')
add('Fig2b', 'Fig 2b', 'Brain module-trait Pearson correlations and P values (GSE33000, n = 459)', REF / 'brain_module_traits.csv')
add('Fig2c', 'Fig 2c', 'One-sided Fisher tests for all 128 lung-brain module pairs (BH FDR)', SRC / 'Fig2c_module_overlap.tsv')
h, r = read(SRC / 'Fig2c_module_overlap.tsv')
r = [x for x in r if (x[0], x[1]) in (('red', 'blue'), ('turquoise', 'turquoise'))]; assert len(r) == 2
add('Fig2d', 'Fig 2d', 'Shared-gene counts of the two pairs defining Axis A (red-blue, 137) and Axis B (turquoise-turquoise, 434)', SRC / 'Fig2c_module_overlap.tsv', h, r)
add('Fig2e', 'Fig 2e', 'Covariate-adjusted standardized case-control effects of the four source modules', REF / 'adjusted_effects.csv')
add('Fig3a', 'Fig 3a', 'Module preservation Zsummary in GSE47460 and GSE84422', SRC / 'Fig3a_preservation.tsv')
add('Fig3b', 'Fig 3b', 'sigQC median pairwise correlation of Axis A members and expression-matched random gene sets (Table C a)', SRC / 'Fig3b_sigqc.tsv')
add('Fig3c', 'Fig 3c', 'Axis A disease effects (SD, 95% CI) by dataset and random-effects pooling', SRC / 'Fig3c_forest.tsv')
add('Fig3d', 'Fig 3d', 'Per-member standardized disease effects in GSE47460 and the pooled AD datasets', SRC / 'Fig3d_member_effects.tsv')
for organ, f in (('lung', SC / 'Lung_UMAP_cells.tsv.gz'), ('brain', SC / 'Brain_UMAP_cells.tsv.gz')):
    h, r, n = umap_summary(f)
    panel = 'Fig 4a, b' if organ == 'lung' else 'Fig 4c, d'
    add(f'Fig4{"ab" if organ == "lung" else "cd"}', panel,
        f'{organ.capitalize()} UMAP cell groups and Axis A UCell scores summarized by group ({n:,} cells or nuclei; per-cell values exceed the spreadsheet row limit and are regenerated by the code repository)',
        f, h, r)
add('Fig4e', 'Fig 4e', 'Donor-mean Axis A UCell by lung myeloid subtype', LM / 'Lung_plotted_points.tsv')
add('Fig4f', 'Fig 4f', 'Donor-mean Axis A UCell by brain immune Supertype', SC / 'Brain_plotted_points.tsv')
add('Fig4g', 'Fig 4g', 'Lung within-donor UCell differences, myeloid minus other cell classes, for Axis A and 5,000 expression-matched random gene sets', LM / 'lung_matched_null_draws.tsv')
add('Fig4h', 'Fig 4h', 'Brain within-donor UCell differences, immune minus other cell classes, for Axis A and 5,000 expression-matched random gene sets', REF / 'brain_null.tsv')
add('Fig4gh_obs', 'Fig 4g, h', 'Observed Axis A within-donor UCell differences, target minus other cell classes (donor-clustered CR2 models)', LM / 'localization.tsv')
add('Fig5a', 'Fig 5a', 'Donor median Axis A UCell in lung myeloid cells (55 donors)', SRC / 'Fig5a_lung_donors.tsv')
add('Fig5b', 'Fig 5b', 'Donor Axis A UCell in brain Microglia 2 and local CPS (79 donors)', SRC / 'Fig5b_brain_donors.tsv')
add('Fig5c', 'Fig 5c', 'Donor pseudobulk fry and CAMERA tests of Axis A', SRC / 'Fig5c_pseudobulk.tsv')
add('Fig5d', 'Fig 5d', 'Axis A disease effects before and after adding the CIBERSORTx target-cell proportion', SRC / 'Fig5d_composition.tsv')
add('Fig5e', 'Fig 5e', 'Per-sample standardized Axis A score and estimated target-cell proportion in eight tissue cohorts', SRC / 'Fig5e_samples.tsv')
add('Fig5f', 'Fig 5f', 'Case-control difference in target-cell proportion (SD, 95% CI), P and r with the Axis A score', SRC / 'Fig5f_proportion_difference.tsv')
add('Fig6a_terms', 'Fig 6a', 'Enriched terms shown (gene ratio, BH q, gene count)', EN / 'axisA_functional_enrichment_complete_family_terms.csv')
add('Fig6a_genes', 'Fig 6a', '18 most-annotated Axis A genes shown', EN / 'axisA_functional_enrichment_selected_genes.csv')
add('Fig6a_links', 'Fig 6a', 'Gene-term membership links', EN / 'axisA_functional_enrichment_complete_family_links.csv')
add('Fig6b', 'Fig 6b', 'Signature effects and Axis A effects after adding each signature (SD, 95% CI)', SRC / 'Fig6b_signatures.tsv')
add('Fig6b_r', 'Fig 6b', 'Correlation of each signature score with the Axis A score (GSE47460; median across six AD datasets)', SRC / 'Fig6b_signature_r.tsv')
add('Fig6c', 'Fig 6c', 'MAGMA coefficients (beta, SE), P and Holm P', SRC / 'Fig6c_magma.tsv')
add('Fig7a', 'Fig 7a', 'Peripheral blood Axis A score differences (SD, 95% CI) and CAMERA P', SRC / 'Fig7a_blood.tsv')
add('Fig7b', 'Fig 7b', 'Axis A disease effects in the target and comparator disease contexts', SRC / 'Fig7b_scope.tsv')
add('Fig8b', 'Fig 8b', 'Shared decreases, within-model response percentile and out-degree for 357 TFs', SRC / 'Fig8b_TF_breadth.tsv')
h, r = top10(CO / 'all357_TF_ranking.tsv')
add('Fig8c', 'Fig 8c', 'Shared-decrease, shared-increase and opposite-direction counts for the top ten TFs', CO / 'all357_TF_ranking.tsv', h, r)
# 11-09: Fig 8d, e show STAT3 and CEBPB only (the LYL1 panel was dropped); former Fig 8g is now Fig 8f
h, r = read(CO / 'Fig6def_source_signed_logratio.tsv')
r = [x for x in r if x[h.index('TF')] in ('STAT3', 'CEBPB')]; assert len(r) == 272
add('Fig8de', 'Fig 8d, e', 'Predicted delta_X (natural-log scale) of 136 Axis A genes after STAT3 and CEBPB knockout', CO / 'Fig6def_source_signed_logratio.tsv', h, r)
add('Fig8f', 'Fig 8f', 'TF activity correlation with the donor Axis A score and disease-association t', SRC / 'Fig8f_TF_activity.tsv')
# Fig A-B in S1 File: the files read by build_S1_S2_Fig.R
ST = ROOT / 'WGCNA_figures/reproduction_materials/figure_revision_20260915_Arial8_PDF/source_data'
AB = CAND / 'AxisB_fixed_score_replication_20261019/results'
add('FigA_lung', 'Fig A a (S1 File)', 'Lung network scale-free fit (signed R2) and mean connectivity by soft-threshold power (GSE76925; selected power 8)', ST / 'lung_soft_threshold.csv')
add('FigA_brain', 'Fig A b (S1 File)', 'Brain network scale-free fit (signed R2) and mean connectivity by soft-threshold power (GSE33000; selected power 12)', ST / 'brain_soft_threshold.csv')
add('FigB_a', 'Fig B a (S1 File)', 'Axis B enrichment terms shown (gene ratio, BH q, gene count)', EN / 'axisB_functional_enrichment_complete_family_terms.csv')
h, r = read(AB / 'axisB_by_cohort.tsv')
r = [x for x in r if x[h.index('role')] == 'independent']; assert len(r) == 7
add('FigB_b_cohorts', 'Fig B b (S1 File)', 'Axis B disease effects (SD, 95% CI) in GSE47460 and the six AD datasets', AB / 'axisB_by_cohort.tsv', h, r)
h, r = read(AB / 'axisB_pooled.tsv')
r = [x for x in r if x[0] in ('AD pooled, DL (primary)', 'AD pooled, + neuronal proportion (secondary)')]; assert len(r) == 2
add('FigB_b_pooled', 'Fig B b (S1 File)', 'Pooled AD Axis B effects without and with the CIBERSORTx neuronal proportion', AB / 'axisB_pooled.tsv', h, r)

wb = openpyxl.Workbook()
ws = wb.active; ws.title = 'README'
ws.append(['S1 Data. Values underlying Figs 2-8 and Fig A-B in S1 File. Fig 1 is a schematic; its thumbnails reuse values of Figs 3-5 and 7.']); ws['A1'].font = Font(bold=True)
# 11-09: the README no longer lists local project paths, which readers cannot open. Those paths are kept for the
# authors in code/supplement/parts/S1_Data_source_map.tsv.
ws.append(['Each sheet holds the values plotted in the named panel. The analyses that produced them are described in the Supplementary Methods in S1 File, the input datasets are listed in Table 1 of the main text and Table D in S2 File, and the analysis code is available at https://github.com/yaofeijjj-sudo/copd-ad-convergence. Fig 4a-d per-cell values are summarized by cell group.'])
ws.append([])
ws.append(['Sheet', 'Panel', 'Content', 'Rows'])
for c in ws[4]: c.font = Font(bold=True)
with open(ROOT / '06_resubmission_20261008/code/supplement/parts/S1_Data_source_map.tsv', 'w', encoding='utf-8') as fh:
    fh.write('sheet\tpanel\tsource_file_project_relative\n')
    for name, panel, desc, header, rows, src in sheets: fh.write(f'{name}\t{panel}\t{src}\n')
for name, panel, desc, header, rows, src in sheets:
    ws.append([name, panel, desc, len(rows)])
    sh = wb.create_sheet(name)
    sh.append(header)
    for c in sh[1]: c.font = Font(bold=True)
    for row in rows:
        sh.append([num(x) if isinstance(x, str) else x for x in row])
wb.save(OUT)
print(OUT, 'sheets:', len(sheets))
