"""10-29: GSE302339 donors cood8 and copd8 are one person; lung donor-level values rerun with 55 donors in
03_candidate_analysis/Lung_sc_donor_merge_20261029. Updates the lung rows of S6 Table (a, b, c, e, g), S12 Table (lung columns)
and S13 Table (lung rows) in parts/ in place; brain rows are untouched. Runs after append_S6g_alveolar_macrophage.py and
update_S6_region_label.py, before build_S1_S2_File.py.
Run from project root: _workspace/envs/celltypist/bin/python -I 06_resubmission_20261008/code/supplement/update_lung_donor_merge_20261029.py
"""
import csv, gzip
from pathlib import Path
import openpyxl

ROOT = Path.cwd(); PARTS = ROOT / '06_resubmission_20261008/code/supplement/parts'
T = ROOT / '03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables'


def tsv(p):
    op = gzip.open if str(p).endswith('.gz') else open
    with op(p, 'rt') as f: return list(csv.DictReader(f, delimiter='\t'))


def header(ws, row=4): return {c.value: c.column for c in ws[row] if c.value is not None}


def setv(ws, r, h, name, v): ws.cell(r, h[name]).value = v


num = float

# ---------------- S6 ----------------
wb = openpyxl.load_workbook(PARTS / 'S6_Table.xlsx')
# (a) lung broad cell types
ws = wb['S6a']; h = header(ws)
tests = {r['broad_celltype']: r for r in tsv(T / 'axisA_donor_COPD_vs_control.tsv')}
dm = [r for r in tsv(T / 'axisA_donor_broad_medians.tsv') if r['eligible'] == 'TRUE']
done = 0
for r in range(5, ws.max_row + 1):
    if not str(ws.cell(r, 1).value).startswith('Lung'): continue
    g = ws.cell(r, h['Cell group (original label)']).value; t = tests[g]; d = [x for x in dm if x['broad_celltype'] == g]
    n = len(d); assert ws.cell(r, h['Median donor UCell, all donors']).value is None   # lung rows report the two group medians only
    assert int(t['n_group1']) + int(t['n_group2']) == n
    for k, v in [('Donors (n)', n), ('Control donors (n)', int(t['n_group1'])), ('COPD donors (n)', int(t['n_group2'])),
                 ('Cells or nuclei (n)', sum(int(x['n_cells']) for x in d)), ('Median donor UCell, control', num(t['median_group1'])),
                 ('Median donor UCell, COPD', num(t['median_group2'])),
                 ('Effect', num(t['delta_group2_minus_group1'])), ('P', num(t['wilcoxon_p'])), ('BH FDR', num(t['wilcoxon_fdr_bh']))]:
        setv(ws, r, h, k, v)
    done += 1
assert done == 6
# (b) lung pseudobulk rows
ws = wb['S6b']; h = header(ws)
pb = {r['model']: r for r in tsv(T / 'pseudobulk_fry_camera.tsv') if r['organ'] == 'lung'}
done = 0
for r in range(5, ws.max_row + 1):
    if ws.cell(r, 1).value != 'lung': continue
    x = pb[ws.cell(r, h['Model']).value]
    for k, c in [('Donors (n)', 'n_donors'), ('Axis A members after filtering (n)', 'axis_after_filter'), ('fry direction', 'fry_direction'),
                 ('fry P', 'fry_p'), ('CAMERA direction', 'camera_direction'), ('CAMERA P', 'camera_p'), ('CAMERA inter-gene correlation', 'camera_corr'),
                 ('fry BH FDR (organ, primary)', 'fry_BH_within_organ'), ('CAMERA BH FDR (organ, primary)', 'camera_BH_within_organ')]:
        v = x[c]; setv(ws, r, h, k, None if v in ('', 'NA') else (v if c.endswith('direction') else (int(v) if c in ('n_donors', 'axis_after_filter') else num(v))))
    done += 1
assert done == 5
# (c) donor counts
ws = wb['S6c']; h = header(ws)
for r in range(5, ws.max_row + 1):
    if ws.cell(r, 1).value != 'lung': continue
    assert ws.cell(r, h['Donors (n)']).value == 56; setv(ws, r, h, 'Donors (n)', 55)
    rule = ws.cell(r, h['Inclusion rule']).value
    rule = rule.replace('all 56 eligible: 33 COPD / 23 control', 'all 55 eligible: 32 COPD / 23 control').replace('all 56 donors retained', 'all 55 donors retained')
    setv(ws, r, h, 'Inclusion rule', rule)
ws.cell(2, 1).value = ('Numbers of cells or nuclei and donors after the inclusion rule of each analysis. GSE302339 subject IDs "cood8" (GSM9102255) '
                       'and "copd8" (GSM9102261) are two lobes of one donor (identical clinical records) and are counted as one donor.')
# (e) flagged libraries
ws = wb['S6e']; h = header(ws)
qc = tsv(T / 'myeloid_donor_test_library_check.tsv')
for r, x in zip((5, 6), qc):
    assert ws.cell(r, 1).value == x['analysis']
    for k, c in [('Donors (n)', 'donors'), ('Control donors (n)', 'control'), ('COPD donors (n)', 'copd'), ('Myeloid cells (n)', 'cells'),
                 ('Median donor UCell, control', 'median_control'), ('Median donor UCell, COPD', 'median_copd'),
                 ('Difference in medians (COPD − control)', 'delta'), ('P', 'P')]:
        setv(ws, r, h, k, int(x[c]) if c in ('donors', 'control', 'copd', 'cells') else num(x[c]))
# (g) alveolar macrophages
ws = wb['S6g']; h = header(ws)
pbr = {r['analysis']: r for r in tsv(T / 'AM_pseudobulk_fry_camera.tsv')}; uc = {r['analysis']: r for r in tsv(T / 'AM_donor_ucell_wilcoxon.tsv')}
rowmap = {'Alveolar macrophages (primary)': 'primary: AM', 'Alveolar macrophages, ≥50 cells per donor': 'sensitivity a: AM, >=50 cells',
          'Alveolar macrophages, 2 flagged libraries removed': 'sensitivity b: AM, flagged libraries removed',
          'Alveolar macrophage family (4 labels)': 'sensitivity c: AM family', 'All myeloid cells (reference; S6 Table a, b)': 'reference: Myeloid (S6 Table)'}
for r in range(5, ws.max_row + 1):
    a = rowmap[ws.cell(r, 1).value]; p, u = pbr[a], uc[a]; assert p['n_donors'] == u['n_donors']
    for k, v in [('Donors (n)', int(p['n_donors'])), ('COPD donors (n)', int(p['n_COPD'])), ('Control donors (n)', int(p['n_control'])),
                 ('Cells (n)', int(p['n_cells'])), ('Axis A members after filtering (n)', int(p['axis_after_filter'])),
                 ('fry direction', p['fry_direction']), ('fry P', num(p['fry_p'])), ('CAMERA direction', p['camera_direction']),
                 ('CAMERA P', num(p['camera_p'])), ('CAMERA inter-gene correlation', num(p['camera_corr'])),
                 ('Median donor UCell, control', num(u['median_control'])), ('Median donor UCell, COPD', num(u['median_COPD'])),
                 ('Difference in medians (COPD − control)', num(u['diff_medians'])), ('Hodges–Lehmann shift (COPD − control)', num(u['HL_shift'])),
                 ('95% CI low', num(u['HL_low'])), ('95% CI high', num(u['HL_high'])), ('Wilcoxon P', num(u['wilcoxon_p']))]:
        setv(ws, r, h, k, v)
sh = {r['disease_group']: r for r in tsv(T / 'AM_share_of_myeloid_summary.tsv')}
f2 = lambda x: f'{num(x):.2f}'
old = 'COPD 0.48 (0.38–0.60), control 0.48 (0.29–0.66).'
leg = ws.cell(2, 1).value; assert old in leg
ws.cell(2, 1).value = leg.replace(old, f"COPD {f2(sh['COPD']['median_AM_share'])} ({f2(sh['COPD']['q25'])}–{f2(sh['COPD']['q75'])}), "
                                       f"control {f2(sh['Control']['median_AM_share'])} ({f2(sh['Control']['q25'])}–{f2(sh['Control']['q75'])}).")
wb.save(PARTS / 'S6_Table.xlsx')

# ---------------- S12 (lung columns of the CellOracle ranking) ----------------
rk = {r['TF']: r for r in tsv(T / 'all357_TF_ranking.tsv')}
order = sorted(rk.values(), key=lambda r: -num(r['mean_abs_delta_lung']))   # same rank rule as build_S1_S13_tables.py
lrank = {r['TF']: i for i, r in enumerate(order, 1)}
wb = openpyxl.load_workbook(PARTS / 'S12_Table.xlsx'); ws = wb['S12a']; h = header(ws); n = 0
for r in range(5, ws.max_row + 1):
    tf = ws.cell(r, h['TF']).value
    if tf is None: continue
    x = rk[tf]; n += 1
    for k in ('rank_shared_down', 'shared_down', 'shared_up', 'opposite', 'lung_only_change', 'brain_only_change'):
        col = {'rank_shared_down': 'Rank (shared decreases)', 'shared_down': 'Decreased in both (n)', 'shared_up': 'Increased in both (n)', 'opposite': 'Opposite direction (n)',
               'lung_only_change': 'Changed in lung only (n)', 'brain_only_change': 'Changed in brain only (n)'}[k]
        assert int(ws.cell(r, h[col]).value) == int(x[k]), (tf, k)   # sign counts are unchanged by the merge
    for col, k in [('Mean signed delta_X, lung', 'mean_signed_delta_lung'), ('Mean |delta_X|, lung', 'mean_abs_delta_lung'),
                   ('Response percentile, lung', 'response_percentile_lung'), ('Shared decreases per 100 out-degree', 'shared_down_per_100_outdegree_mean')]:
        setv(ws, r, h, col, num(x[k]))
    setv(ws, r, h, 'Within-organ response rank, lung', lrank[tf])
assert n == 357
ws.cell(2, 1).value = ws.cell(2, 1).value.rstrip() + (' Lung delta_X is averaged within donors and then equally across the 55 donors of the subsample '
                                                       '(GSE302339 "cood8" and "copd8" counted as one donor).')
wb.save(PARTS / 'S12_Table.xlsx')

# ---------------- S13 (lung rows) ----------------
tf = {(r['model'], r['network'], r['TF']): r for r in tsv(T / 'tf_activity_all.tsv.gz') if r['organ'] == 'lung'}
wb = openpyxl.load_workbook(PARTS / 'S13_Table.xlsx'); ws = wb['S13a']; h = header(ws); n = 0
cols = [('Targets measured (n)', 'n_targets_measured', int), ('Activity logFC (exposure)', 'logFC', num), ('Mean activity', 'AveExpr', num), ('t', 't', num),
        ('P', 'P.Value', num), ('BH FDR, top-10 CellOracle TFs', 'FDR_top10', num), ('BH FDR, all TFs', 'FDR_all', num), ('Correlation with donor Axis A score (r)', 'r_with_axisA_score', num)]
for r in range(5, ws.max_row + 1):
    if ws.cell(r, 1).value != 'lung': continue
    x = tf[(ws.cell(r, h['Model']).value, ws.cell(r, h['Regulatory prior']).value, ws.cell(r, h['TF']).value)]; n += 1
    for col, k, fn in cols: setv(ws, r, h, col, None if x[k] in ('', 'NA') else fn(x[k]))
assert n == len(tf) == 6630, (n, len(tf))
wb.save(PARTS / 'S13_Table.xlsx')
print('S6 a/b/c/e/g, S12 lung columns and S13 lung rows updated')
