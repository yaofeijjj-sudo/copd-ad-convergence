"""Rebuild the two supplementary tables from the official-singscore results, keeping the row order,
labels and number formatting of the 2026-10-03 tables (AxisA_disease_specificity_redesign_20261003/results)."""
import csv, math
OLD = 'AxisA_disease_specificity_redesign_20261003/results/'
NEW = 'AxisA_singscore_official_rerun_20261004/results/'
def rd(p): return list(csv.DictReader(open(p), delimiter='\t'))
def r2(x): return '' if x in ('', 'NA', None) else f"{round(float(x), 2):g}" if abs(round(float(x),2))>=0.01 or float(x)==0 else f"{round(float(x),2):.2f}"
def s3(x):
    if x in ('', 'NA', None): return ''
    v = float(x); return f"{v:.3g}"
# disease scope
old = list(csv.reader(open(OLD + 'Supplementary_Table_AxisA_disease_scope.csv')))
new = rd(NEW + 'axisA_specificity_results.tsv'); assert len(old) - 1 == len(new) == 20
h = old[0]; col = {k: i for i, k in enumerate(h)}
out = [h]
for row, r in zip(old[1:], new):
    assert row[col['Cases']] == r['n_case'] and row[col['Controls']] == r['n_control']
    row = list(row)
    row[col['Effect (SD)']] = r2(r['effect_SD']); row[col['95% CI low']] = r2(r['ci_low']); row[col['95% CI high']] = r2(r['ci_high'])
    row[col['P']] = s3(r['p']); row[col['FDR (family)']] = s3(r['fdr_family'])
    row[col['cameraPR direction']] = r['camera_direction']; row[col['cameraPR P']] = s3(r['camera_p'])
    row[col['Matched-random empirical P']] = s3(r['random_empirical_p'])
    out.append(row)
csv.writer(open(NEW + 'Supplementary_Table_AxisA_disease_scope.csv', 'w', newline='')).writerows(out)
# AD reference meta
old = list(csv.reader(open(OLD + 'Supplementary_Table_AD_reference_meta.csv')))
new = rd(NEW + 'AD_reference_meta.tsv'); assert len(old) - 1 == len(new) == 7
h = old[0]; col = {k: i for i, k in enumerate(h)}
out = [h]
for row, r in zip(old[1:], new):
    row = list(row)
    for k, v in [('Effect (SD)', r2(r['effect_SD'])), ('SE', s3(r['se'])), ('95% CI low', r2(r['ci_low'])), ('95% CI high', r2(r['ci_high'])),
                 ('P', s3(r['p'])), ('cameraPR direction', r.get('camera_direction', '')), ('cameraPR P', s3(r.get('camera_p'))),
                 ('Matched-random empirical P', s3(r.get('random_empirical_p'))), ('I2 (%)', '' if r.get('I2') in ('', 'NA', None) else f"{float(r['I2']):.1f}"), ('tau2', s3(r.get('tau2'))), ('Heterogeneity P', s3(r.get('Q_p')))]:
        if k in col: row[col[k]] = v
    out.append(row)
csv.writer(open(NEW + 'Supplementary_Table_AD_reference_meta.csv', 'w', newline='')).writerows(out)
print('written')
