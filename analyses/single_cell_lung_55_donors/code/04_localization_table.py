# Two-organ localization table (REF/localization.tsv layout) with the 55-donor lung row; Holm across the two organs as in model_and_null.R.
# Run from project root: _workspace/envs/celltypist/bin/python -I 03_candidate_analysis/Lung_sc_donor_merge_20261029/code/04_localization_table.py
import pandas as pd
T = '03_candidate_analysis/Lung_sc_donor_merge_20261029/results/tables/'
ref = pd.read_csv('05_github_release/copd-ad-convergence/validation/reference/localization.tsv', sep='\t')
lung = pd.read_csv(T + 'lung_broad_primary_contrast.tsv', sep='\t')
out = pd.concat([lung, ref[ref.cohort == 'SEA-AD'].drop(columns='Holm_two_organs')], ignore_index=True)
p = out.p.values; o = p.argsort(); h = [0.0, 0.0]; run = 0.0
for k, i in enumerate(o): run = max(run, min(1.0, (2 - k) * p[i])); h[i] = run
out['Holm_two_organs'] = h
assert out.n_donors.tolist() == [55, 80]
out.to_csv(T + 'localization.tsv', sep='\t', index=False); print(out[['cohort', 'effect', 'ci_low', 'ci_high', 'p', 'n_donors', 'Holm_two_organs']])
