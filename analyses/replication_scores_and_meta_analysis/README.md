# Axis A scores recomputed with the official singscore package (2026-10-04)

**Why:** The author asked that the published algorithm be used directly rather than an in-house equivalent. The earlier branches computed "singscore-equivalent" mean within-sample ranks (ties = average, scaled by n + 1). This branch replaces only the scoring step with the official Bioconductor implementation, `singscore::rankGenes()` + `singscore::simpleScore()` (singscore 1.30.0; Foroutan et al. 2018).

**Unchanged:** Axis A (137 Entrez IDs), datasets, sample rules, covariates, linear models, random seeds, the DerSimonian–Laird pooling and leave-one-out, and BH families. Inputs are read from `AxisA_disease_specificity_redesign_20261003` (read-only).

**Expected difference:** After per-dataset standardisation, the two scores differ only in tie handling. Official `rankGenes` uses ties = "min".

**Run (project root):**
```
R_LIBS=_workspace/r_libs LANG=en_US.UTF-8 _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/02_run_specificity.R
R_LIBS=_workspace/r_libs_affy:_workspace/r_libs LANG=en_US.UTF-8 _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_singscore_official_rerun_20261004/code/04_ad_reference_meta.R
```

## Results (vs the 2026-10-03 in-house score)

- **Microarray datasets:** identical results (GSE47460 both platforms, GSE84422, GSE118553, GSE48350, GSE36980, GSE53987, GSE10072, GSE33000).
- **RNA-seq datasets:** only rounding-level changes; no direction or call changed.
  - GSE150910 CHP: 0.64 → 0.63, FDR 0.0015 → 0.0017
  - GSE140841: 0.79 → 0.76
  - GSE313307: P 0.024 → 0.025
- **AD pooled:** 0.46 (0.22–0.70), P = 1.7e-4, I² = 21% → 0.46 (0.22–0.69), P = 1.5e-4, I² = 20%.
- **Leave-one-out:** unchanged at 0.37–0.53, all P < 0.003.

## Files

- **Figures:** `figures/` holds Figure 4d (`Figure4_AxisA_disease_scope_forest.pdf`) and Supplementary Fig. S1 (`SuppFigure_AD_reference_meta.pdf`), redrawn with the original scripts (paths only changed). Text collisions: 0.
- **Supplementary tables:** `results/Supplementary_Table_*.csv`, rebuilt by `code/06_supplementary_tables.py` with the same rows and format.
- **Status:** adopted by the author on 2026-10-04. The manuscript now cites these numbers.

## 2026-10-27 补记：GPL6480 次要复现写入正文

- 作者 10-27 决定把预先设定的次要复现 GSE47460-GPL6480（COPD 75 vs 对照 17，Axis A 134/137 可测）写入正文：1.17 SD（95% CI 0.58–1.76，P = 1.8 × 10⁻⁴），数值取自 `results/axisA_specificity_results.tsv`，未重新计算。
- 供体核查：GPL6480 的 153 名供体与 GPL14550 的 429 名供体无重叠（按 LT 编号）。与肺发现队列 GSE76925 比较，5 对供体年龄、性别和 FEV₁ 占预计值百分比（支气管扩张剂前或后）相同，其表达相关在 92 名 GPL6480 供体中排第 19–91 位，均非同一人（[qa/GPL6480_vs_GSE76925_fingerprint_pairs_20261027.tsv](qa/GPL6480_vs_GSE76925_fingerprint_pairs_20261027.tsv)；方法同 `AxisA_member_direction_robustness_20261018/code/05_followup_checks.R`）。
- 写入位置：结果第 2 部分、方法、Table 1、S4 Table a 与 d。
