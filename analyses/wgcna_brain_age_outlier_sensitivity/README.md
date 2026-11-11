# T04 AD age-outlier fixed-module sensitivity

T04 restored the seven GSE33000 samples excluded by the within-group 1.5×IQR age rule without rebuilding the AD network. The current probes, 8,807 genes, module assignments, preprocessing reference and module-eigengene loadings were frozen; `GSM1424097` remained excluded as the independent connectivity outlier.

## Decision

**GO.** The sensitivity cohort contained 466 samples (156 controls and 310 AD cases). AD blue remained positive (`β=0.05280`) and AD turquoise remained negative (`β=−0.05003`); both retained their current ranks of 2 and 3 among the 16 non-grey modules. The non-grey disease-effect vectors were highly concordant (`r=0.99758`). No preregistered upgrade trigger fired, so a full WGCNA sensitivity was not launched and T05 is unlocked.

## Start here

- Full Chinese report: `docs/T04_age_outlier_sensitivity_report_zh_CN_20260724.md`
- Preregistered contract: `docs/T04_preregistered_analysis_contract_20260724.md`
- Comparison table: `tables/t04_fixed_module_disease_effect_comparison_20260724.csv`
- Acceptance checks: `tables/t04_acceptance_checks_20260724.csv`
- Go/no-go record: `tables/t04_go_no_go_decision_20260724.csv`
- Figure and source data: `figures/` and `tables/t04_figure_source_data_20260724.csv`

