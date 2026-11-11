# Axis A/B MAGMA 基因窗口敏感性分析（2026-10-27）

- **目的：** 补报 2026-07-15 预分析方案（`GWAS和外周血/analysis/05_disease_specificity_gwas/results/current_axis_full135_full418_multiancestry_gwas_magma_20260715/preanalysis_protocol_20260715.md` 第 71–72、85 行）预先设定的两种注释窗口敏感性分析：±10 kb；上游 35 kb、下游 10 kb。各窗口内对两条轴作 Holm 校正。原方案在旧版 135/418 基因集上运行过，本次只对现行固定基因集（Axis A 137 / Axis B 434，Entrez）重跑基因集步骤。
- **冻结输入：** 同一 AD GWAS（GCST90704648，不含大型生物样本库，13,648,198 个 SNP）的各窗口 `genes.raw`；基因集 `GWAS和外周血/SYN01_Entrez_downstream_20261001/results/magma_inputs/locked_axes.sets`。MAGMA v1.10。COPD 不在该方案的窗口敏感性范围内，且其窗口输入与主分析的去重输入不同，未运行。
- **代码：** `code/run_windows.sh`；日志 `logs/run.log`；结果 `results/AD_*.gsa.out`。

| 窗口 | Axis A（130 个基因）P | Axis B（409 个基因）P | Holm（两轴） |
|---|---|---|---|
| 仅基因内（主分析，复现） | 7.1 × 10⁻⁶ | 0.135 | A 1.4 × 10⁻⁵；B 0.135 |
| ±10 kb | 8.5 × 10⁻⁴ | 0.116 | A 0.0017；B 0.116 |
| 上游 35 kb / 下游 10 kb | 0.051 | 0.0086 | A 0.051；B 0.017 |

- **解读：** Axis A 的 AD 富集在仅基因内和 ±10 kb 窗口下显著，窗口扩到上游 35 kb 时减弱至 P = 0.051；此时 Axis B 达到 Holm P = 0.017。已写入正文结果、方法和 S10 Table（f）。
- **输入与主分析的差别（10-28 核实）：** 三个窗口的 `genes.raw` 都来自 07-15 方案，每个 SNP 的样本量为有效样本量（如 rs10 为 120,805）；主分析（S10 Table a–e、图 4c）改用病例数加对照数（rs10 为 145,092），并排除重复 rsID，SNP 与 P 值相同（`GWAS和外周血/SYN01_Entrez_downstream_20261001/qa/AD_N_input_identity.json`）。因此表中“仅基因内”一行是窗口比较的同口径参照，与主分析略有差别（Axis A P = 7.12 × 10⁻⁶ 对 7.18 × 10⁻⁶）；S10 Table f 中该行标为 “Genic only (window-analysis input)”，表注写明原因。
