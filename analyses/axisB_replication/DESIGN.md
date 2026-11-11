# Axis B 固定评分的独立复现（2026-10-19，运行前写定）

依据：修改意见 5.1。原文："Axis B 保留为次级发现，固定 434 个成员，补充独立肺、脑队列评分方向、覆盖率与区间。验证稳定时作为伴随代谢改变简述；不稳定时移入补充材料并如实说明。"作者 10-18 决定照此执行。

## 问题与估计量

- **问题**：固定的 434 个 Axis B 成员，在独立组织队列中是否随疾病降低？
- **估计量**：标准化 Axis B 评分（singscore simpleScore，在本队列分析样本内 z 化）的调整后病例–对照差异，单位为 SD。
- **独立单位**：供体。

## 输入（冻结）

- **成员**：`WGCNA_original_modules_Entrez_correction_20261001/results/cross/AxisB_confirmed_genes.csv`（434 个，两侧 turquoise 交集）。
- **队列与模型**：`AxisA_member_direction_robustness_20261018/source_data/cohorts.rds`，以及同一分支 `02_member_effects_and_sensitivity.R` 中的 `prep()`/`score_eff()`，原样复用。只把基因集换成 Axis B，与 Axis A 的复现完全同一流程。
- **样本排除**：GSE48350 的 GSM1176257（作者 10-19 决定，见 `adopted_sample_exclusions_20261019.tsv`）。

## 分析

1. **主分析（确认性）**：独立队列 GSE47460（肺），以及 6 个 AD 队列用 DerSimonian–Laird 合并。与 Axis A 的复现相同。
2. **描述**：发现队列 GSE76925、GSE33000（存在循环性，只作描述）；每个队列的可测成员数。
3. **敏感性**：
   - REML + Hartung–Knapp；
   - 只用所有验证平台共同可测的成员；
   - 不去掉 GSM1176257（便于与 10-02 的结果对照）。
4. **次要**：AD 队列中加入 CIBERSORTx（版本 B）估计的神经元比例（Excitatory + Inhibitory 之和，标准化），再合并。依据是 10-16 的 S7 Table：脑侧 Axis B 来源模块的关联约 94% 随神经元比例变化。

## 判定规则（写定，不按结果调整）

- **稳定**：GSE47460 与 AD 合并效应都为负，且 95% CI 都不含 0。按修改意见，作为伴随的代谢改变在正文中简述。
- **其他情况为不稳定**：移入补充材料，并如实写明哪一侧没有复现。
- 神经元比例调整只用于解读，不参与判定。

## 风险

- Axis B 富含氧化磷酸化/神经元相关基因，脑侧受神经元丢失影响强，可能主要反映细胞组成。
- 10-02 的结果显示，GSE47460 中 Axis B 效应约为 0；脑侧结果随读出方式改变。
