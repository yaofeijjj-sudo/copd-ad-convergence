# 方法学核查（2026-10-19）

依据：修改意见 5.1。**稿件未修改，候选结果。** 本分支只做核对和一项背景敏感性分析，没有新的疾病推断。

## 1. 模块重叠的 Fisher 检验背景（[fisher_128_pairs_background_comparison.tsv](results/fisher_128_pairs_background_comparison.tsv)）

**核对**：原 128 个模块对的检验中，每个模块集合都限定在同一个 3,623 个共享非灰基因的背景内。本分支逐一重算，128 对的交集数、模块大小和 P 值与原表一致，脚本中有断言。

**敏感性分析**：背景改为两侧共同通过过滤的 4,832 个基因（含灰模块），模块定义不变，同样做 128 对单侧 Fisher 检验，在 128 对内 BH 校正。

| 模块对 | 原背景 3,623：交集 / OR / FDR | 含灰背景 4,832：交集 / OR / FDR |
|---|---|---|
| Axis A（肺 red × 脑 blue） | 137 / 4.95 / 5.2 × 10⁻²⁹ | 137 / 4.29 / 1.7 × 10⁻²⁶ |
| Axis B（turquoise × turquoise） | 434 / 1.76 / 3.1 × 10⁻¹² | 434 / 1.72 / 1.6 × 10⁻¹² |

- 显著模块对：原背景 18 对，含灰背景 17 对，两者共有 16 对。
- 判定改变的 3 对都在 FDR 0.05 附近，且都不是 Axis A 或 Axis B：
  - brown × brown：0.33 → 0.035；
  - green × black：0.018 → 0.073；
  - yellow × magenta：0.033 → 0.052。
- 两种背景下，red × blue 都是 P 值最小的一对。

**结论**：只用非灰基因作背景，不影响 Axis A 与 Axis B 的选取，也不影响它们的显著性。

## 2. 单细胞供体流转（[single_cell_donor_flow.tsv](results/single_cell_donor_flow.tsv)）

**肺 GSE302339**：各步骤都是 56 名供体。
- 定位：211,464 个质控后细胞。
- 供体层面关联与伪 bulk：116,558 个髓系细胞，33 名 COPD、23 名对照，每名供体都 ≥ 20 个细胞。
- CellOracle：从 64,542 个肺泡巨噬细胞中固定抽取 30,000 个，56 名供体全部保留。

**脑 SEA-AD DFC**（2026-06-22 版，`Brain Region` = DFC）：

| 分析 | 细胞核 | 供体 | 规则 |
|---|---|---|---|
| 定位 | 1,357,145 | 80 | 发布的 DFC 全部细胞核 |
| CellOracle（Microglia 2） | 25,056 | 80 | DFC 全部 Micro-PVM_2，不分测序方法 |
| CPS 关联与伪 bulk（Microglia 2） | 23,198 | 79 | `method` = 10Xv3.1、`experiment_component_failed` = False、`Used in analysis` = True，且每名供体 ≥ 20 个细胞核 |

**80 → 79 的原因**（[run_02.log](qa/run_02.log)）：
- CPS 分析去掉了 10x Multiome 的 1,763 个细胞核，以及失败实验组分中的 95 个细胞核。
- 供体 H21.33.020 在 DFC 的 93 个 Micro-PVM_2 细胞核全部来自失败的实验组分，过滤后为 0，因此不进入 CPS 分析。
- 这组条件是 CPS 分析当时写定的主筛选，记录在定位分支的 filter audit（`细胞定位/analysis/06_single_cell_projection/results/SEAAD_DFC_axisA_localization_20260831_candidate/tables/SEAAD_DFC_Immune_precheck_filter_audit.tsv`）中。

**CellOracle 缺少的 1 个 Axis A 成员**：FER1L4（Entrez 80307）。肺、脑模型中都没有它唯一且稳定的 Ensembl 特征，标为缺失，没有补零（`CellOracle_AxisA_Entrez137_20261001/README.md`）。因此两侧可测成员都是 136 个。

**写入稿件时**：方法或补充表注明脑 CellOracle 与 CPS 分析使用的细胞核范围不同，并写明上面的筛选条件。

## 3. 脑发现队列的年龄极端值（已有结果，未重跑）

修改意见要求“不剔除年龄极端供体、合理调整年龄”的敏感性分析。这项分析已在 2026-07-24 完成并预先登记（`WGCNA_全部图件/复现资料/analysis/02_brain_wgcna/results/T04_AD_age_outlier_fixed_module_sensitivity_20260724/`）。
- 做法：恢复按组内 1.5 × IQR 规则剔除的 7 名供体，共 466 个样本（156 对照、310 AD）；固定原模块和特征向量载荷，以疾病 + 年龄 + 性别建模。
- 结果：blue β 由 0.0518 变为 0.0528，turquoise 为 −0.0500；两者在 16 个非灰模块中的排名不变。所有模块的疾病效应向量相关 r = 0.998。
- 结论：剔除年龄极端值不影响结果。稿件可以在方法中加一句并引用该结果。

## 运行

```
_workspace/envs/r/bin/Rscript 03_candidate_analysis/Methods_checks_20261019/code/01_fisher_background.R
_workspace/envs/celloracle/bin/python -I 03_candidate_analysis/Methods_checks_20261019/code/02_brain_donor_recount.py
```

均从项目根目录运行。
