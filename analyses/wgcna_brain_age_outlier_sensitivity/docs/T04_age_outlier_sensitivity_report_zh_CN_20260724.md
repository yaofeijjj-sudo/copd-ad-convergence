# T04：GSE33000 年龄异常样本固定模块敏感性报告

日期：2026-07-24

## 结论

**T04 通过，判定为 GO。** 在保持 current 8,807 个基因、模块归属、参考表达尺度、模块 PC1 载荷和年龄/性别调整模型不变的条件下，将原先按组内 1.5×IQR 规则排除的 7 个年龄样本重新纳入后，AD blue 和 AD turquoise 的疾病效应方向、效应量、置信区间、FDR 和模块排序均保持。预先锁定的 Axis 主判据和全模块辅助判据全部未触发，因此不启动完整 WGCNA sensitivity，current 459 例网络继续作为主网络，T05 解锁。

这项结果支持的精确结论是：**当前 Axis A/Axis B 脑侧疾病方向不是由这 7 个年龄样本的排除所驱动。** 它不证明组内 1.5×IQR 是唯一或最优的年龄处理方法，也不消除 PMI、pH、RIN、细胞组成或其他未测混杂。

## 为什么使用固定模块投影

若直接用 466 例重新建网，模块成员、切树、合并和 eigengene 载荷都会同时变化，无法判断差异来自“恢复 7 个样本”还是“重估网络”。T04 因此只检验 current 结论对年龄样本选择的敏感性：

1. HD 样本继续排除；
2. 连接度异常样本 `GSM1424097` 继续排除；
3. 参考标准化目标由原年龄筛选后、连接度筛选前的 460 例学习；
4. current selected probes、8,807 个基因及其 17 个颜色标签完全冻结；
5. 每个模块的基因中心、标准差、PC1 载荷和方向只由 current 459 例学习；
6. 7 个恢复样本作为样本外观测投影到这些固定模块；
7. 对 current 459 例和敏感性 466 例分别拟合相同的 `ME ~ Disease_State + age_z + sex_male` 模型。

该设计是“固定网络内的样本纳入敏感性”，不是第三条发现网络。

## 样本流

| 阶段 | 对照 | AD | 合计 |
| --- | ---: | ---: | ---: |
| 排除 HD 后、年龄筛选前 | 157 | 310 | 467 |
| 组内 1.5×IQR 年龄筛选后 | 154 | 306 | 460 |
| current 网络（再排除 `GSM1424097`） | 153 | 306 | 459 |
| T04 固定模块敏感性（恢复 7 例，仍排除 `GSM1424097`） | 156 | 310 | 466 |

恢复样本为 3 名对照和 4 名 AD：

- 对照：`GSM1424124`（22 岁）、`GSM1424092`（95 岁）、`GSM1424154`（106 岁）；
- AD：`GSM1424080`（53 岁）、`GSM1423917`（54 岁）、`GSM1423926`（54 岁）、`GSM1423809`（55 岁）。

## 技术复现验收

全部 18 项验收通过：

- current 表达矩阵重建最大绝对误差：`0`；
- 17 个保存 eigengene 与固定 PC1 投影的最小 Pearson 相关：`1.000000`；
- eigengene 最大绝对误差：`1.61×10⁻¹⁵`；
- current 疾病系数最大复现误差：`2.36×10⁻¹⁶`；
- current 标准误最大复现误差：`1.82×10⁻¹⁷`；
- 样本数、基因数、模块数、恢复样本及连接度异常样本均与锁定合同一致。

因此，敏感性结果不是由不同探针、不同基因集、不同模块标签或不同 eigengene 定义造成。

## Axis A/B 主要结果

| 脑侧模块 | current n=459，β（95% CI） | current FDR | 恢复 7 例 n=466，β（95% CI） | 敏感性 FDR | 绝对效应变化 | 非灰模块绝对效应排序 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Axis A · AD blue | 0.05177（0.04333–0.06021） | 5.04×10⁻²⁸ | 0.05280（0.04485–0.06075） | 4.02×10⁻³² | 增强 1.98% | 2 → 2 |
| Axis B · AD turquoise | −0.05005（−0.05911 至 −0.04098） | 8.49×10⁻²⁴ | −0.05003（−0.05853 至 −0.04154） | 1.21×10⁻²⁶ | 衰减 0.029% | 3 → 3 |

blue 的正向疾病效应略增强，turquoise 的负向效应几乎不变。二者均未发生方向反转，95% CI 均未跨 0，BH FDR 均远低于 0.05，绝对效应衰减未超过 25%，排序也未跌出前 5。

## 全模块结果

- 16 个非灰模块的 current 与敏感性疾病系数相关为 `r=0.99758`。
- 对 current 中 `|β|≥0.01` 的 12 个模块，方向一致率为 `100%`。
- 若不设效应量门槛，16 个非灰模块中有 14 个方向一致。发生符号变化的 purple 和 midnightblue 在 current 和敏感性分析中均接近 0 且不显著，因此属于零附近的数值摆动，不构成稳定生物学方向的反转。
- blue、turquoise 及效应量前 7 的非灰模块排序均未改变。

## 预设升级判据

| 判据 | blue | turquoise | 结论 |
| --- | ---: | ---: | --- |
| 方向反转 | 否 | 否 | PASS |
| 95% CI 跨 0 | 否 | 否 | PASS |
| BH FDR ≥0.05 | 否 | 否 | PASS |
| 绝对效应衰减 >25% | 否 | 否 | PASS |
| current 前 3 降至前 5 之外 | 否 | 否 | PASS |

全模块辅助判据亦通过：`r=0.99758≥0.90`，重点模块方向一致率 `100%≥80%`。

## Go/no-go 决策

- **当前网络：保留为主分析。**
- **完整 WGCNA sensitivity：不启动。** 这是预设判据未触发后的主动决定，不是分析缺失。
- **T05：解锁。**
- **稿件边界：** 可以写“恢复全部 7 个年龄异常样本的固定模块敏感性未改变 Axis A/B 脑侧方向或效应排序”；不能写“年龄没有影响”“年龄异常值排除已被证明正确”或“已排除所有混杂”。

## 稿件建议表述

### 中文 Results

作为与结局相关样本选择的敏感性检验，我们冻结 current AD 网络的探针、8,807 个基因、模块归属及由 459 例学习的模块特征向量载荷，并将 7 个按组内 1.5×IQR 规则排除的年龄样本投影回固定模块；连接度异常样本仍予排除，敏感性队列共 466 例（156 例对照、310 例 AD）。AD blue 的年龄/性别调整疾病效应由 β=0.0518 变为 0.0528，AD turquoise 由 β=−0.0500 变为 −0.0500，二者方向、95% CI、FDR 及非灰模块绝对效应排序均保持。全部 16 个非灰模块的疾病系数高度一致（r=0.998），说明 Axis A/B 脑侧方向并非由这 7 例的排除所驱动。

### English Results

As a sensitivity analysis for outcome-related sample selection, we froze the probes, 8,807-gene feature space, module assignments and module-eigengene loadings of the current AD network and projected the seven samples excluded by the within-group 1.5×IQR age rule into the fixed modules. The connectivity outlier remained excluded, yielding 466 samples (156 controls and 310 AD cases). The age- and sex-adjusted disease effect changed from β=0.0518 to 0.0528 for AD blue and from β=−0.0500 to −0.0500 for AD turquoise; direction, 95% confidence intervals, FDR and absolute-effect ranks among non-grey modules were retained. Disease-effect estimates across all 16 non-grey modules were highly concordant (r=0.998), indicating that the brain-side directions of Axes A and B were not driven by exclusion of these seven samples.

## 补充图图注草案

**Supplementary Fig. Sx | Fixed-module sensitivity to restoration of seven GSE33000 age outliers.** **a**, Age- and sex-adjusted AD-status coefficients and 95% confidence intervals for the fixed AD blue and turquoise module eigengenes in the current age-filtered cohort (n=459) and after projection of the seven restored samples (n=466). **b**, Correspondence of AD-status coefficients across all 16 non-grey fixed modules; the dashed line denotes identity, and Axes A and B are highlighted. Module definitions, reference preprocessing parameters and eigengene loadings were learned only from the current cohort. P values were obtained from linear models adjusted for age and sex and corrected across 17 fixed modules using the Benjamini–Hochberg method. Source data are provided in `t04_figure_source_data_20260724.csv`.

## 主要产物

- 运行前合同：`docs/T04_preregistered_analysis_contract_20260724.md`
- 全报告：本文件
- 固定模块效应比较：`tables/t04_fixed_module_disease_effect_comparison_20260724.csv`
- 恢复样本表：`tables/t04_restored_age_outlier_samples_20260724.csv`
- Axis 与全模块触发器：`tables/t04_axis_upgrade_triggers_20260724.csv`、`t04_global_upgrade_triggers_20260724.csv`
- 验收与决策：`tables/t04_acceptance_checks_20260724.csv`、`t04_go_no_go_decision_20260724.csv`
- 补充图：`figures/t04_ad_age_outlier_fixed_module_robustness.{svg,pdf,tiff,png}`
- 图形源数据：`tables/t04_figure_source_data_20260724.csv`
- 投影参数与 eigengene：`t04_fixed_module_projection_bundle_20260724.rds`
- 输入/代码与输出哈希：`tables/t04_input_code_manifest_20260724.tsv`、`t04_output_hashes_20260724.tsv`
