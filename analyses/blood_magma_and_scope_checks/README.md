# Axis A 辅助证据核查：外周血、MAGMA、疾病范围（2026-10-19）

依据：修改意见 3.1–3.3；设计见 [DESIGN.md](DESIGN.md)，运行前写定（SHA-256 见 `DESIGN.sha256`）。**稿件未修改，候选结果。**

## 3.1 外周血（[blood_camera_and_effects.tsv](results/blood_camera_and_effects.tsv)）

**复现核对**：6 组重建的 limma t 值与 SYN01 冻结值逐基因相同，最大差 5.3 × 10⁻¹⁵（[blood_baseline_gate.tsv](qa/blood_baseline_gate.tsv)）。

| 数据集 | 疾病 | 样本 | n | fgsea NES（P） | CAMERA 方向 / P / BH | 受试者级差异，SD（95% CI） |
|---|---|---|---|---|---|---|
| GSE54837 | COPD | 全血 | 220 | 2.35（2.4 × 10⁻¹¹） | 上调 / 0.068 / 0.19 | 0.29（0.01–0.57） |
| E-MTAB-5278 | COPD | 全血 | 109 | 2.19（2.8 × 10⁻⁹） | 上调 / 0.11 / 0.19 | 0.32（−0.06 至 0.70） |
| GSE248493 | COPD | PBMC | 39 | 1.58（0.0017） | 上调 / 0.39 / 0.39 | 0.34（−0.30 至 0.98） |
| GSE63060 | AD | 全血 | 116 | 1.96（1.6 × 10⁻⁷） | 上调 / 0.13 / 0.19 | 0.02（−0.36 至 0.40） |
| GSE63061 | AD | 全血 | 112 | 2.53（4.5 × 10⁻¹⁵） | 上调 / 0.013 / 0.077 | 0.49（0.11–0.87） |
| GSE248417 | AD | 全血 | 97 | 1.81（1.0 × 10⁻⁴） | 上调 / 0.31 / 0.37 | 0.62（0.20–1.03） |

**解读**：
- 6 组方向一致，都是患者偏高。
- 考虑集合内基因相关性之后，证据明显减弱：CAMERA 估计的集合内平均相关为 0.08–0.13，只有 GSE63061 的 P < 0.05，BH 校正后没有一组 < 0.05。受试者级评分差异的 CI 在 3 组中不含 0。
- 两种方法的零假设不同。fgsea（预排序、基因置换）把基因当作相互独立；Axis A 成员彼此相关，所以 fgsea 的 P 值偏小。CAMERA 正是为此设计（Wu & Smyth 2012）。
- 因此血液结果只能写成“非配对外周队列中方向一致的表达偏高”，不能用 fgsea 的极小 P 值作为强证据。

**AddNeuroMed 重复供体**（[addneuromed_cross_dataset_matches.tsv](results/addneuromed_cross_dataset_matches.tsv)，[run_02.log](qa/run_02.log)）：
- 在 17,344 个共同基因上比较 GSE63060（329 个样本）与 GSE63061（388 个样本）。按预设规则标出 2 对，但它们的相关只有 0.59–0.62，与各自数据集内不同个体之间的最近邻相关相当（中位数 0.52–0.54，最大 0.80–0.85）。
- 第一对性别不同（男/女）、年龄不同（75/84 岁）；第二对年龄相差 18 岁。
- 结论：未发现跨数据集的重复供体。

## 3.2 MAGMA

**GWAS 文件**（来自 GWAS Catalog 协调文件的元数据）：

| 用途 | 编号 | 表型 | 祖源 | N | 基因组版本 |
|---|---|---|---|---|---|
| AD | GCST90704648 | 阿尔茨海默病（共识荟萃分析，不含大型生物库；medRxiv 10.1101/2025.10.20.25338060） | 欧洲 | 147,103 | GRCh38 |
| COPD | GCST90244098 | 机器学习肺量图推断的 COPD 易感性评分 | 欧洲 | 325,027 | GRCh38 |

**背景敏感性**（[magma_background_sensitivity.tsv](results/magma_background_sensitivity.tsv)）：原 4 项检验逐值重现。背景限定为能进入网络或交集分析的 4,832 个基因后，在新背景内重新做 Holm 校正：

| 检验 | 原背景：β（95% CI），P | 4,832 背景：β（95% CI），P |
|---|---|---|
| **AD × Axis A** | 0.36（0.20–0.53），7.2 × 10⁻⁶ | **0.46（0.25–0.66），6.4 × 10⁻⁶** |
| AD × Axis B | 0.05（−0.04 至 0.15），0.13 | −0.03（−0.15 至 0.09），0.70 |
| COPD × Axis A | 0.07（−0.13 至 0.28），0.25 | −0.08（−0.34 至 0.18），0.72 |
| COPD × Axis B | 0.01（−0.11 至 0.12），0.44 | −0.10（−0.25 至 0.05），0.91 |

结论不随背景改变：只有 AD × Axis A 富集。

**以 Patir 2019 核心小胶质基因为条件**（[magma_patir_conditional.tsv](results/magma_patir_conditional.tsv)）：
- 参照集为表 S3 中标记 Core 的 249 个基因：247 个映射到 Entrez，227 个进入基因分析，与 Axis A 重叠 61 个。
- 条件化后，Axis A 的 β = 0.30（SE 0.089），P = 4.0 × 10⁻⁴，两轴 Holm P = 7.9 × 10⁻⁴；无条件时 β = 0.36，P = 7.2 × 10⁻⁶。
- Patir 核心集本身也富集 AD 信号（P = 2.2 × 10⁻⁴），与 Axis A 同模型时为 P = 0.015。
- Axis A 中不属于 Patir 核心集的 74 个基因，单独检验仍富集（P = 0.0013）。
- 与 10-16 的 BRETIGEA 版本结论一致：Axis A 的 AD 富集不能由“属于一般小胶质基因集”完全解释。

## 3.3 疾病范围

**慢性过敏性肺炎（GSE150910）的批次**（[GSE150910_diagnosis_by_batch.tsv](results/GSE150910_diagnosis_by_batch.tsv)）：
- 15 个处理批次中，没有一个同时含病例和对照：82 例 CHP 全部在 2020 年 6–7 月的 10 个日期处理，103 名对照在 12 月的 4 个日期及 "other" 中处理。
- 测序板（plate）在两组间是均衡的。
- 结论：疾病与处理批次完全共线，CHP 与对照的差异无法与批次分开，加批次协变量或 ComBat 都不能补救。按修改意见，这项比较应移入补充材料，作为描述性结果，正文不再写成“同样由免疫介导的慢性过敏性肺炎也升高”。

**需要改的稿件句子**（只列出，未改稿；行号为 [论文全文_简体中文.md](../../06_resubmission_20261008/论文全文_简体中文.md) 当前版本）：

| 行 | 原文要点 | 建议 |
|---|---|---|
| 16（摘要） | “升高还见于慢性过敏性肺炎……而未见于特发性肺纤维化……” | 删去慢性过敏性肺炎；“未见于”改为“未检出统计学证据” |
| 89 | 慢性过敏性肺炎 0.63 SD 作为正文结果 | 移入补充材料，注明疾病与处理批次完全共线 |
| 89 | IPF、曾吸烟者“则未见升高”；精神分裂症、抑郁症“未见升高” | 改为“未检出统计学证据”，并给出效应与 CI |
| 91 | “也不太可能仅由吸烟引起” | 删除。只有一个曾吸烟者阴性队列，不能据此排除吸烟作用 |
| 117（讨论） | 同 16 | 同 16 |

## 运行

从项目根目录依次运行：

```
LC_ALL=en_US.UTF-8 Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/01_blood_camera_effects.R
LC_ALL=en_US.UTF-8 Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/02_addneuromed_followup.R
_workspace/envs/celloracle/bin/python -I 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/03b_extract_patir.py
LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/03_magma_background_patir.R
R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/04_chp_batch_crosstab.R
```

日志在 `qa/` 和 `logs/`。环境中没有 R 的 readxl 包，所以 Patir 表格先用 Python 提取为 `source_data/patir2019_core_symbols.txt`；原 Excel 文件也保存在 `source_data/`。
