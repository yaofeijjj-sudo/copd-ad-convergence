# Axis A 的疾病情境范围（2026-10-03）

这一部分回答：Axis A 的升高是 COPD 与 AD 特有，还是同组织其他疾病也会出现？按作者决定，结果放在外周血部分（Figure 4 区域）。

**设计与变更记录：**
- [DESIGN.md](DESIGN.md)：主设计，运行前冻结，哈希见 `DESIGN.sha256`。
- [DESIGN_ADDENDUM_AD_reference.md](DESIGN_ADDENDUM_AD_reference.md)：AD 参照的检索、入选与合并规则，在下载新数据前写定，哈希见 `DESIGN_ADDENDUM_AD_reference.sha256`，需在项目根目录核验。
- [DEVIATIONS.md](DEVIATIONS.md)：全部偏离（第 1–14 项），以及每项确定于看结果之前还是之后。

## 复现（项目根目录）

```
_workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/01_parse_metadata.R
R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/02_run_specificity.R
R_LIBS=_workspace/r_libs_affy:_workspace/r_libs LANG=en_US.UTF-8 _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/04_ad_reference_meta.R
LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/03_draw_forest.R
LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_disease_specificity_redesign_20261003/code/05_draw_ad_meta.R
```

- **原始数据：** `02_data/raw/disease_specificity_tissue_20261003/`，取自 GEO 官方文件。另含 GSE48350 额上回 CEL、Ensembl 116 转录本→基因对应表。
- **R 包：** `_workspace/r_libs_affy` 只用于 GSE48350 的 RMA，见 DEVIATIONS 第 11 项。

## 文件

| 内容 | 文件 |
|---|---|
| 主图（Figure 4 区域） | `figures/Figure4_AxisA_disease_scope_forest.pdf` |
| AD 六数据集合并图（补充图） | `figures/SuppFigure_AD_reference_meta.pdf` |
| 对照疾病补充表 | `results/Supplementary_Table_AxisA_disease_scope.csv` |
| AD 参照补充表 | `results/Supplementary_Table_AD_reference_meta.csv` |
| 完整结果 | `results/axisA_specificity_results.tsv`、`results/AD_reference_meta.tsv`、`results/AD_reference_leave_one_out.tsv` |

两张图均为可编辑矢量 PDF（Arial 10 pt），文字碰撞检查为 0。

## 结果

### 目标疾病

**COPD：** GSE47460-GPL14550，145/91，0.44 SD（0.15 至 0.72），P = 0.0026；与正文独立队列分析的数值逐项一致。

**AD：** 6 个独立额叶数据集合并。

| 数据集 | 区域 | AD/对照 | 效应（SD） | 95% CI | P |
|---|---|---|---|---|---|
| GSE84422 | DLPFC | 41/16 | 0.15 | −0.49 至 0.79 | 0.64 |
| GSE118553 | 额叶皮质 | 40/21 | 0.46 | −0.13 至 1.04 | 0.12 |
| GSE313307 | 前额叶皮质（RNA-seq，调整临终感染） | 113/89 | 0.33 | 0.04 至 0.62 | 0.024 |
| GSE140841 | BA9（RNA-seq） | 37/27 | 0.79 | 0.18 至 1.40 | 0.012 |
| GSE48350 | 额上回 | 21/20 | 0.92 | 0.35 至 1.49 | 0.0024 |
| GSE36980 | 额叶皮质 | 15/18 | 0.07 | −0.84 至 0.99 | 0.87 |
| **合并（随机效应）** | | **267/191** | **0.46** | **0.22 至 0.70** | **0.00017** |

- **方向：** 六个数据集均为正。
- **异质性低：** I² = 21%，异质性 P = 0.28。
- **逐一剔除：** 合并效应在 0.37–0.53 之间，P 均 < 0.003，结论不依赖任何单个数据集。

### 对照疾病（家族内 BH 校正）

| 组织 | 对照疾病 | 数据集 | 病例/对照 | 效应（SD） | 95% CI | FDR | 判定 |
|---|---|---|---|---|---|---|---|
| 肺 | 慢性过敏性肺炎 | GSE150910 | 81/103 | 0.64 | 0.28 至 0.99 | 0.0015 | 检出升高 |
| 肺 | IPF | GSE150910 | 102/103 | −0.14 | −0.44 至 0.15 | 0.52 | 未检出 |
| 肺 | 曾吸烟 vs 从不吸烟（无肺病） | GSE150910 | 51/45 | 0.09 | −0.32 至 0.50 | 0.66 | 未检出 |
| 脑 | 亨廷顿病 | GSE64810 BA9 | 19/49 | 0.80 | 0.30 至 1.30 | 0.011 | 检出升高 |
| 脑 | 帕金森病 | GSE68719 BA9 | 29/44 | 0.35 | −0.04 至 0.75 | 0.13 | 未检出 |
| 脑 | 精神分裂症 | GSE53987 BA46 | 15/19 | 0.18 | −0.55 至 0.92 | 0.61 | 未检出 |
| 脑 | 双相障碍 | GSE53987 BA46 | 17/19 | −0.67 | −1.19 至 −0.15 | 0.032 | 检出降低 |
| 脑 | 抑郁症 | GSE53987 BA46 | 17/19 | −0.36 | −0.99 至 0.27 | 0.32 | 未检出 |

### 次要情境（只作描述）

- **GSE47460 内部：**
  - IPF/UIP vs 对照：两个平台均接近零；
  - 其他间质性肺病：0.47（P = 0.005）与 0.52（P = 0.22）；
  - COPD 直接高于 IPF/UIP：0.44 与 0.91。
- **GSE33000 亨廷顿病：** 1.28 SD。对照参与过 AD 建网，不作独立证据。
- **GSE10072 癌旁肺组织，当前吸烟 vs 从不吸烟：** 0.46（−0.48 至 1.40）。

### 注释性读出（不改变判定）

- **cameraPR：** 除效应接近零的 GSE47460-GPL14550 IPF 外，方向与主要读出一致；6 个 AD 数据集均为上调。
- **与表达匹配随机基因集比较的经验 P：**

| 情境 | 经验 P |
|---|---|
| COPD（GPL14550） | 0.038 |
| COPD（GPL6480） | 0.002 |
| 慢性过敏性肺炎 | 0.12 |
| 亨廷顿病 | 0.062 |
| 帕金森病 | 0.030 |
| 双相障碍 | 0.10 |
| AD，GSE140841 | 0.010 |
| AD，GSE48350 | 0.008 |
| AD，GSE313307 | 0.22 |
| AD，其余数据集 | 0.19–0.95 |

## 按预先写定的规则解释

- **肺：** 只有慢性过敏性肺炎检出升高。对应规则："与免疫介导的肺部炎症相关，而非纤维化本身"。吸烟暴露未检出。
- **脑：** AD（合并）与亨廷顿病升高；精神疾病均未检出升高，双相障碍为降低。对应规则："多种神经退行性疾病共有的髓系反应，在非退行性精神疾病中未检出"。
- **综合表述：** 在已检验的情境中，Axis A 升高见于 COPD、慢性过敏性肺炎、AD（6 个独立额叶数据集合并）和亨廷顿病；在 IPF、无肺病吸烟者、精神分裂症和抑郁症中未检出升高，在双相障碍中降低。因此，Axis A 代表免疫介导的慢性肺部炎症与神经退行所共有的髓系反应，COPD 肺与 AD 脑是其中方向一致的两例。

## 解释时必须保留的限制

- **AD 参照的变更过程：** 冻结设计中的 AD 参照是 GSE84422。看到它（0.15 SD）和 GSE118553（0.46 SD）的结果后，作者要求先改用 GSE118553，再联网补充数据集。最终采用在下载新数据前写定的入选标准，纳入全部 6 个合格数据集合并，GSE84422 也在其中。变更顺序如实记录于 DEVIATIONS 第 10–14 项。
- **跨数据集比较：** 平台、区域、诊断方式各不相同。I² 低，但只比较方向，不比较效应大小。
- **潜在供者重叠：** GSE313307 与 GSE118553 均来自英国脑库，可能有供者重叠，无法核验。
- **数据处理：**
  - GSE48350 从原始 CEL 重新做 RMA；对照限于不小于 AD 组最小年龄（74 岁）者，剔除 28 名较年轻对照。
  - GSE140841 由 kallisto 转录本计数汇总，92% 的计数映射到 Entrez。
- **过敏性肺炎与处理批次混杂：** 见 DEVIATIONS 第 1 项。
- **共享对照：** 亨廷顿病与帕金森病两项对比共享 44 名对照。
- **死亡情况：** 精神疾病与神经退行性疾病供者的死亡情况不同。GSE313307 已调整临终呼吸道感染，其余数据集无此信息。
- **未检出不等于不存在：** 精神疾病组每组 15–19 人，置信区间较宽。
- **正文的 AD 独立队列仍为 GSE84422：** Figure 1 和正文独立队列复现用的是 GSE84422，此处 AD 参照为 6 个数据集合并，写稿时要分别说明。

## 核验

- 输入哈希：`qa/input_sha256.tsv`、`qa/input_sha256_AD_reference.tsv`。
- 运行日志：`qa/run.log`、`qa/run_AD_reference.log`。
- 环境：`qa/sessionInfo.txt`、`qa/sessionInfo_AD_reference.txt`。
- 样本规则：`qa/sample_groups.tsv`、`qa/AD_reference_sample_rules.tsv`。
- 加入 GSE118553 后，对照疾病各项的数值与第一轮完全一致；第一轮结果备份在 `qa/axisA_specificity_results_before_AD_reference_change.tsv`。
- GSE84422 与 GSE118553 在合并脚本中的效应值与主脚本一致（0.149、0.455）。
- README 旧版本备份在 `qa/README_before_AD_meta.md`。
