# Axis A 在独立数据集中是否仍是一个协同变化的基因集（sigQC，2026-10-14）

**状态：候选分析，已完成；未写入论文。** 设计在运行前冻结：[DESIGN.md](DESIGN.md)（SHA-256 见 `DESIGN.sha256`，运行后核验一致）。

## 问题

Axis A 固定为 137 个 Entrez 基因。此前已有：来源模块（肺 red、脑 blue）的结构保留（Zsummary），以及 Axis A 评分与疾病的独立复现。缺的是：**这 137 个基因本身**在未参与定义它的数据中是否仍一起变化，以及哪些成员稳定。

## 方法

- sigQC 0.1.24（Dhawan 2019 *Nat Protoc* 14:1377），CRAN 官方最终版；该包于 2025-12-19 因依赖 biclust 被归档，从 CRAN 归档安装到独立库 `_workspace/r_libs_sigqc`（安装脚本与日志：`_workspace/sigqc_install_20261014/`）。
- 8 个独立数据集，样本和预处理与复现分析（`AxisA_singscore_official_rerun_20261004`）完全一致；`qa/matrix_build_check.tsv` 核对病例数、对照数和 Axis A 可测数，8/8 一致。
- 主指标：sigQC 紧密度 `autocor_median`（成员间 Spearman 相关矩阵的中位数）。主比较：表达匹配随机基因集 1,000 次（20 个平均表达分箱）；另报告 sigQC 自带的随机基因对照与置换对照（各 200 次）。
- 基因层面：每个成员与其他成员 Spearman 相关的中位数，逐数据集计算；与 sigQC 输出的相关矩阵逐元素一致（最大差 5.6 × 10⁻¹⁶）。
- 未使用疾病标签（RNA-seq 低表达过滤沿用复现分析的 `filterByExpr` 分组参数，仅为复现同一矩阵）。

## 结果

### 集合层面：8/8 个数据集中 Axis A 均成团

| 数据集 | 组织 | 样本 | 可测成员 | 紧密度 | 表达匹配随机：均值 / 95% 分位 | 经验 P |
|---|---|---:|---:|---:|---|---:|
| GSE47460-GPL14550 | 肺 | 236 | 137 | 0.379 | −0.003 / 0.002 | 0.001 |
| GSE47460-GPL6480 | 肺 | 92 | 134 | 0.361 | 0.000 / 0.007 | 0.001 |
| GSE84422 | DLPFC | 57 | 135 | 0.201 | 0.043 / 0.066 | 0.001 |
| GSE118553 | 额叶皮质 | 61 | 135 | 0.047 | 0.003 / 0.007 | 0.001 |
| GSE313307 | 前额叶 RNA-seq | 202 | 119 | 0.392 | 0.011 / 0.027 | 0.001 |
| GSE140841 | BA9 RNA-seq | 64 | 127 | 0.370 | 0.025 / 0.038 | 0.001 |
| GSE48350 | 额上回 | 41 | 135 | 0.214 | −0.002 / 0.008 | 0.001 |
| GSE36980 | 额叶皮质 | 33 | 132 | 0.443 | 0.024 / 0.046 | 0.001 |

经验 P = 0.001 为 1,000 次抽样的下限（观测值高于全部随机值）。sigQC 自带随机基因与置换对照的 97.5% 分位均在 0.006–0.096，全部低于观测值。三种评分方式（平均、中位数、第一主成分）高度一致：平均 vs 中位数 ρ = 0.80–0.98。GSE118553 的紧密度最低（0.047），虽高于随机，但明显弱于其他数据集。

### 基因层面：99 个一致，37 个不一致，1 个覆盖不足

- **一致（99）**：在可测的每个数据集中，与其他成员的中位相关都 > 0，且至少在 6 个数据集中可测。最稳定的是髓系通用基因（NCKAP1L、LAPTM5、ITGB2、SYK、CD68、FCER1G，跨数据集中位相关约 0.5）。补体、TREM2 通路和 AD 风险成员全部为一致：TREM2、C1QA、C1QB、FCER1G、SYK、HCK、LYN、CSF1R、CD33、MS4A6A、APOC1、C3AR1。
- **不一致（37）**：其中 23 个只在 1 个数据集中 ≤ 0，且大多是 GSE118553。在 ≥ 2 个数据集中 ≤ 0 的有 14 个：SPIB、SEL1L3、DNASE2B、SLC46A3、IGLL1、TRAF1、ASIP、UBE2C、PLA2G7、SERPINF2、CLEC10A、CD209、SERINC2、IL16。
- **覆盖不足（1）**：CARD9（8 个数据集中仅 5 个可测）。

这些标签只是描述，不重新定义 Axis A；所有评分仍使用 137 个基因。

### 解读边界

组织数据中成员间的相关，部分来自样本间髓系细胞含量的差异。因此结果支持“Axis A 在独立组织中是一个共同变化的髓系程序”，不支持“单个细胞内的共调控”。脑 RNA-seq 中 sigQC 的“表达高于数据集中位数的比例”很低，反映髓系基因在全脑组织中本身表达低，不影响紧密度结论。

## 偏离与技术说明

1. 首次运行时，sigQC 紧密度步骤缺少 `gtools`，随后又因非交互设备关闭了显示列表，导致 `gridGraphics` 截图返回 NULL；sigQC 在 tryCatch 中把 `autocor_median` 记为 0。处理：补装 `gtools`/`caTools`，并在 `02_run_sigqc.R` 中令 `dev.new()` 打开的设备启用显示列表；删除不完整输出后全部重跑。sigQC 代码未改动。
2. sigQC 的可选包 RankProd 需要系统 GMP 库，未安装；sigQC 因此只跳过跨数据集秩积表。基因层面的跨数据集汇总由预设的 `gene_level_coherence.tsv` 提供。
3. 原 `02_run_sigqc.R` 拆分为 `02_run_sigqc.R`（只运行 sigQC）和 `03_summarise.R`（汇总、表达匹配随机和基因层面），计算定义未变。

## 文件

- `results/set_level_coherence.tsv`：集合层面紧密度及表达匹配随机
- `results/gene_level_coherence.tsv`：每个成员在 8 个数据集中的中位相关及分类
- `results/sigqc_radar_metrics.tsv`、`results/sigqc_control_summaries_autocor.tsv`：sigQC 官方指标与对照
- `sigqc_output/`：sigQC 原始输出（雷达图、相关热图、对照箱线图等 PDF）
- `qa/`：输入哈希、矩阵核对、运行日志、sessionInfo

## 复现（项目根目录）

```sh
R_LIBS=_workspace/r_libs_affy:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_coherence_20261014/code/01_build_matrices.R
R_LIBS=_workspace/r_libs_sigqc:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_coherence_20261014/code/02_run_sigqc.R
R_LIBS=_workspace/r_libs_sigqc:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_sigQC_coherence_20261014/code/03_summarise.R
```
