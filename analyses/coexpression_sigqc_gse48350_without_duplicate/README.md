# GSE48350 共表达检查去除重复芯片后的重算（2026-10-30）

**状态**：已完成，结果已写入 S2 File 的 Table C。作者同意“都改”。

**原因**：
- 10-14 的 sigQC 共表达检查（`AxisA_sigQC_coherence_20261014`）在 GSE48350 中用了全部 41 张芯片。
- 自 10-19 起，复现分析已剔除 GSM1176257：它的 CEL 文件与 GSE11882 的对照 GSM350078 逐字节相同，见 S2 File 的 Table D (b)。
- 因此，Table D (b) 中“从全部分析中剔除”与 Table C 的样本数互相矛盾。

**做法**：
- 定义与 10-14 的 `03_summarise.R` 完全相同：
  - 观测值为 sigQC 的 autocor_median，即成员两两 Spearman 相关矩阵的中位数；10-14 的运行已核对它与 sigQC 官方输出相同。
  - 对照为 1,000 个按平均表达（20 个箱）匹配的随机集，随机种子为 20261014。
- 从缓存的 RMA 矩阵（`_workspace/sigqc_axisA_20261014/matrices.rds`）中去掉该样本。复现分支也是在预处理之后剔除它的。
- 其余 7 个数据集不变，直接沿用 10-14 的结果。

**结果**（`qa/old_vs_new_GSE48350.tsv`）

| 项目 | 41 张芯片 | 40 张芯片 |
|---|---|---|
| 一致性（autocor_median） | 0.214 | 0.206 |
| 随机集最大值 | 0.022 | 0.035 |
| 经验 P（单侧） | 0.001 | 0.001 |

- 成员层面，GSE48350 中位相关的新旧值 r = 0.993，最大变化 0.042。
- 2 个成员的分类改变：一致 99 → 97，不一致 37 → 39，覆盖不足 1 个。

**输出**：
- `results/set_level_coherence.tsv`、`results/gene_level_coherence.tsv`：8 个数据集，其中 GSE48350 为新值。
- `06_resubmission_20261008/code/supplement/update_parts_20261030.py`：把新值写入 Table C 的部件 S3_Table.xlsx。

**代码**：`code/01_rerun_GSE48350.R`，运行命令见文件开头。
