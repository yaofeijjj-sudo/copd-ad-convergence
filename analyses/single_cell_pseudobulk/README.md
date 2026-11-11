# Axis A 供体伪 bulk 疾病比较（2026-10-14，候选）

**结论：两个器官中，髓系细胞内部的 Axis A 整体都没有随疾病显著改变。** 肺髓系在 COPD 中偏低（P = 0.12），脑 Microglia 2 随局部病理略偏高（P = 0.32）。这与图 5 的 UCell 结果方向和结论一致。按事先写定的解释规则：在当前样本量下，没有证据表明组织层面 Axis A 的升高来自细胞内表达升高。本分析不检验细胞组成。

设计在运行前冻结，见 [DESIGN.md](DESIGN.md)。本分析不是盲法：图 5 和 2026-09-25 诊断性核对的方向已知。

## 主要结果（[fry_summary.tsv](results/fry_summary.tsv)）

| 器官 | 细胞类别 | 供体 | 模型 | 可检验成员 | fry 方向 | P | 成员上调 / 下调 | 成员 log2FC 中位数 |
|---|---|---:|---|---:|---|---:|---:|---:|
| 肺 GSE302339 | Myeloid | 33 COPD / 23 对照 | 主要：+ 年龄、性别、吸烟状态 | 135 | Down | 0.115 | 39 / 96 | −0.079 |
| 肺 | Myeloid | 同上 | 次要：未调整 | 135 | Down | 0.134 | 40 / 95 | −0.090 |
| 脑 SEA-AD DFC | Micro-PVM_2（Microglia 2） | 79 | 主要：每 1 SD 局部 CPS，+ 年龄、性别、PMI、RIN、研究 | 120 | Up | 0.316 | 72 / 48 | +0.020 |

- 每个器官一个主要检验，按稿件做法不跨器官校正；两器官 BH 参考值为 0.23 和 0.32。
- 两侧都没有任何 Axis A 成员达到全基因组 FDR < 0.05。逐基因结果见 [axisA_member_logFC.tsv](results/axisA_member_logFC.tsv)，仅作描述。
- `fry` 的 `PValue.Mixed`（不区分方向，检验是否有成员发生任一方向的变化）不是事先指定的读出：肺 0.37，脑 0.047。它不说明整组朝某一方向变化，不作解释。
- 与图 5 的对应：肺 UCell 中位差为 −0.016（FDR 0.062），脑 β = +0.0012/SD（FDR 0.42）。两种读出的方向一致，结论都是阴性。

## 方法

- 伪 bulk：同一供体、同一细胞类别的原始 UMI 直接求和。肺从 65 个 GEO 原始 10x 矩阵中，按 SYN01 已评分、通过 QC 的 `Myeloid` 细胞求和（116,558 个细胞；同一供体多个肺叶文库合并），基因按 Ensembl ID 汇总。脑复用 2026-09-25 的计数表：与图 5 UCell 相同的 23,198 个 Micro-PVM_2 细胞核，按原 h5ad 的 `gene_ids` 一一映射回 Ensembl ID。
- 检验：edgeR v4 官方流程（Chen 等 2025 *NAR*，DOI 10.1093/nar/gkaf018），依次为 `filterByExpr(design)`、TMM（`normLibSizes`）、`estimateDisp`，然后用 `fry` 对预定义基因集做自包含方向性检验。逐基因 log2FC 来自 `glmQLFit` + `glmQLFTest`。
- 基因集：当前 Axis A 137 基因中，在该数据里可测的 136 个，即 SYN01 的 `lung_axis_measured.txt` / `brain_axis_measured.txt`；表达过滤后剩余的成员数见上表。
- 软件：R 4.5.3，edgeR 4.8.2，limma 3.66.0（[会话信息](logs/R_sessionInfo.txt)）。

## 核验（[verification.json](qa/verification.json)）

- 细胞和供体数与 S3 Table 一致：肺 116,558 个细胞 / 56 名供体，脑 23,198 个细胞核 / 79 名供体。
- 肺侧每名供体的伪 bulk 总 UMI，与细胞 QC 时独立记录的每细胞总 UMI 加和完全相等（最大差值为 0）。
- 脑计数表的 SHA-256 与 2026-09-25 的记录一致（b997c1f9…）。
- `fry` 与 99,999 次旋转的 `roast` 对照：肺 Down，0.115 vs 0.116；脑 Up，0.316 vs 0.315（[fry_vs_roast.tsv](qa/fry_vs_roast.tsv)）。
- 首次运行时 `fry` 要求先估计离散度，因此按官方要求补上 `estimateDisp`，这一步发生在看到任何结果之前，记录在 DESIGN.md。之后为结果表补充基因符号列，重跑后统计结果逐字节一致。

## 局限

- 肺 COPD 组年龄更大、男性更多、包年更多；吸烟状态与诊断部分重叠。
- 肺 `Myeloid` 是混合类别，内部亚群比例的变化也会表现为伪 bulk 表达变化。
- 脑侧是单细胞核数据，可能漏检部分小胶质细胞激活基因。
- “组织升高来自髓系细胞增多”仍是推测，需要组成分析（例如 propeller）才能检验。

## 文件

- 代码：[01 肺伪 bulk](code/01_build_lung_pseudobulk.R)、[02 脑基因映射](code/02_brain_gene_map.py)、[03 edgeR + fry](code/03_edger_fry.R)、[04 核验](code/04_verify.R)
- 计数与供体：[source_data/](source_data/)、[results/lung_donors_used.tsv](results/lung_donors_used.tsv)、[results/brain_donors_used.tsv](results/brain_donors_used.tsv)
- 运行命令：从项目根目录执行，并设置 `LC_ALL=en_US.UTF-8`，依次运行 `Rscript …/code/01_build_lung_pseudobulk.R`、`_workspace/envs/celloracle/bin/python -I …/code/02_brain_gene_map.py`、`Rscript …/code/03_edger_fry.R`、`Rscript …/code/04_verify.R`。
