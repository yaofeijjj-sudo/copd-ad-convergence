# 肺单细胞文库层面质控核查（2026-10-28，已完成，已写入正文与 S6 Table d–f）

## 做了什么

现行质控在每个文库内按 3 倍 MAD 剔除离群细胞（scuttle quickPerCellQC）并以 scDblFinder 剔除双细胞。若整个文库质量都差，文库内的相对阈值无法识别。按 OSCA（Amezquita 等 2020）比较各批次质控结果的建议，在 65 个文库之间比较质控后线粒体比例的中位数：对数尺度上高于全部文库中位数 3 倍 MAD（阈值 17.5%）的文库标记为离群。

## 结果

- 标记 2 个文库：GSM9102301（ua72，中位数 67.7%）和 GSM9102263（copd4 的两个文库之一，20.7%）。
- 髓系细胞供体层面 COPD vs 对照（与 `lung_rank_and_scores.R` 相同的聚合方式：文库中位数 → 供体内中位数，≥20 个细胞），现行 Axis A（137 个基因）：全部文库复现 S6 Table a（56 名供体，中位数差 −0.0159，P = 0.0103）；去除 2 个文库后 55 名供体，差值 −0.0156，P = 0.0137（[results/myeloid_donor_test_library_check.tsv](results/myeloid_donor_test_library_check.tsv)）。结论不变。
- 伪 bulk 检验未重跑。

## 10-28 更正

`code/run.R` 的髓系检验读取的是 9 月旧版 Axis A（135 个基因）的细胞评分（`细胞定位/analysis/06_single_cell_projection/results/GSE302339_axisA_official_ucell_20260903_candidate/`）。其校验只比对了 P 值，而秩检验的 P 在新旧版本间相同，所以没有发现；中位数差 −0.0154 与 S6 Table a 的 −0.0159 不一致。`code/run_axisA137.R` 改用现行 SYN01 的文库中位数（`细胞定位/SYN01_Entrez_downstream_20261001/results/lung/tables/axisA_library_broad_medians.tsv`）重跑，并校验全部报告值。旧结果保留为 `results/myeloid_donor_test_library_check_old135gene_superseded.tsv`。文库标记（`library_qc_with_flag.tsv`）和 CellTypist 标签只取决于质控指标，与基因集无关，沿用 `run.R` 的输出。S6 Table e 已按新结果更新。

## 复现

```sh
Rscript 03_candidate_analysis/Lung_sc_library_QC_check_20261028/code/run.R          # 文库标记与 CellTypist 标签
Rscript 03_candidate_analysis/Lung_sc_library_QC_check_20261028/code/run_axisA137.R # 髓系供体检验（现行 137 个基因）
```
