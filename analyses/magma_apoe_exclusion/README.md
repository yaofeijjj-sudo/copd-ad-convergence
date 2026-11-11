# Axis A AD MAGMA：排除 APOE 区域的敏感性分析（2026-10-16）

回应替代解释 R7 的一半（“富集完全由 APOE 区域驱动”）；另一半（以一般小胶质基因集作条件）仍待做。

## 设计

- 问题：现行 Axis A（137 个基因）在 AD GWAS 中的 MAGMA 竞争性富集，去掉 APOE 区域后是否仍然存在。
- 统计单位：基因。检验与主分析相同：MAGMA v1.10 竞争性基因集检验，单侧正向，按默认设置校正基因大小、密度、逆 MAC 和样本量。
- 冻结输入：SYN01 分支的 `AD_total_N.genes.raw`（GCST90704648，不含生物样本库）和 `locked_axes.sets`（Axis A 137 / Axis B 434）；NCBI37.3 基因坐标。
- 区域：APOE（Entrez 348）基因边界上下游各 1 Mb，GRCh37 chr19:44,409,039–46,412,650。定义沿用 2026-07-15 预设方案 §8 和 `locus_exclusion_region_manifest_20260715.csv`，运行前已固定。Jansen 2019 *Nat Genet*（[DOI](https://doi.org/10.1038/s41588-018-0311-9)）在 MAGMA 中使用的较窄 APOE 区域（chr19:45,020,859–45,844,508）完全落在这个窗口内。
- 两种做法：
  - V1 只移出基因集：区域内的 Axis 基因从基因集中去掉，背景不变（7 月方案的实现方式）。
  - V2 移出全部分析：区域内 69 个基因用 `--settings gene-exclude` 从整个分析中去掉（MAGMA v1.10 手册中的 `--settings gene-include/gene-exclude`）。
- 校正：每种做法内对两个轴作 Holm 校正，与 7 月方案的 variant family 相同。

## 结果

区域内的 Axis 基因：Axis A 只有 APOC1（基因 P = 1.1 × 10⁻¹⁴）；Axis B 有 ZNF222、ZNF226。APOE 本身不在任一轴中。

| 做法 | 轴 | 基因数 | β | SE | P（单侧） | Holm P（两轴） |
|---|---|---|---|---|---|---|
| 完整基因集（SYN01 参照） | Axis A | 130 | 0.364 | 0.084 | 7.18 × 10⁻⁶ | 1.44 × 10⁻⁵ |
| V1 只移出基因集 | Axis A | 129 | 0.336 | 0.084 | 3.22 × 10⁻⁵ | 6.44 × 10⁻⁵ |
| V2 移出全部分析 | Axis A | 129 | 0.343 | 0.083 | 1.78 × 10⁻⁵ | 3.56 × 10⁻⁵ |
| 完整基因集（SYN01 参照） | Axis B | 409 | 0.052 | 0.047 | 0.132 | 0.132 |
| V1 只移出基因集 | Axis B | 407 | 0.049 | 0.047 | 0.149 | 0.149 |
| V2 移出全部分析 | Axis B | 407 | 0.053 | 0.047 | 0.127 | 0.127 |

参照行的 Holm P 是两轴内校正。正文中的 2.87 × 10⁻⁵ 是 AD/COPD × 两轴共四项检验的校正结果。

结论：排除 APOE 区域后，Axis A 的 β 下降约 6–8%，两种做法下仍通过校正；Axis B 仍未检出。这与 7 月 135 基因版的结果方向一致（β 0.358 → 0.330，P 1.27 × 10⁻⁵ → 5.59 × 10⁻⁵）。这只是关联层面的敏感性分析，不涉及因果方向。

## 文件

- 代码：[code/run_apoe_exclusion.py](code/run_apoe_exclusion.py)（用 `python3 -I` 运行）
- 结果表：[results/AD_MAGMA_APOE_exclusion.csv](results/AD_MAGMA_APOE_exclusion.csv)
- 区域基因清单：[inputs/APOE_1Mb_region_gene_manifest.csv](inputs/APOE_1Mb_region_gene_manifest.csv)；排除列表：`inputs/APOE_1Mb_region_genes.entrez.txt`
- MAGMA 输出：`magma/`；运行日志：`logs/`；输入哈希：[results/provenance.json](results/provenance.json)

## 稿件

2026-10-16 已写入中文稿（`06_resubmission_20261008/论文全文_简体中文.md`），结果遗传段和方法 MAGMA 段各一句。正文只报告 V2（“排除区域”的字面做法），V1 记在 `paper/潜在暴雷点汇总_20261016.md` B5.4，留作回复材料。英文 Word 尚未同步，拟用措辞：*The enrichment persisted after excluding the APOE region (±1 Mb), the strongest AD risk locus (P = 1.78 × 10⁻⁵); APOC1 was the only Axis A gene in this region.*
