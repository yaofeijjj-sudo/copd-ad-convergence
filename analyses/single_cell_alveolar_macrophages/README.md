# 肺泡巨噬细胞内 Axis A 的供体层面疾病检验（2026-10-28，候选结果，未入稿）

- 依据：修改意见 2.4。
- 设计：[DESIGN.md](DESIGN.md)，运行前写定，SHA-256 见 `DESIGN.sha256`。
- 独立单位：供体（56 名，COPD 33 名，对照 23 名）。
- 细胞群：CellTypist 标签为 `Alveolar macrophages` 的细胞，共 64,542 个，与 CellOracle 肺模型的来源细胞群相同。

## 结论

**肺泡巨噬细胞内部未检出 Axis A 随 COPD 的升高，两种读出的方向都是偏低。** 按设计第 6 节规则 1，肺侧"细胞内未检出升高"成立。

- **伪 bulk**：fry 为下调，P = 0.25；CAMERA 为下调，P = 0.13。可测成员 134 个，其中 85 个 logFC < 0。
- **供体 UCell**：COPD 患者的评分低于对照。供体中位数分别为 0.170 和 0.189，Hodges–Lehmann 位置差为 −0.023（95% CI −0.044 至 −0.001），Wilcoxon P = 0.036。
- **两种读出显著性不同**（设计规则 4）。核查结果显示，差别不是协变量造成的：在供体 UCell 中位数上用与伪 bulk 相同的协变量做线性回归，调整后为 −0.72 SD（95% CI −1.33 至 −0.11），P = 0.021。差别来自检验对象不同：
  - UCell 在每个细胞内按基因表达排序计分，反映 Axis A 基因相对其他基因的位次；
  - fry 检验成员 logFC 整体是否偏离 0（自包含检验）；
  - CAMERA 比较 Axis A 成员与其他基因（竞争性检验）。

  两种读出方向一致，只是显著性不同。
- **与全部髓系细胞的结果相同**：髓系 UCell 的 P = 0.010，伪 bulk 的 P 为 0.12 和 0.056，方向都是偏低。
- **髓系内部的组成不能解释这一结果**：肺泡巨噬细胞占髓系细胞的比例，COPD 与对照的中位数分别为 0.48 和 0.48（四分位距 0.38–0.60 和 0.29–0.66）。这一项只作描述，未作检验。

## 结果表

| 分析 | 供体 | 伪 bulk fry | 伪 bulk CAMERA | UCell HL 位置差（95% CI） | Wilcoxon P |
|---|---|---|---|---|---|
| **主分析：肺泡巨噬细胞** | 56 | 下调，0.25 | 下调，0.13 | −0.023（−0.044 至 −0.001） | 0.036 |
| a. 每名供体 ≥ 50 个细胞 | 54 | 下调，0.13 | 下调，0.12 | −0.023（−0.045 至 −0.001） | 0.044 |
| b. 去除 2 个线粒体离群文库 | 55 | 下调，0.34 | 下调，0.25 | −0.022（−0.043 至 0.000） | 0.051 |
| c. 肺泡巨噬细胞家族（4 个标签） | 56 | 下调，0.24 | 下调，0.19 | −0.023（−0.044 至 −0.002） | 0.034 |
| 参照：全部髓系细胞（S6 Table） | 56 | 下调，0.115 | 下调，0.056 | −0.021（−0.042 至 −0.006） | 0.010 |

文件：
- [pseudobulk_fry_camera.tsv](results/pseudobulk_fry_camera.tsv)
- [donor_ucell_wilcoxon.tsv](results/donor_ucell_wilcoxon.tsv)
- [donor_ucell_lm_covariate_check.tsv](results/donor_ucell_lm_covariate_check.tsv)
- [axisA_member_logFC.tsv](results/axisA_member_logFC.tsv)
- [AM_share_of_myeloid_by_donor.tsv](results/AM_share_of_myeloid_by_donor.tsv)

## 核验

- **原始计数**：用同一读取流程把全部髓系细胞按供体相加，结果与 10-14 冻结的伪 bulk 逐元素相同（`code/01_build_pseudobulk.R` 中有断言）。
- **参照行复现 S6 Table**：
  - fry P 与 0.1151 的差为 2 × 10⁻¹²；
  - CAMERA P 与 0.0562 的差为 6 × 10⁻¹²；
  - Wilcoxon P 与表中 0.01032 的差为 4 × 10⁻⁶。这是表中数值的舍入误差，断言阈值为 5 × 10⁻⁵，DESIGN 写的是 1e-6。

## 对稿件的含义（待作者决定是否写入）

- 结果"细胞内"一节的肺部分，可以在髓系细胞结果后补一句肺泡巨噬细胞的结果，写法为"未检出升高、方向偏低"，并写明 UCell 与伪 bulk 的 P 值。
- 髓系 UCell 检验新增了 Hodges–Lehmann 区间：−0.021（−0.042 至 −0.006）。这正好补上了修改意见 2.4 指出的"阴性结果缺区间"。
- S6 Table 可以加一个子表，放上面的结果表。

## 运行

```sh
_workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_AM_pseudobulk_20261028/code/01_build_pseudobulk.R
_workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_AM_pseudobulk_20261028/code/02_tests.R
_workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_AM_pseudobulk_20261028/code/03_ucell_covariate_check.R
```

R 4.5.3，edgeR 4.8.2，limma 3.66.0。会话信息见 `logs/R_sessionInfo.txt`。
