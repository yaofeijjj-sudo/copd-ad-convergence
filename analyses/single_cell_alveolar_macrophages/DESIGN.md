# 肺泡巨噬细胞内 Axis A 的供体层面疾病检验：设计（2026-10-28，运行前写定）

依据：修改意见 2.4，"肺侧以肺泡巨噬细胞为主要对象"。作者于 2026-10-28 授权。

## 1. 问题与估计量

- **问题**：在 COPD 患者与对照的肺泡巨噬细胞内部，Axis A 的表达是否不同。
- **为什么要做**：现有的肺侧细胞内检验，单位是全部髓系细胞（116,558 个，包括肺泡巨噬细胞、单核细胞、单核来源巨噬细胞和树突状细胞）。GSE302339 的原作者报告，肺气肿型 COPD 中巨噬细胞比例下降、CD14 单核细胞比例上升（Zhang 2025）。如果髓系内部的亚群比例变了，混合类别的伪 bulk 结果就会同时反映亚群组成和细胞内表达两种变化。按细胞类型分别汇总伪 bulk，是 Squair 2021 给出的做法。CellOracle 的肺侧模型用的也是同一群肺泡巨噬细胞。
- **估计量**：
  - 伪 bulk：调整年龄、性别和吸烟状态后，COPD 相对对照的 Axis A 集合方向及其 P 值，分别用 fry（自包含检验）和 CAMERA（竞争性检验）。
  - 供体 UCell：先取每个文库内细胞评分的中位数，再取每名供体的中位数；比较 COPD 与对照，报告中位数差、Wilcoxon P，以及 Hodges–Lehmann 位置差及其 95% CI。

## 2. 独立单位与细胞群

- **独立单位**：供体。GSE302339 每个文库只对应一名供体。
- **主分析细胞群**：CellTypist（Human_Lung_Atlas，best match）标签为 `Alveolar macrophages` 的通过质控的细胞，共 64,542 个，与 CellOracle 肺模型的来源细胞群相同。
- **供体纳入标准**：每名供体至少 20 个细胞，沿用原标准。可行性检查显示 56 名供体均满足（33 名 COPD，23 名对照）。

## 3. 冻结输入（SHA-256）

- 细胞表 `细胞定位/SYN01_Entrez_downstream_20261001/results/lung/objects/axisA_cell_scores.tsv.gz`：e7adbf66ce14d56a28c41b7fada835762857d2ee977ae9813b7115d18900738c。包含 barcode、供体、CellTypist 标签、现行 137 基因的 UCell 评分和协变量。
- 可测 Axis A 成员 `细胞定位/SYN01_Entrez_downstream_20261001/results/lung_axis_measured.txt`（136 个 Ensembl ID）：f5728b9b22607e2eb35753e81ba516167980660b5e652fe86a9d93c314ffb661。
- 样本台账 `01_current_analysis/.../GSE302339_axisA_official_ucell_20260903_candidate/tables/sample_ledger.tsv`：d54ca41729913c6f84df0ed2ef081982bb6ff540ea703f49506b777c82632d89。
- 文库标记 `03_candidate_analysis/Lung_sc_library_QC_check_20261028/results/library_qc_with_flag.tsv`：5169fe1263cfb9a6437bb81b6495c6624fc24081f175846199c9705c90813e59。
- 原始计数：`_workspace/gse302339_axisA_20260903/10x_h5/` 中的 65 个过滤后 10x h5，与 10-14 髓系伪 bulk 的来源相同。

## 4. 流程

1. **构建伪 bulk**：逐个读取 65 个 h5，按 barcode 取出主分析细胞，把同一供体的原始 UMI 计数直接相加，得到 Ensembl 层面的计数矩阵。读取代码沿用 `AxisA_pseudobulk_disease_20261014/code/01_build_lung_pseudobulk.R`，只改细胞筛选条件。断言：计数为非负整数；细胞总数与细胞表一致。
2. **伪 bulk 检验**：edgeR v4 依次做 filterByExpr（含设计矩阵）、TMM 和 estimateDisp，再做 fry；同一对象经 voom 转换后做 CAMERA（`inter.gene.cor = NA`）。设计为 `~ disease + age_z + sex + smoke`，与髓系主分析相同。另报告每个成员的 glmQLF logFC，用于描述成员方向。
3. **供体 UCell 检验**：聚合方式与 `lung_rank_and_scores.R` 相同，先取文库中位数，再取供体中位数。用双侧 Wilcoxon 秩和检验（`exact = FALSE`），并用 `conf.int = TRUE` 得到 Hodges–Lehmann 位置差及 95% CI。
4. **参照行**：用同一脚本对全部髓系细胞重算一遍，报告值须与 S6 Table 一致。断言如下：
   - fry P 与 0.1151 的差 < 1e-6；
   - CAMERA P 与 0.05617 的差 < 1e-6；
   - Wilcoxon P 与 0.01032 的差 < 1e-6。

   髓系的 Hodges–Lehmann 区间是新增的报告项。

## 5. 检验家族与敏感性分析

- **主检验**：肺泡巨噬细胞 COPD vs 对照这一个假设。fry 与 CAMERA 是对同一假设的两种检验，与正文现有做法一致，两者都报告，不再校正。供体 UCell 检验作为与之对照的另一读出。
- **敏感性分析（只作描述）**：
  - a. 每名供体至少 50 个细胞（预计剩 54 名供体）；
  - b. 去除 2 个线粒体离群文库（GSM9102301、GSM9102263）；
  - c. 扩大到肺泡巨噬细胞家族，即 `Alveolar macrophages`、`Alveolar Mph CCL3+`、`Alveolar Mph proliferating`、`Alveolar Mph MT-positive`。
- **描述性读出（不检验）**：每名供体肺泡巨噬细胞占全部髓系细胞的比例，报告两组的中位数和四分位距。可行性检查时已看到这一比例（COPD 0.48，对照 0.48）。没有安装 β-二项回归的 R 包，所以这一项不作推断。

## 6. 预先写定的解读规则

1. 主检验中 fry 与 CAMERA 的 P 均 ≥ 0.05：结论为"肺泡巨噬细胞内部未检出 Axis A 随 COPD 的升高"，结果节的细胞内检验写入这一结果并报告方向。
2. fry 与 CAMERA 均 P < 0.05 且方向为上调：肺侧"细胞内未检出"不再成立，结果节、讨论组成段和摘要相应句都要改写，交作者决定。
3. 方向为下调且 P < 0.05：如实报告"肺泡巨噬细胞内偏低"，组织升高与细胞内结果方向相反，在结果中直接写出。
4. 伪 bulk 与 UCell 结论不一致：两者并列报告，并核查评分尺度、过滤后成员数和协变量的差别。
5. 无论哪种结果，都不据此修改 Axis A 基因集、细胞数阈值或协变量。

## 7. 主要风险

- 细胞注释用的是 CellTypist best match（不做 majority voting），肺泡巨噬细胞的边界取决于参考图谱。敏感性分析 c 部分对应这一风险。
- 统计功效有限：33 名 COPD 对 23 名对照，无法排除小幅变化。因此报告 Hodges–Lehmann 区间。
- 解离和捕获会影响细胞回收，回收比例不等于组织中的真实比例；描述性比例只用来理解髓系类别与肺泡巨噬细胞结果的差别。
- 设计写定前，可行性检查已看过每名供体的细胞数和肺泡巨噬细胞比例，但没有看过任何疾病检验结果。
