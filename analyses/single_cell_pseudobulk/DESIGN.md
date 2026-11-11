# Axis A 供体伪 bulk 疾病比较：执行前设计（2026-10-14）

状态：候选分析，作者 2026-10-14 授权（“先做伪 bulk”）。不替换图 5 / S3 Table 已接受的 UCell 结果；结果无论方向和显著性都报告。

## 问题

在图 5 定位到 Axis A 的细胞类别内，固定的 Axis A 成员作为一个整体，是否在供体层面随疾病改变表达？这是对图 5 疾病检验读出方式的补充：图 5 读的是“每名供体 UCell 分数的中位数”，本分析读的是“每名供体原始计数加和后的基因集表达”。

## 不是盲法

已知信息：图 5 肺髓系 UCell 为 COPD 偏低（FDR = 0.062），脑 Microglia 2 与 CPS 无关（FDR = 0.42）；2026-09-25 诊断性核对（肺泡巨噬细胞、旧 135 基因、脑 ACT 内 High vs Not AD，limma-voom + fry）两侧方向均为 Down，P = 0.178 / 0.729。因此本分析不称为独立验证或预注册，只是按标准读出方式重新报告同一问题。

## 方法依据（作者官方）

- 伪 bulk：同一供体、同一细胞类别的原始 UMI 计数直接求和（Crowell 2020 *Nat Commun* muscat，DOI 10.1038/s41467-020-19894-4；Squair 2021 *Nat Commun*，DOI 10.1038/s41467-021-25960-2）。
- 检验：edgeR v4 官方流程（Chen 等 2025 *Nucleic Acids Res*，DOI 10.1093/nar/gkaf018）。原文写明单细胞数据可汇总为伪 bulk 后进行 edgeR 分析；对“来自先前独立研究”的预定义基因集，用 roast/mroast/fry 做自包含检验，且“在 edgeR 中 fry 总是优于 roast 或 mroast”。本分析按此使用 `edgeR::fry`（DGEList 方法，z-score 转换）。
- 版本：R 4.5.3，edgeR 4.8.2，limma 3.66.0。

## 固定设计

| 项目 | 肺 | 脑 |
|---|---|---|
| 数据 | GSE302339，65 个 GEO 原始 10x 矩阵 | SEA-AD DFC，10Xv3.1，通过原始 QC |
| 细胞类别 | CellTypist 大类 `Myeloid`（与图 5 一致） | Supertype `Micro-PVM_2`（Microglia 2） |
| 细胞集 | SYN01 已评分、通过 QC 的细胞（116,558 个） | 与图 5 UCell 完全相同的 23,198 个细胞核（复用 2026-09-25 计数表，SHA-256 b997c1f9…） |
| 独立单位 | 供体（`subject_id`）；同一供体多个肺叶文库求和 | 供体（`Donor ID`） |
| 纳入 | 每名供体该类别 ≥ 20 个细胞（图 5 规则）：COPD 33，对照 23 | 同规则：79 名供体 |
| 暴露 | COPD vs 对照（病例减对照） | 标准化局部 CPS（每 1 SD） |
| 主要模型 | `~ disease + age + sex + smoke`（与 GSE47460 批量复现的协变量一致：年龄、性别、吸烟状态） | `~ CPS_Local_z + age_z + sex + PMI_z + RIN_z + study`（与图 5 一致） |
| 次要模型 | `~ disease`（与图 5 未调整的 Wilcoxon 对应） | 无 |
| 基因集 | 当前 Axis A 137 基因（Entrez）在该数据中可测的 136 个 Ensembl ID（SYN01 `lung_axis_measured.txt`） | 同上（`brain_axis_measured.txt`，136 个；按 h5ad `gene_ids` 把基因符号行映射回 Ensembl） |

流程（两器官相同）：`DGEList` → `filterByExpr(y, design)` → `normLibSizes`（TMM）→ `fry(y, index = Axis A 中通过过滤的成员, design, contrast = 暴露系数)`。

> 执行记录（2026-10-14，看结果前）：`fry` 的 DGEList 方法要求先估计离散度，首次运行因此报错。按 edgeR 官方要求在 `fry` 前加入 `estimateDisp(y, design)`，检验和其他设定不变。

- **主要读出**：fry 的方向（Up/Down）和双侧方向性 P 值。每个器官一个主要检验；按稿件既有做法，多重校正家族在器官内，因此不再跨器官校正（同时列出两器官 BH 值供参考）。
- **描述性读出**：`estimateDisp` + `glmQLFit` + `glmQLFTest` 给出每个成员的 log2FC；报告上调成员比例和成员 log2FC 中位数。不对单基因作结论；单基因名义显著不能推翻整组检验（项目工作流规则）。

## 解释规则（看结果前写定）

- 方向为 Up 且 P < 0.05：髓系细胞内部 Axis A 随疾病升高，可补充“组织升高至少部分来自细胞内表达”。
- 方向为 Down 且 P < 0.05：细胞内部随疾病降低，与组织层面升高方向相反；如实写出。
- P ≥ 0.05：在该样本量下没有细胞内改变的证据；与图 5 一致，组织层面升高不能由细胞内表达解释。不写成“证明无差异”。
- 本分析不检验细胞组成；“组织升高来自细胞比例增加”仍需单独的组成分析（propeller）才能作为结论。

## 不做的事

不按结果更换细胞类别（如改用肺泡巨噬细胞）、基因子集、协变量、过滤阈值或检验函数；不新增其他敏感性分析。

## 主要风险

- 肺：COPD 组更年长、男性比例更高、包年更多；吸烟状态与诊断部分重叠（对照 9 名从不吸烟，COPD 仅 1 名）。
- 肺 `Myeloid` 是混合类别（肺泡巨噬、间质巨噬、单核细胞等）：其内部亚群比例的变化也会表现为伪 bulk 表达变化。这与图 5 检验的层级相同，但解释时需注明。
- 脑：单细胞核测序可能漏检部分小胶质细胞激活基因；脑侧计数沿用 CellOracle 输入对象中同一批细胞核的原始 UMI。
