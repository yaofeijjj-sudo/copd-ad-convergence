# 分析设计：发现队列中细胞组成对模块–疾病关联的影响

状态：**暂缓（2026-10-02 作者决定）**。投稿前不运行，留待审稿人提出细胞组成问题时使用。已完成的只有标志基因重叠计数（`qa/00_marker_overlap_counts.txt`）；第 4.2 节主分析版本（A 或 B）尚未决定。运行前须重新确认并冻结本文件。
日期：2026-10-02。对应[审稿意见清单](../../../paper/审稿意见与修改清单_20261002.md)第 1 项。

## 1. 要回答的问题

Axis A 两个来源模块（肺 red、脑 blue）与疾病的调整后关联，有多少可以由组织中髓系细胞（肺）或小胶质细胞（脑）比例的差异解释？

**估计目标：**
- 主要：在现有 Figure 2d 模型中加入估计的细胞比例后，疾病效应 β 的变化。β 以模块特征向量的标准差为单位。
- 辅助：模块特征向量与估计比例的 Pearson 相关，沿用 Johnson 2021 的报告方式。

**解释边界：**
- 细胞比例可能处在因果链中间：疾病 → 髓系细胞增多 → 模块升高。所以调整后关联减弱，说明关联主要通过细胞组成体现，而不是说明原关联是“假的”。
- 估计比例有误差，调整不完全，残余的关联也不能全部归为细胞内表达变化。

## 2. 统计单位与样本

独立单位是参与者，样本完全沿用 Figure 2d：
- 肺 GSE76925：147 人（COPD 与对照）；
- 脑 GSE33000：459 人（306 AD，153 对照）。

不增减样本。

## 3. 冻结输入

| 输入 | 来源 |
|---|---|
| 模块特征向量与样本协变量 | 原 WGCNA 结果：`WGCNA_全部图件/复现资料/analysis/01_lung_wgcna/results/COPD_WGCNA/`（`module_eigengenes.csv`、`sample_table_used.csv`）；脑侧 `01_Results_and_Methods/01_Figure1_WGCNA/analysis/02_brain_wgcna/results/AD_WGCNA_no_residual_after_age_outlier/` 中的同名文件 |
| 现有调整后效应（对照基准） | `cross WGCNA/forest_revision_20260917/source_data/` 中的四行源表 |
| 组织表达矩阵 | 用原脚本 `15_final_wgcna_copd_ad.R` 中的原预处理函数（quantile 归一化、方差最大探针合并），重建 MAD 过滤前的全基因矩阵；探针身份按 Entrez 校正台账转为 Entrez ID。不重建网络，不改动模块 |
| 肺单细胞参考 | GSE302339，211,464 个质控后细胞，沿用 Figure 5 的 CellTypist 注释与六大类划分 |
| 脑单核参考 | SEA-AD DFC final-nuclei（2026-06-22 发布），官方 `Subclass` 注释 |
| BRETIGEA 标志基因 | 包内自带的人类脑细胞共识标志基因 |
| Axis A 成员 | `AxisA_confirmed_genes.csv`（137 个基因，SHA-256 `dd81a1af…6188`） |

所有输入在运行开始时记录 SHA-256；基因一律按 Entrez ID 对应，不用符号或别名。

## 4. 方法（均为已发表方法）

### 4.1 细胞比例估计

| 器官 | 主方法 | 敏感性方法 |
|---|---|---|
| 脑 GSE33000 | **BRETIGEA** `findCells`，用包内标志基因估计小胶质细胞等 6 类相对比例（McKenzie 2018 *Sci Rep*；Johnson 2021 *Sci Rep* 在 AD 芯片数据上的用法） | **dtangle**，标志基因由 SEA-AD DFC 单核参考按 dtangle 官方 `find_markers` 选取（Hunt 2019 *Bioinformatics*；Sutton 2022 *Nat Commun* 推荐用于脑组织） |
| 肺 GSE76925 | **dtangle**，标志基因由 GSE302339 单细胞参考按 `find_markers` 选取（Hunt 2019，支持芯片数据与跨平台参考） | **CIBERSORTx** absolute 模式，同一单细胞参考（Newman 2019 *Nat Biotechnol*；*BMC Genomics* 2024 LTRC 研究在 COPD 肺组织中的用法）。需要注册账号；拿不到则不做，并如实记录 |

**参考数据的构建：** 按细胞类型取归一化表达的均值。所有方法的参数都用官方默认值，不按结果调参。

**MuSiC 不用：** 它要求线性尺度的计数，芯片 log 强度数据不符合其前提（Avila Cobos 2020 *Nat Commun*）。

### 4.2 标志基因与 Axis A 重叠的处理（待作者确认）

用同一批基因估计比例、再去调整含这些基因的模块，存在循环问题。我没有找到专门讨论这一点的已发表文献，因此预先规定报告三个版本：

| 版本 | 标志基因 | 角色 |
|---|---|---|
| A | 方法默认的全部标志基因 | 敏感性 |
| B | 去掉与 Axis A 137 个基因重叠的标志基因 | **主分析（建议）** |
| C | 去掉与整个来源模块（肺 red 487 个、脑 blue 1,651 个已确认 Entrez 基因；240 和 853 是限于两侧共享背景的数目）重叠的标志基因 | 敏感性；若剩余标志基因太少导致无法估计，如实报告 |

### 4.3 调整模型

在 Figure 2d 现有模型中只加入估计比例（标准化后），其余不变：

- 肺：`ME_SD ~ COPD + age_z + sex_male + packyears_z + 髓系比例_z`
- 脑：`ME_SD ~ AD + age_z + sex_male + 小胶质比例_z`

**主要模块：** 肺 red、脑 blue（Axis A 来源模块）。
**次要模块：** 两侧 turquoise（Axis B 来源模块）。肺调整同样的髓系比例；脑另以神经元比例做一次，用来检验氧化磷酸化下降是否反映神经元丢失。

**敏感性：** 同时加入多类细胞的比例，并去掉一类作为参照，避免各类比例之和为 1 造成的共线性。这是 Patrick 2020 *PLoS Comput Biol* 的做法。

### 4.4 报告内容（预先规定，结果不论好坏都报告）

每个模块、每个方法、每个标志基因版本，都报告：
- 加入比例前后的 β、95% CI、P；
- β 的变化百分比 = 1 − β调整后 / β原始；
- 模块特征向量与估计比例的相关；
- 估计比例在病例与对照间的差异。

**判定规则：** 只看调整后 β 的 95% CI 是否排除 0、方向是否不变。不另设百分比阈值。变化百分比只作描述。

**结论的写法：**
- 调整后 CI 仍排除 0：写“关联在调整细胞比例后仍然存在，但幅度减弱 X%”。
- CI 包含 0：写“关联主要通过细胞组成体现”，正文核心表述相应改为“髓系细胞程序的共同增多”。
- 主方法和敏感性方法结论不一致时，如实报告不一致，以主方法 + 版本 B 为主，不挑选结果较好的那一个。

### 4.5 质量检查

- 各细胞类型实际使用的标志基因数；
- 不同方法估计的同一类细胞比例之间的相关（例如 BRETIGEA 与 dtangle 的小胶质比例）；
- 估计比例的分布是否合理（是否有大量 0 或 1）；
- 模型残差诊断。

结果中没有真实细胞计数作金标准，估计比例的绝对值不作解读。

## 5. 主要风险

1. **跨平台差异：** 参考是 10x 单细胞或单核数据，组织数据是芯片，平台差异会降低估计精度。Sutton 2022 的结论是参考数据的影响比算法大，所以脑侧同时报告不依赖单核参考的 BRETIGEA。
2. **过度调整：** 如果细胞比例处在因果链中间，调整会移除真实的疾病效应。这一点在结果解释中明确写出。
3. **调整不足：** 估计误差会使调整不完全。
4. **肺样本量：** 147 人，调整后的置信区间会变宽。CI 包含 0 要写成“不能排除”，不能写成“无关联”。
5. **资源冲突：** SEA-AD 单核文件 40 GB，构建参考时内存占用大，必须等 CellOracle 脑侧跑完再运行，不并行。

## 6. 不在本分析范围内

- 不重建 WGCNA，不改变模块、Axis 成员或其他已有结果；
- 不修改稿件，结果出来后由作者决定如何写入；
- 外周血的白细胞组成调整（清单第 5 项）另做设计。

## 7. 参考文献

- McKenzie AT et al. Brain cell type specific gene expression and co-expression network architectures. *Sci Rep* 2018. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5995803/>
- Johnson TS et al. Combinatorial analyses reveal cellular composition changes have different impacts on transcriptomic changes of cell type specific genes in Alzheimer's Disease. *Sci Rep* 2021. <https://doi.org/10.1038/s41598-020-79740-x>
- Hunt GJ et al. dtangle: accurate and robust cell type deconvolution. *Bioinformatics* 2019. <https://doi.org/10.1093/bioinformatics/bty926>
- Sutton et al. Comprehensive evaluation of deconvolution methods for human brain gene expression. *Nat Commun* 2022. <https://pmc.ncbi.nlm.nih.gov/articles/PMC8924248/>
- Patrick E et al. Deconvolving the contributions of cell-type heterogeneity on cortical gene expression. *PLoS Comput Biol* 2020. <https://pmc.ncbi.nlm.nih.gov/articles/PMC7451979/>
- Newman AM et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nat Biotechnol* 2019. <https://doi.org/10.1038/s41587-019-0114-2>
- Computational deconvolution of cell type-specific gene expression in COPD and IPF lungs reveals disease severity associations. *BMC Genomics* 2024. <https://pmc.ncbi.nlm.nih.gov/articles/PMC11654147/>
- Avila Cobos F et al. Benchmarking of cell type deconvolution pipelines for transcriptomics data. *Nat Commun* 2020. <https://pmc.ncbi.nlm.nih.gov/articles/PMC7648640/>
