# CIBERSORTx 细胞组成调整（2026-10-18）

依据：修改意见第 2.3 节；作者 10-18 决定使用 CIBERSORTx，肺用 GSE76925、GSE47460，脑用 6 个 AD 队列，调整前后分别合并。10-02 的设计（`Deconvolution_composition_adjustment_20261002/protocol/DESIGN.md`）已把 CIBERSORTx 列为肺侧方法，当时因没有账号未运行。

## 问题与估计量

- **问题**：组织中 Axis A 评分与疾病的关联，有多少可由估计的细胞组成解释？
- **估计量**：同一批样本中，两个模型的疾病系数（标准化 Axis A 评分，单位 SD）。
  - 基础模型：Axis A 评分 ~ 疾病 + 协变量（与已采用的复现模型相同）
  - 组成模型：在基础模型上加入估计的细胞比例
- **独立单位**：供体。

## 输入（冻结）

- **组织表达**：`AxisA_member_direction_robustness_20261018/source_data/cohorts.rds`，样本与协变量和已采用的复现完全相同（该分支已逐位重现）。GSE33000 是双色比值数据，不做 CIBERSORTx，保留 BRETIGEA。
- **肺参考**：GSE302339 质控后细胞，沿用图 4 的 CellTypist 大类：Myeloid、T/NK、B/Plasma、Epithelial、Endothelial、Stromal。Other/Progenitor 只有 14 个细胞，不纳入。
- **脑参考**：SEA-AD DFC（2026-06-22 版，`UMIs` 层原始计数）。按 SEA-AD Class/Subclass 合并为 8 类：Excitatory（Glutamatergic）、Inhibitory（GABAergic）、Astrocyte、Oligodendrocyte、OPC、Immune（小胶质细胞与血管周围巨噬细胞）、Endothelial、VLMC & Perivascular。Sutton 2022 *Nat Commun* 的结论是：主要大类的去卷积准确，细分亚型之间共线性强。

## 参考抽样（不使用疾病信息）

- 网站限制每个矩阵最多 5,000 万个数据点。每类细胞在供体之间均衡抽取：先在每名供体内随机抽，再按细胞类型补足到固定数目。肺每类 300 个细胞，脑每类 250 个细胞，随机种子 20261018。
- 只保留在抽样细胞中至少有一个计数、且有唯一官方基因名的基因。
- 数值为原始 UMI 计数，即线性尺度。CIBERSORTx 会自行归一化。

## 版本

- **版本 B（主分析）**：建参考前去掉 137 个 Axis A 基因，使它们不能进入特征矩阵。这样做是为了避免"用 Axis A 基因估计的髓系比例去调整 Axis A"，与 10-15 dtangle 和 BRETIGEA 的主分析（版本 B）一致。
- **版本 A（敏感性）**：保留全部基因。

## CIBERSORTx 设置（网页版，参照官方教程）

1. **建特征矩阵**：Create Signature Matrix → scRNA-Seq，参数取网站默认值（kappa 999、q 值 0.01、每类标志基因数取默认范围、replicates 和 sampling 取默认）。不启用非造血基因过滤。
2. **估计比例**：Impute Cell Fractions → Custom，相对模式。
   - 启用 S 模式批次校正（单细胞来源的特征矩阵），提供同一参考文件。
   - 芯片数据启用分位数归一化，RNA-seq 不启用（网站说明）。
   - 置换次数取网站默认值。
3. **组织数据**：芯片的 log2 值取 2^x 还原为线性尺度。RNA-seq 用 CPM（GSE313307 为原始计数，GSE140841 为 kallisto 估计计数）。基因用 Entrez 对应的当前官方基因名。

## 下游模型

- **主调整**：基础模型加入目标髓系比例（肺 Myeloid，脑 Immune），标准化后作为一个协变量。这与 dtangle 和 BRETIGEA 的主模型相同。
- **全组成调整（敏感性）**：加入除参考类以外的全部比例。参考类为平均比例最高的一类，因为比例之和为 1。
- **报告内容**：
  - 每个队列调整前后的疾病系数、95% CI 和人数；
  - 6 个 AD 队列分别用 DerSimonian–Laird 合并调整前和调整后的效应；
  - 髓系比例的病例–对照差异；
  - 髓系比例与 Axis A 评分的相关。
- **肺发现队列 GSE76925**：另报告来源模块 red 特征向量的同一调整（图 2d 模型），与 dtangle 结果直接对比。
- **解释**：系数减弱说明关联对组成调整敏感，不计算中介比例。系数不变也不能单独证明细胞内改变，需要结合单细胞结果。

## 风险

- 单细胞（核）参考与芯片平台差异大，S 模式只能部分校正。
- 肺参考中结构细胞偏少；脑参考按 NeuN 分选，神经元比例偏高。这些影响参考中的细胞数量，不影响每类细胞的表达谱。
- 网页版的数据点上限决定了每类抽取的细胞数。
- 结果与 dtangle、BRETIGEA 不一致时两者都报告，不按显著性取舍。
