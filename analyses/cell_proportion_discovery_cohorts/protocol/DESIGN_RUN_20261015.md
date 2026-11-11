# 去卷积组成调整：运行前确认与实施细节（2026-10-15）

作者 2026-10-15 授权运行（“进行去卷积分析。分析后的结果我们讨论之后再进行下一步”）。本文件在运行前写定，补充并具体化 [DESIGN.md](DESIGN.md)。两者冲突时以本文件为准，冲突之处逐条列在下面。结果只交给作者讨论，不改稿件。

## 运行时已知的信息（不是盲法）

- 同日伪 bulk 分析：肺髓系细胞内 Axis A 偏低（P = 0.12），脑 Microglia 2 内无变化（P = 0.32）。
- 原作者的单细胞比例分析（scCODA）：肺气肿型和终末期 COPD 的巨噬细胞在免疫细胞中的相对比例降低；SEA-AD 中 Micro-PVM_3 随 AD 增多。
- 组织学文献：COPD 肺组织中巨噬细胞绝对数增加（Retamales 2001；Hogg 2004）。

## 确认的选择

1. **标志基因版本**：主分析为版本 B（去掉与 Axis A 137 个 Entrez ID 重叠的标志基因）；A（全部）、C（去掉整个来源模块：肺 red 487 个、脑 blue 1,651 个已确认 Entrez ID）作为敏感性分析。
2. **CIBERSORTx 不做**：需要注册账号，本次不创建，按原设计“拿不到则不做”记录。
3. **脑侧 dtangle 敏感性不做（偏离原设计）**：GSE33000 是双色芯片相对参考池的 log10 比值（数值约 −1 到 1.2）。dtangle 的文档输入是 log 尺度表达强度（官方示例为 log2 RMA），比值数据超出其文档用法。BRETIGEA 只对标志基因做中心化后的 SVD，不依赖绝对尺度，作为脑侧唯一方法。
4. **脑侧不用 SEA-AD 参考**：因为第 3 条，不需要构建脑单核参考。

## 表达矩阵重建

- 原始文件：`02_data/raw/bulk_lung_copd/GSE76925_series_matrix.txt.gz`、`02_data/raw/bulk_brain_ad/GSE33000_series_matrix.txt.gz`。
- 预处理：原样复制原脚本 `15_final_wgcna_copd_ad.R` 中的 `normalize_microarray_matrix`、`collapse_probes_to_genes`、`first_token`，以及脑侧排除亨廷顿病、按组 1.5 × IQR 剔除年龄离群者的步骤。样本顺序和集合与原脚本一致；只是跳过 MAD 过滤（保留全部基因）。最终样本限定为 `sample_table_used.csv` 中的 147 / 459 人。
- 核验：重建矩阵在原网络 8,873 / 8,807 个特征上应与保存的 `dat_expr` 逐值相等（最大绝对差 < 1e-8）。不相等则停止并报告。
- 基因 ID：脑侧行名本来就是平台 EntrezGeneID。肺侧行名是平台 Symbol，用 org.Hs.eg.db 3.22.0 按官方符号唯一映射到 Entrez（不用别名，一对多的丢弃），与 `code/00_count_marker_overlap.R` 的规则相同。BRETIGEA 标志基因的符号用同一规则映射。

## 细胞比例估计

**肺：dtangle 2.0.10**（Hunt 2019）
- 参考：GSE302339 的 211,464 个 QC 细胞，按图 5 的 6 个大类分组。每个细胞先换算为 CPM，再在大类内取均值，最后取 log2(均值 + 1)。基因从 Ensembl 经 org.Hs.eg.db 唯一映射到 Entrez。
- 调用：`dtangle(Y = 肺 log2 矩阵（样本 × 共同基因）, references = 6 × 共同基因参考, data_type = "microarray-gene")`，其余参数取默认值（`marker_method = "ratio"`，`n_markers` 默认）。
- 版本 B / C：在调用前从共同基因中删去相应的 Entrez ID，由 dtangle 在剩余基因中自行选择标志基因。
- 主要协变量：`Myeloid` 估计比例。

**脑：BRETIGEA 1.0.3**（McKenzie 2018；Johnson 2021 用于 AD 芯片数据）
- 调用：`findCells(inputMat = 脑矩阵（Entrez 转为唯一官方符号）, markers = markers_df_human_brain, nMarker = 50, method = "SVD", scale = TRUE)`，即包的默认值。
- 版本 B / C：调用前从标志基因表中删去相应 Entrez ID 对应的标志基因，函数按原排序顺延取足 50 个。
- 主要协变量：`mic`（小胶质细胞）替代比例。

## 模型（与图 2d 相同，只增加比例项）

- 基线，须先复现图 2d（肺 red β = 0.6855，脑 blue β = 1.1079，误差 < 1e-6）：
  - 肺：`scale(ME) ~ copd_status + age_z + sex_male + packyears_z`
  - 脑：`scale(ME) ~ Disease_State + age_z + sex_male`
- **主要**：肺 red 加 `scale(Myeloid)`（dtangle，版本 B）；脑 blue 加 `scale(mic)`（BRETIGEA，版本 B）。
- **次要**：肺 turquoise 加同一髓系比例；脑 turquoise 加 `scale(neu)`（神经元）。
- **敏感性**：
  - 版本 A、C；
  - 多细胞类型同时调整：dtangle 纳入除参考类型以外的全部类型，参考类型定为平均估计比例最高的一类；BRETIGEA 纳入全部 6 个替代比例。

## 报告（结果不论好坏都报告）

每个模型报告：
- 疾病 β、95% CI、P（调整前后）；
- 变化百分比 = 1 − β调整后 / β基线；
- ME 与估计比例的 Pearson 相关；
- 估计比例在病例与对照间的差异：`scale(比例) ~ 疾病 + 同一协变量`；
- 实际使用的标志基因数。

判定规则沿用 DESIGN.md 第 4.4 节：只看主要模型调整后 β 的 95% CI 是否排除 0、方向是否不变。

## 风险（新增）

- 肺：单细胞参考是解离后的细胞，结构细胞严重不足（上皮约 7%）。参考谱本身可靠，但 dtangle 估计的是相对比例，绝对值不作解读。
- 版本 B 在肺侧仍会保留 Axis A 以外、来自同一髓系程序的标志基因，所以估计的髓系比例与 red 模块天然相关。这是“比例在因果链中间”的风险，不是循环；解释时按 DESIGN.md 第 1 节的边界处理。
