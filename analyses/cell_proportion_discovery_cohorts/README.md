# 去卷积组成调整（2026-10-15 运行，候选，待与作者讨论）

设计：[protocol/DESIGN.md](protocol/DESIGN.md)（10-02）及运行前补充 [protocol/DESIGN_RUN_20261015.md](protocol/DESIGN_RUN_20261015.md)，两者的 SHA-256 记录在 `protocol/DESIGN_RUN_20261015.sha256`。主分析为版本 B（去掉与 Axis A 重叠的标志基因）。稿件未修改。

## 主要结果（[adjusted_models.tsv](results/adjusted_models.tsv)）

| 模型 | 基线疾病 β（95% CI） | 加入比例后 β（95% CI） | 变化 | 判定 |
|---|---|---|---|---|
| **肺 red（Axis A 来源）+ 髓系比例**（dtangle，B） | 0.686（0.276–1.095） | **0.774（0.436–1.112）** | −12.9%（变大） | CI 排除 0，方向不变：关联不能由髓系比例解释 |
| **脑 blue（Axis A 来源）+ 小胶质比例**（BRETIGEA，B） | 1.108（0.927–1.288） | **0.602（0.455–0.750）** | 45.6% | CI 排除 0，方向不变：约一半随小胶质比例变化 |
| 肺 turquoise（Axis B 来源）+ 髓系比例（次要） | −0.713 | −0.720（−1.123 至 −0.317） | −1.0% | 不受影响 |
| 脑 turquoise（Axis B 来源）+ 神经元比例（次要） | −1.071 | **−0.067（−0.130 至 −0.004）** | 93.7% | 几乎全部随神经元比例变化 |

**敏感性分析：**
- 版本 A、C 与 B 结论相同：肺 red 加入髓系比例后 β 为 0.75–0.80；脑 blue 为 0.54–0.61。
- 多类型同时调整：
  - 肺 red：0.44–0.51，CI 仍排除 0。
  - 脑 blue：约 0（−0.05 至 −0.01），CI 包含 0。但脑侧六个 BRETIGEA 替代比例高度共线（神经元与 OPC r = −0.89，方差膨胀因子最高 10.7），合起来能解释 AD 状态 58% 的方差（[collinearity_versionB.json](qa/collinearity_versionB.json)）。这个模型实际上是在控制疾病严重度本身，属于 DESIGN.md 第 5 节所说的过度调整，不能解释为“组成解释了全部关联”。

## 估计比例在病例与对照间的差异（[proportion_case_control.tsv](results/proportion_case_control.tsv)，版本 B）

- **肺**（dtangle，调整年龄、性别、包年）：
  - 髓系比例**无差异**：病例减对照 −0.16 SD，P = 0.45；平均估计值对照 0.193、COPD 0.185。
  - **T/NK（+0.47 SD，P = 0.03）和 B/浆细胞（+0.49 SD，P = 0.02）在 COPD 中升高**，与原作者单细胞结果及 Hogg 2004 的淋巴滤泡增多方向一致。
- **脑**（BRETIGEA，调整年龄、性别）：
  - 神经元 −1.06 SD；小胶质 +0.89 SD；OPC、内皮、星形胶质均升高（P 均 < 1e-12）；少突胶质无差异。

## 与模块特征基因的相关（[ME_proportion_correlations.tsv](results/ME_proportion_correlations.tsv)，版本 B）

- 肺 red 与髓系比例 r = 0.53，与 T/NK 比例 r = 0.70，与 B/浆细胞比例 r = 0.48。
- 脑 blue 与小胶质比例 r = 0.82；脑 turquoise 与神经元比例 r = 0.98。

## 方法与核验

- **表达矩阵**：用原脚本的同一组函数从原始 series matrix 重建，平台注释用存档的提交者平台表（离线）。重建矩阵在原网络的 8,873 / 8,807 个特征上与保存的 `dat_expr` 逐值相同（最大差 0）。基线模型复现图 2d 的效应（最大差 4e-15）。
- **肺参考**：GSE302339 的 211,450 个属于图 5 六个大类的 QC 细胞，各大类取平均 CPM 后再取 log2；在 13,560 个共同基因（版本 B）上运行 dtangle 2.0.10，`data_type = "microarray-gene"`，其余参数为默认。
- **脑**：BRETIGEA 1.0.3 `findCells`，参数为默认（nMarker 50，SVD，scale）。版本 B 去掉了 14 个与 Axis A 重叠的小胶质标志基因。
- **实际使用的标志基因数**：[lung_dtangle_markers.tsv](results/lung_dtangle_markers.tsv)、[brain_BRETIGEA_markers.tsv](results/brain_BRETIGEA_markers.tsv)。
- **未执行**：CIBERSORTx（需要注册账号）；脑侧 dtangle（GSE33000 是双色芯片 log10 比值数据，超出 dtangle 的文档输入）。原因见运行前补充设计。

## 局限

- 估计的是相对比例，绝对值不作解读。肺参考来自解离后的单细胞，结构细胞不足；芯片数据与 10x 参考之间存在跨平台差异。
- 脑侧 BRETIGEA 给出的是替代比例（标志基因的第一奇异向量）。AD 中神经元丢失强烈，会同时牵动所有类型的估计。
- 细胞比例可能处在因果链中间（疾病 → 细胞组成改变 → 模块变化），所以调整后效应变小不等于原关联是假的。

## 文件

- 代码：`code/01_rebuild_bulk.R` → `02_lung_reference.R` → `03_deconvolve_and_model.R` → `04_collinearity_qa.R`，均从项目根目录运行，需设置 `LC_ALL=en_US.UTF-8`。
- 日志和会话信息：`logs/`。
