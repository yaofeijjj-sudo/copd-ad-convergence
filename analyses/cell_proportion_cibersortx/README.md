# CIBERSORTx 细胞组成调整（2026-10-19 运行）

依据：修改意见 2.3；设计见 [DESIGN.md](DESIGN.md)。**稿件未修改，状态为候选结果，待作者决定是否入稿。**

**样本**：按作者 10-19 决定，去掉 GSE48350 的重复芯片 GSM1176257（与对照 GSM350078 的 CEL 文件逐字节相同，见 `AxisA_member_direction_robustness_20261018/results/adopted_sample_exclusions_20261019.tsv`）。它仍在上传的混合文件中，只是不进入模型。去掉前的结果保存在 [before_GSM1176257_removal/](results/before_GSM1176257_removal/)，两者结论相同（AD 合并调整后 0.09 vs 0.10）。

主分析为版本 B：建参考前去掉 137 个 Axis A 基因，两个特征矩阵中 Axis A 基因数均为 0。版本 A 保留全部基因，其特征矩阵中肺含 67 个、脑含 61 个 Axis A 基因，作为敏感性分析。基础模型与已采用的复现模型完全相同；脚本逐一核对了 7 个独立队列的基础效应，与已采用值（GSE48350 为去掉重复芯片后的值）的差异 < 10⁻⁸。

## 结论（版本 B）

调整方式：在基础模型中加入目标髓系比例，肺为 Myeloid，脑为 Immune（小胶质细胞 + 血管周围巨噬细胞），标准化后作为一个协变量。效应单位为 SD，括号内为 95% CI。

| 队列 | 基础效应 | 加入髓系比例后 | 减弱 | 加入全部比例后 |
|---|---|---|---|---|
| 肺 GSE76925（发现） | 0.54（0.12–0.96） | 0.50（0.23–0.76） | 8% | 0.31（0.09–0.54） |
| 肺 GSE47460（独立） | 0.44（0.15–0.72） | **0.25（0.06–0.44），P = 0.012** | 43% | 0.10（−0.02 至 0.22） |
| AD 6 队列合并（DL） | 0.49（0.21–0.76） | **0.10（−0.01 至 0.21），P = 0.073** | 80% | 0.00（−0.09 至 0.09） |

- **肺发现队列的来源模块 red**（图 2d 模型）：0.69 → 加入 Myeloid 后 0.65（0.35–0.95）。10-15 的 dtangle 结果为 0.77，方向和结论一致。
- **AD 各队列**（[adjusted_models_by_cohort.tsv](results/adjusted_models_by_cohort.tsv)）：加入 Immune 比例后，6 个队列的效应都降到 0.19 以下，CI 都包含 0。GSE48350 由 1.04 降到 −0.07。
- **版本 A 与 B 结论相同**：GSE47460 为 0.22（0.04–0.41）；AD 合并为 0.09（−0.00 至 0.19）。

**估计比例与 Axis A 评分高度相关**（[proportion_case_control.tsv](results/proportion_case_control.tsv)）：
- 肺 Myeloid：r = 0.76（GSE76925）、0.74（GSE47460）。
- 脑 Immune：r = 0.73–0.96。

**髓系比例的病例–对照差异**（调整基础模型协变量）：
- 肺 Myeloid：GSE76925 +0.06 SD（P = 0.78）；GSE47460 +0.26 SD（P = 0.075）。
- 脑 Immune：GSE48350 +1.25 SD（P = 2.1 × 10⁻⁵）、GSE140841 +0.82 SD（P = 0.011）；其余 4 个队列 +0.02 至 +0.32 SD，P ≥ 0.09。

**判读**（按 DESIGN.md 预先写定的规则）：
- 系数减弱说明组织中 Axis A 与疾病的关联对组成调整敏感，不能据此计算中介比例。
- 肺独立队列在调整后仍有关联；AD 侧的关联在调整后不再能与 0 区分。
- 脑侧 Immune 比例与 Axis A 评分的相关高达 0.73–0.96，两者在整块组织数据中几乎不可分。因此，调整后效应接近 0，只说明整块组织数据无法把"免疫细胞变多"和"每个免疫细胞内 Axis A 升高"分开，不能证明细胞内没有变化。后一问题需要依靠单细胞或伪 bulk 结果。

## 拟合质量（[cibersortx_fit_quality.tsv](qa/cibersortx_fit_quality.tsv)）

- 16 次运行中，每个样本的置换检验 P 都 < 0.05（100 次置换）。
- 各队列观测与重建表达的中位相关为 0.64–0.91。最低的是 GSE118553（0.64）。

## 实际运行设置

- **网站**：CIBERSORTx 网页版（cibersortx.stanford.edu），2026-10-19 运行（服务器时间 10-18 22:40 起）。账号同时最多运行 5 个任务。
- **建特征矩阵**（4 个任务）：
  - scRNA-Seq 模式；关闭分位数归一化（网站默认）；kappa 999、q 值 0.01、每类 300–500 个标志基因；未过滤非造血基因。
  - Min. Expression（fraction）0.75、Replicates 5、Sampling 0.5，均为网站默认值。运行日志记录 `fraction: 0.75`；任务参数页显示为 "Min. Expression: 1"，以日志为准。
  - 网站说明建议：液滴数据若因稀疏导致入选基因过少，可把阈值降到 0.5。本次特征矩阵含 2,268–2,654 个基因，未触发这一条件，保持默认值。
- **估计比例**（16 个任务）：
  - Custom、相对模式、S 模式批次校正，提供建矩阵时使用的同一单细胞参考文件（网站在 S 模式下要求此文件）。
  - 芯片数据启用分位数归一化，RNA-seq（GSE313307、GSE140841）不启用。
- **文件核验**：
  - 上传文件超过工具的 10 MB 限制，先在本地切块，在浏览器中拼回。拼回后与本地文件的 SHA-256 一致，服务器上显示的文件大小也一致（[upload_sha256.tsv](qa/upload_sha256.tsv)）。
  - 结果从网站取回本地时逐个核对 SHA-256（[downloaded_outputs_sha256.tsv](qa/downloaded_outputs_sha256.tsv)）。
  - 特征矩阵见 `qa/signature_matrix_*.txt`。

## 与设计的差异

- **置换次数**：设计写的是"取网站默认值"（none）。网站说明建议 ≥ 100，本次设为 100。置换只影响拟合 P 值，不改变比例估计。
- **GSE118553 样本数**：混合文件含 62 个样本，模型使用其中 61 个，因为 1 个样本协变量缺失。这一处理与已采用的复现模型相同。

## 限制

- 单细胞（核）参考与芯片平台差异大，S 模式只能部分校正。
- 肺参考中结构细胞偏少；脑参考按 NeuN 分选，神经元比例偏高。估计值是相对比例，绝对值不作解读。
- 网页版有数据点上限，参考每类只抽 250–300 个细胞。
- 未做：用留出供体构建模拟混合样本检查可分辨性（修改意见 2.3）。

## 运行

```
LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/CIBERSORTx_composition_adjustment_20261018/code/05_models.R
```

从项目根目录运行。输入为 `results/raw/<版本>_<队列>.tsv`（CIBERSORTx 输出原文件）；代码 01–04 生成上传文件。R 运行时出现的 "stack imbalance" 警告不影响结果，基础效应与已采用值逐位一致。

## 文件

- [adjusted_models_by_cohort.tsv](results/adjusted_models_by_cohort.tsv)：每个队列三种模型。
- [adjusted_models_AD_pooled.tsv](results/adjusted_models_AD_pooled.tsv)：AD 合并。
- [GSE76925_red_module_adjusted.tsv](results/GSE76925_red_module_adjusted.tsv)：发现队列来源模块。
- [proportion_case_control.tsv](results/proportion_case_control.tsv)：比例差异及与评分的相关。
- [run_05.log](qa/run_05.log)：运行日志。
