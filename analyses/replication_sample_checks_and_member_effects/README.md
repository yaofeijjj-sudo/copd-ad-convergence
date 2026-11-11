# Axis A 成员方向、复现稳健性与供体核查（2026-10-18）

依据：购买的修改意见第 2.1、2.2 节；设计见 [DESIGN.md](DESIGN.md)。**稿件未修改，状态为候选结果，待作者决定是否入稿。**

已采用的复现结果（`AxisA_singscore_official_rerun_20261004`）只读未改。本分支先完整重现这些结果：7 个队列的评分效应与已采用值的差异均 < 10⁻¹⁵，可测成员数完全一致（[qa/check_reproduce_adopted_effects.tsv](qa/check_reproduce_adopted_effects.tsv)）。

## 结论

1. **成员方向（意见 2.1）**
   - 137 个成员中，108 个在独立肺队列（GSE47460）与 6 个 AD 队列合并后都升高，5 个都降低，24 个两侧方向相反。
   - 两侧 95% CI 都在 0 以上的有 34 个。
   - 发现队列中，134 个可比成员里 113 个两侧都升高。
   - 全部成员与来源模块特征向量正相关（kME 中位数：肺 0.63，脑 0.67）。两个来源模块特征向量都与疾病正相关，方向约定一致。
   - 方向相反的 24 个多数是"一侧接近 0"。差异较大的有：PLA2G7、MR1、PTPN7、CLEC10A、FCMR（肺升、脑降），以及 IL16、SLC7A7、SLC15A3、HLA-DOA（肺降、脑升）。
   - 结果见 [member_effects_wide.tsv](results/member_effects_wide.tsv)、[member_direction_counts.tsv](results/member_direction_counts.tsv)，图见 [heatmap](figures/member_direction_heatmap.pdf)、[scatter](figures/member_lung_brain_scatter.pdf)。
2. **逐一剔除成员**：评分效应不依赖少数基因。GSE47460 为 0.426–0.443，AD 合并为 0.446–0.466，所有剔除情形 P ≤ 2.1 × 10⁻⁴（[leave_one_member_out.tsv](results/leave_one_member_out.tsv)）。
3. **复现稳健性（意见 2.2）**，AD 合并效应（SD，95% CI，P）：

   | 分析 | 效应 | 95% CI | P |
   |---|---|---|---|
   | 主分析（DerSimonian–Laird） | 0.46 | 0.22–0.69 | 1.5 × 10⁻⁴ |
   | REML | 0.46 | 0.22–0.69 | 1.3 × 10⁻⁴ |
   | REML + Hartung–Knapp | 0.46 | 0.15–0.77 | 0.013 |
   | 只用 108 个共同可测成员 | 0.46 | 0.21–0.71 | 2.9 × 10⁻⁴ |
   | GSE84422 去掉 possible AD | 0.47 | 0.25–0.69 | 2.7 × 10⁻⁵ |
   | GSE313307 只用非感染供体 | 0.51 | 0.28–0.74 | 1.1 × 10⁻⁵ |
   | 不含 GSE313307 | 0.53 | 0.22–0.84 | 8.1 × 10⁻⁴ |

   GSE313307 的分层结果（[GSE313307_strata.tsv](results/GSE313307_strata.tsv)）：
   - 四组人数为 AD/对照 × 无/有感染 = 55/52 和 58/37。
   - 无感染者中 AD 高 0.45 SD（P = 0.021），有感染者中高 0.12 SD（P = 0.61）。
   - AD × 感染交互项不显著（−0.24，P = 0.41）。
4. **供体与样本核查**
   - **GSE48350 存在重复芯片。** 标为 AD（女，86 岁）的 GSM1176257，其 CEL 文件解压后与标为对照（男，80 岁，来自 GSE11882）的 GSM350078 逐字节相同（MD5 f8ca2cef…，扫描时间同为 2006-03-30 17:02:23）。性别基因显示这张芯片来自男性，即属于 GSM350078。这是 GEO 提交错误，原分析把同一张芯片同时计为病例和对照。去掉 GSM1176257 后，GSE48350 效应为 1.04 SD，AD 合并为 0.49 SD（P = 5.1 × 10⁻⁴）。
   - **性别基因与标注不一致。** 以 Y 连锁基因减 XIST 判断性别：GSE33000 有 7/459 例，GSE47460 5/236，GSE118553 6/62，GSE313307 1/202，GSE48350 4/41（含上述重复芯片），GSE76925、GSE140841、GSE36980 为 0。GSE84422（GPL96）的性别基因信号不呈双峰，无法判断。去掉不一致样本后，GSE47460 为 0.46 SD，AD 合并为 0.49 SD（P = 6.6 × 10⁻⁴）。
   - **跨队列重复供体：未检出。** doppelgangR 只在跨数据集比较中标出一个"枢纽"样本 GSE84422 GSM2234299：它与各数据集多名年龄、性别都不同的样本高度相关，属于离群样本，不是重复。GSE118553 与 GSE140841（均由 King's College London 提交）之间没有表达重复。肺侧有 4 对元数据完全相同（年龄、性别、FEV₁）：一对诊断不同，其余三对的表达相关在 236 个候选中排第 72–231 位，不是同一供体。
   - **GSE47460 内部有两对对照表达近乎相同**（LT058156/LT075094，r = 0.997；LT115840/LT079487，r = 0.991），无法确定是否为同一人。去掉这 4 人后效应为 0.42 SD（P = 0.0045）。
   - 详见 [duplicate_impact.tsv](results/duplicate_impact.tsv)、[sex_check_mismatches.tsv](results/sex_check_mismatches.tsv)、[duplicate_doppelgangR_summary.tsv](results/duplicate_doppelgangR_summary.tsv)。
5. **队列核查表**：[cohort_audit.tsv](results/cohort_audit.tsv) 列出来源、脑区、诊断、人数、协变量、可测成员数、表达匹配随机基因集经验 P 和 CAMERA P。不可测成员逐个列出原因（[member_coverage_summary.tsv](results/member_coverage_summary.tsv)）：
   - RNA-seq 队列主要因低表达过滤：GSE313307 去掉 17 个，GSE140841 去掉 9 个；另各有 1 个（LILRB3）无唯一 Ensembl→Entrez 映射。
   - 芯片队列因平台没有对应探针（例如 GPL96 和 GPL570 都没有 CARD9、CD68）。

## 运行

```
R_LIBS=_workspace/r_libs_affy:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/01_prepare_cohorts.R
R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/02_member_effects_and_sensitivity.R
R_LIBS=_workspace/r_libs_doppel:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/03_duplicate_donors.R
R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/04_figures.R
R_LIBS=_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/AxisA_member_direction_robustness_20261018/code/05_followup_checks.R
```

从项目根目录运行。doppelgangR 1.38.0 及其依赖（sva 3.58.0、genefilter、mgcv、mnormt）安装在 `_workspace/r_libs_doppel`；编译时需要 `SDKROOT=$(xcrun --show-sdk-path)`、把 `_workspace/envs/r/bin` 加入 PATH，并用 `Makevars_noomp` 关闭 OpenMP。

## 限制与说明

- 发现队列的成员效应有循环性，只作描述。
- 137 个成员相互相关，方向计数不作二项检验。
- doppelgangR 在跨平台比较中检出能力有限，"未检出"不能完全排除重复。
- 性别不一致可能来自元数据错误或样本调换，本分支不作区分。
- GSE76925 的来源（Temple/Brigham）取自 GEO 贡献者名单，未与原论文核对。
- R 运行时偶尔出现 "stack imbalance" 警告。它不影响结果：评分效应与已采用值逐位一致。
- 未做：单细胞供体流转（意见 5.1，56/80/79 名供体的差异）。
