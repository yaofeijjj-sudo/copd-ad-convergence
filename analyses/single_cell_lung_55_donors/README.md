# 肺单细胞供体合并重跑（GSE302339 cood8 = copd8，2026-10-29）

**状态：已采纳**（作者 10-29 决定全面更新）。图 3、图 6、S2 File 的 Table F/L/M、S1 Data 和中英文正文已改用本目录的结果。

## 问题

GEO 样本表中受试者编号 `cood8`（GSM9102255，标题 s2，实验室编号 "COPD 8 LUL"）与 `copd8`（GSM9102261，s79，"COPD8 LLL"）的临床记录完全相同：68 岁男性，GOLD 4，FEV₁ 27，FEV₁/FVC 74，40 包年，曾吸烟，戒烟时间 "November"，肺气肿型 COPD。同一数据集中其余 8 名多肺叶供体（copd2、copd4、copd6、copd10、dnaz1、dnaz4、dnaz6、dnaz7）各文库共用一个编号，因此 `cood8` 是录入错误，两者是同一供体两个肺叶的文库。原流程按 `subject_id` 合并供体，把这名供体算作两人：56 名供体（33 COPD / 23 对照），实际为 55 名（32 / 23）。

## 设计

- 估计量、检验、协变量、基因集（136 个可测 Axis A 成员）、阈值（每名供体每类 ≥ 20 个细胞）和随机种子均与原分析相同，只把 `cood8` 映射为 `copd8`。
- 供体内合并沿用各原脚本的规则：定位分析为文库均值的等权平均；疾病检验为文库中位数再取中位数；伪 bulk 为计数相加；CellOracle 为供体内细胞均值，再在供体间等权平均。
- 每个脚本先用原编号复现冻结结果（关卡），再用合并后的编号重跑。

## 代码与结果

| 脚本 | 内容 | 复现的冻结来源 |
|---|---|---|
| `code/01_localization_and_donor_tests.R` | 肺六大类定位（CR2）与 5,000 个匹配随机集合；髓系亚群 CR2 与图 3e 供体点；六大类 Wilcoxon（BH）、肺气肿分型、HL 位置差；去除 2 个离群文库 | SYN01 `lung_rank_and_scores.R`、`model_and_null.R`、`lung_myeloid.R`；`Lung_sc_library_QC_check_20261028` |
| `code/02_pseudobulk_AM_TF.R` | 肺泡巨噬细胞与髓系伪 bulk fry/CAMERA、UCell Wilcoxon、协变量线性模型、AM 占比；肺髓系 5 个伪 bulk 模型与 decoupleR TF 活性（脑侧行原样复制） | `AxisA_AM_pseudobulk_20261028`；`SingleCell_pseudobulk_supplements_20261019` |
| `code/03_celloracle_lung_reaverage.py` | CellOracle 肺侧 delta_X 按 55 名供体重新等权平均（网络与模拟不变），357 个 TF 排序 | `CellOracle_AxisA_Entrez137_20261001/finalize.py`；`CellOracle_ranking_scale_20261019/01_ranking.py` |
| `code/04_localization_table.py` | 两器官定位表（肺 55 名 + 脑原值，Holm） | REF `localization.tsv` |

逐项新旧对照：`results/frozen_vs_merged_01.tsv`、`results/frozen_vs_merged_02.tsv`、`logs/03.log`。

关卡说明：六大类 Wilcoxon 的供体中位数与冻结表逐值一致；冻结表中 B/Plasma 的 P 值在第 5 位小数上不同（0.44216 vs 0.44218，冻结运行中秩次并列的数值差异），因此该关卡容差为 1e-4。其余关卡容差 ≤ 1e-8。

## 主要变化（冻结 56 名 → 合并 55 名）

- 肺髓系定位效应 0.1104 → 0.1104；5,000 个随机集合均无超过。
- 肺髓系 COPD vs 对照：P 0.0103 → 0.0107，FDR 0.062 → 0.064；HL −0.021（−0.042 至 −0.006）→ −0.022（−0.042 至 −0.005）。去除离群文库后：55 名 P = 0.0137 → 54 名 P = 0.0143。
- 肺泡巨噬细胞：Wilcoxon P 0.036 → 0.035；HL −0.023（−0.044 至 −0.001）→ −0.023（−0.045 至 −0.001）；伪 bulk fry P 0.25 → 0.18，CAMERA P 0.13 → 0.14。AM 占髓系比例中位数：COPD 0.48 → 0.49，对照 0.48。
- 髓系伪 bulk（调整年龄、性别、吸烟）：fry P 0.115 → 0.101，CAMERA P 0.056 → 0.052。
- CellOracle：357 个 TF 的共同下降、共同上升和方向相反的计数全部不变；KLF4、PPARG、MEF2C、BHLHE41 的器官内排名不变。
- decoupleR：STAT3、C/EBPβ 的 r 由 0.46–0.57 变为 0.45–0.57（去除 Axis A 靶基因后由 0.45–0.56 变为 0.44–0.56）。去除 Axis A 靶基因的先验下，肺髓系相关最强前 10 位的边界处 SPI1（r 0.648 → 0.636，第 10 → 11 位）与 TEAD1（r 0.646 → 0.640，第 11 → 10 位）交换，"三个细胞类别均居前 10 位"的 TF 由 STAT1、SPI1 变为 STAT1、TEAD1。正文中 PU.1 相关表述的处理由作者决定。
- 肺六大类中，Stromal 的合格供体由 29 名变为 30 名，T/NK 与 Stromal 的 P 降至 0.043（FDR 0.086）；其余各类无实质变化。

## 未重跑的部分

- CIBERSORTx 肺参考（网页版，手工构建）：1,800 个细胞覆盖 55 名真实供体，其中这名供体贡献 65 个细胞，约为平均的 2 倍；参考仍是真实细胞谱，未重建。正文与 Table F 的覆盖数改为 55。
- 单细胞注释、质控、UCell 评分与 UMAP 与供体编号无关，不受影响。
- `05_github_release/` 中的参考表（localization、lung_null、lung_subtype、CellOracle 等）仍为 56 名供体版本，上传 GitHub 前需同步。
