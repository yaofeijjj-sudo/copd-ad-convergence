# 单细胞补充分析与 decoupleR 转录因子活性（修改意见 2.4、2.5、第 4 节；2026-10-19，运行前写定）

稿件未修改；结果为候选。独立单位：供体。本轮只使用已保存的供体伪 bulk 原始计数，不重新读取原始单细胞矩阵。

## 输入（冻结）

| 器官与类别 | 计数来源 | 供体 | 主要暴露 |
|---|---|---|---|
| 肺 Myeloid（全部髓系） | `AxisA_pseudobulk_disease_20261014/source_data/lung_myeloid_pseudobulk_counts.tsv.gz` | 56（33 COPD / 23 对照） | COPD vs 对照 |
| 脑 Microglia 2（Micro-PVM_2） | `虚拟敲除/AxisA_疾病标尺核对_20260925/source_data/brain_donor_pseudobulk_counts.tsv.gz` | 79 | 每 1 SD 的 CPS_Local |
| 脑小胶质总体（Micro-PVM_2、2_1、2_3、3、4） | `T06E_SEAAD_independent_state_evaluation_20260827_candidate/inputs/SEAAD_DFC_microglia_pseudobulk_counts.mtx.gz` | 79 | 每 1 SD 的 CPS_Local |

- 两个脑类别的细胞核筛选相同：10Xv3.1、非失败组分。
- **CPS_Local**：来自 SEA-AD 官方表 `sea-ad_global_and_local_cps_050126.xlsx`（sheet "Travaglini 2026"，`Brain Region` = DFC），是 DFC 局部的病理进展评分。
- **诊断（次要）**：SEA-AD 供体元数据中的 Overall AD neuropathological Change（ADNC）。Intermediate 或 High 记为 AD 神经病理改变，Not AD 或 Low 记为对照；与 CPS 分开报告，不视为同一暴露。
- **肺表型**：取自 GEO 样本信息台账（原作者分组）——`emphysema_group`（E-COPD 12 / NE-COPD 22 / 对照 24）和 `fev1`（FEV₁ % 预计值）。
- **CollecTRI**：从 OmniPath 下载（2026-10-19），按 decoupleR 2.16.0 `get_collectri(split_complexes = FALSE)` 的官方逻辑处理（`code/00_collectri.R`），SHA-256 已记录。

## 2.4 伪 bulk 稳健性

- **复现**：先重现 10-14 的 edgeR fry 主结果（肺 P = 0.115，脑 Microglia 2 P = 0.316）。
- **CAMERA**：同一 DGEList、filterByExpr 和 TMM，经 voom 后用 limma CAMERA 检验固定的 Axis A 集合，`inter.gene.cor = NA`，与 fry 并列报告。
- **协变量**：
  - 肺：年龄、性别、吸烟状态。
  - 脑：年龄、性别、PMI、RIN、研究。
- **细胞数阈值**：主阈值为每名供体 ≥ 20 个细胞；敏感性阈值事先固定为 ≥ 50 个细胞，报告因此损失的供体数。
- **多重校正**：器官内，对预先列出的主要类别（肺 1 个，脑 2 个）的主检验做 BH 校正。
- **不做**：肺泡巨噬细胞单独成类的伪 bulk。这需要重新读取 65 个原始 10x 矩阵，本轮不做，作为局限说明。

## 2.5 临床分层（肺 Myeloid，fry 与 CAMERA）

- **肺气肿分组**：设计为 3 组（对照、E-COPD、NE-COPD）+ 年龄 + 性别 + 吸烟，直接检验 E-COPD vs NE-COPD 的对比（即两组差异本身），同时报告两组各自与对照的比较。不以“一组显著、另一组不显著”作推断。
- **FEV₁**：连续变量，每 1 SD FEV₁ % 预计值，全部 56 名供体，同样调整年龄、性别、吸烟；不同时放入诊断变量（两者高度相关）。

## 第 4 节：decoupleR 转录因子活性

- **方法**：decoupleR 2.16.0，`run_ulm`（官方推荐），网络为 CollecTRI，`minsize = 5`。输入为每个类别的 log-CPM（filterByExpr + TMM）。
- **活性检验**：TF 活性矩阵按与上面相同的设计用 limma 检验主暴露。
  - **主家族**：CellOracle 前十名（STAT3、FOS、LYL1、JUNB、CEBPB、ELF1、RREB1、KLF6、KLF7、MEF2A），在每个类别内 BH 校正。
  - **次家族**：全部可评估的 TF，BH 校正。
- **与 Axis A 的关系**：同时报告每个 TF 的活性与供体 Axis A 评分（singscore，同一 log-CPM）的相关。
- **敏感性**：从网络中去掉 Axis A 靶基因后重算活性，以免共同基因造成机械相关。去掉后可评估靶基因 < 5 的 TF 标为不可评估。

## 判读

- 依据修改意见第六节。
- 阴性结果结合效应量与区间讨论：FDR 0.06 不写成显著激活，FDR 0.4 也不能证明没有效应。
- 不按 P 值调整阈值、供体或协变量。
