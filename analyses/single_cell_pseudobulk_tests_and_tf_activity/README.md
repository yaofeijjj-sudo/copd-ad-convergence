# 单细胞伪 bulk 补充分析与 decoupleR 转录因子活性（2026-10-19）

依据：修改意见 2.4、2.5、第 4 节；设计见 [DESIGN.md](DESIGN.md)，运行前写定（SHA-256 见 `DESIGN.sha256`）。**稿件未修改，候选结果。**

本轮只用已保存的供体伪 bulk 原始计数，独立单位为供体。脚本先重现了 10-14 的 edgeR fry 主结果（肺 P = 0.115，脑 Microglia 2 P = 0.316，差异 < 10⁻¹⁰，脚本中有断言）。

## 结论

1. **细胞内的 Axis A 没有随疾病升高**，换方法、换类别、换阈值都一样。肺髓系细胞中方向反而略偏低。这再次支持“组织升高来自组成，而非细胞内改变”。
2. **肺气肿型与非肺气肿型 COPD 没有可靠差异**（直接比较 P = 0.07–0.18）；FEV₁ 与 Axis A 无关。
3. **CellOracle 前十名 TF 的活性都不随疾病或病理改变**（主家族 FDR ≥ 0.22）。但在供体之间，**STAT3 和 CEBPB 的活性与 Axis A 评分明显同步**（r = 0.45–0.57；全部 TF 的中位数为 0.14–0.18），从网络中去掉 Axis A 靶基因后不变。这与“维持该程序的候选上游因子”的定位一致，不支持“疾病中被激活的调控因子”这种说法。

## 2.4 伪 bulk 稳健性（[pseudobulk_fry_camera.tsv](results/pseudobulk_fry_camera.tsv)）

| 器官 / 类别 | 暴露 | 供体 | fry：方向，P | CAMERA：方向，P |
|---|---|---|---|---|
| 肺 Myeloid | COPD vs 对照（主） | 56 | 下调，0.115 | 下调，0.056 |
| 脑 Microglia 2 | 每 1 SD CPS_Local（主） | 79 | 上调，0.32 | 上调，0.37 |
| 脑 Microglia 2 | ADNC 中/高 vs 无/低 | 79 | 上调，0.31 | 上调，0.41 |
| 脑 Microglia 2 | CPS，≥ 50 个细胞核 | 78 | 上调，0.32 | 上调，0.37 |
| 脑小胶质总体 | 每 1 SD CPS_Local（主） | 79 | 上调，0.21 | 上调，0.31 |
| 脑小胶质总体 | ADNC 中/高 vs 无/低 | 79 | 上调，0.27 | 上调，0.35 |
| 脑小胶质总体 | CPS，≥ 50 个细胞核 | 79 | 上调，0.21 | 上调，0.31 |

- **CAMERA 与 fry 结论一致**：CAMERA 估计的集合内基因相关为 0.06–0.13。
- **50 个细胞核阈值**：Microglia 2 损失 1 名供体，小胶质总体和肺都没有损失；结果不变。
- **小胶质总体**：合并 Micro-PVM_2、2_1、2_3、3、4，结论与只用 Microglia 2 相同。
- **CPS_Local 的定义已核实**：来自 SEA-AD 官方表 `sea-ad_global_and_local_cps_050126.xlsx`（sheet "Travaglini 2026"，`Brain Region` = DFC），是 DFC 局部的病理进展评分，名称无需更正。
- **诊断（ADNC）**：取自 SEA-AD 供体元数据，与 CPS 分开报告。
- **未做**：肺泡巨噬细胞单独成类，需要重新读取 65 个原始 10x 矩阵。现有肺类别为全部髓系细胞。

## 2.5 临床分层（肺 Myeloid）

表型来自 GEO 样本信息台账（原作者分组）。参与伪 bulk 的 56 人中，E-COPD 12 人、NE-COPD 21 人、对照 23 人。

| 比较（同一 3 组模型 + 年龄、性别、吸烟） | fry：方向，P | CAMERA：方向，P |
|---|---|---|
| **E-COPD vs NE-COPD（直接检验两组差异）** | 上调，0.072 | 上调，0.18 |
| E-COPD vs 对照 | ≈ 0，0.97 | 下调，0.82 |
| NE-COPD vs 对照 | 下调，0.034 | 下调，0.018 |
| 每 1 SD FEV₁ % 预计值（全部 56 人） | 上调，0.91 | 上调，0.91 |

- 两个亚组之间的差异没有达到显著，不能说肺气肿型有特殊的细胞内 Axis A 改变。
- NE-COPD 与对照的比较显示细胞内略偏低，方向与总体结果一致。但这是 3 个探索性比较之一，未作多重校正，不宜单独强调。
- 原作者报告两种表型的巨噬细胞相对丰度不同（Zhang 2025 *Commun Biol*），这属于组成差异，与这里的细胞内检验是两个问题。

## 第 4 节：decoupleR 转录因子活性（[tf_activity_all.tsv.gz](results/tf_activity_all.tsv.gz)）

**方法**：decoupleR 2.16.0 `run_ulm`，CollecTRI 网络（43,062 条相互作用，1,186 个 TF），`minsize = 5`；输入为每个类别的 log-CPM；TF 活性按与上面相同的设计用 limma 检验。

| TF（CellOracle 名次） | 肺 Myeloid：疾病 t / P / 与 Axis A 的 r | 脑 Microglia 2：CPS t / P / r | 脑小胶质总体：CPS t / P / r |
|---|---|---|---|
| STAT3（1） | −0.58 / 0.56 / 0.53 | 0.94 / 0.35 / 0.55 | 1.53 / 0.13 / 0.56 |
| FOS（2） | 0.26 / 0.79 / 0.33 | 0.46 / 0.65 / 0.46 | 1.10 / 0.28 / 0.49 |
| LYL1（3，只有 10–12 个靶基因） | 0.80 / 0.43 / −0.25 | −2.07 / 0.042 / −0.32 | −2.34 / 0.022 / −0.30 |
| JUNB（4） | 0.99 / 0.33 / 0.27 | −0.47 / 0.64 / 0.28 | −0.09 / 0.93 / 0.29 |
| CEBPB（5） | 0.19 / 0.85 / 0.46 | 0.07 / 0.95 / 0.57 | 0.66 / 0.51 / 0.57 |

- **主家族**：在 3 个类别中，前十名 TF 的疾病或病理关联在家族内 FDR 都 ≥ 0.22。LYL1 在脑中的未校正 P 为 0.02–0.04，但方向为负，而且只有 10 个可测靶基因。
- **去掉 Axis A 靶基因后**，上述结果和相关系数基本不变（例如 STAT3 的 r 为 0.50–0.55），说明这种同步不是共同基因造成的机械相关。
- **全部 TF（次家族）**：
  - 肺中 664 个 TF 没有一个 FDR < 0.05；Microglia 2 中只有 1 个（HOXB7）。
  - 小胶质总体中有 15 个 TF 随 CPS 变化（例如 ZEB2、NR1H4、TP63、KLF8），都不在 CellOracle 前十名中。这些结果只作描述。

## 对稿件的含义

- **结果“细胞内”一节**：可以补一句“改用 CAMERA、小胶质总体及更严格的细胞数阈值后结论不变”。
- **肺气肿亚组**：放入补充材料，作为探索性结果，写明两组差异未检出统计学证据。
- **CellOracle 段**：可以写“在供体伪 bulk 中，STAT3 与 C/EBPβ 的推断活性与 Axis A 评分同步（r ≈ 0.5），但不随疾病或病理改变”。这与“维持该程序的候选上游因子”一致，也正好回应了修改意见第 4 节，说明这一判断有独立数据支持。

## 限制

- 肺类别为全部髓系细胞，没有单独分析肺泡巨噬细胞。
- 肺供体中 COPD 33 人（E-COPD 12、NE-COPD 21）、对照 23 人，亚组统计功效有限。
- **CollecTRI 的获取方式**：OmnipathR 3.18.4 的官方下载函数在 `unnest_evidences` 中出错，因此从 OmniPath 接口直接下载同一数据集，再按 decoupleR 2.16.0 `get_collectri` 的源码逻辑处理（`code/00_collectri.R`）。官方函数会追加 TF–miRNA 相互作用，本分析没有加入。原始文件和处理结果的 SHA-256 记录在 `source_data/collectri_sha256.txt`。
- TF 活性是基于先验网络的推断，同一数据换用其他调控先验只能作为计算交叉检查。

## 运行

从项目根目录依次运行：

```
R_LIBS=_workspace/r_libs_decoupler_20261019:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019/code/00_collectri.R
LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs_decoupler_20261019:_workspace/r_libs _workspace/envs/r/bin/Rscript 03_candidate_analysis/SingleCell_pseudobulk_supplements_20261019/code/01_pseudobulk_tf.R
```

decoupleR 2.16.0 与 OmnipathR 3.18.4 装在 `_workspace/r_libs_decoupler_20261019`（Bioconductor 3.22）。编译时需要 `SDKROOT=$(xcrun --show-sdk-path)`，并把 `_workspace/envs/r/bin` 加入 PATH。ADNC 由 SEA-AD 供体元数据导出到 `source_data/seaad_donor_ADNC.tsv`。
