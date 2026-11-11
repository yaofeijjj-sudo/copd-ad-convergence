# R7 条件 MAGMA：以一般小胶质细胞基因集为条件（2026-10-16，候选，已完成）

设计在运行前冻结：[DESIGN.md](DESIGN.md)（SHA-256 见 `DESIGN.sha256`）。作者 10-16 授权。

## 问题

Axis A 富集 AD 遗传信号，是否只是因为它属于小胶质/髓系基因（替代解释 R7：“髓系基因本来就是 AD 风险基因”）？

## 方法

- MAGMA v1.10 竞争性基因集分析，`--model condition=MicrogliaRef`，把“是否属于一般小胶质基因集”作为协变量（手册 v1.10 的 Conditional and joint analysis；de Leeuw 2018 *Nat Commun*；Jansen 2019 以 APOE 基因集为条件的同一用法）。
- 基因层面输入：SYN01 冻结的 `AD_total_N.genes.raw`（18,224 个基因）；基因集：SYN01 `locked_axes.sets`。
- 参照集：原计划用 Patir 2019 *Glia* 的 249 基因核心签名，但 Wiley 补充文件被 Cloudflare 验证拦截（未绕过），按作者同意改用 BRETIGEA 1.0.3 的人类小胶质标志基因（McKenzie 2018 *Sci Rep*）。主分析取前 250 个（与原计划规模一致），敏感性分析用全部 1,000 个；符号按 org.Hs.eg.db 官方符号唯一映射到 Entrez。

参照集（[reference_set_manifest.csv](results/reference_set_manifest.csv)）：

| 参照集 | 映射到 Entrez | 进入基因分析 | 与 Axis A 重叠（进入分析） |
|---|---:|---:|---:|
| 前 250 | 217 | 202 | 47 |
| 全部 1,000 | 775 | 652 | 66 |

## 结果（[conditional_magma_results.csv](results/conditional_magma_results.csv)）

| 模型 | 参照集 | Axis A β（SE） | P | Holm P（两轴） | 小胶质参照集 P |
|---|---|---|---:|---:|---:|
| 无条件 | — | 0.364（0.084） | 7.2 × 10⁻⁶ | 1.4 × 10⁻⁵ | 0.042（前 250）；0.0040（1,000） |
| 以小胶质集为条件（主） | 前 250 | 0.351（0.087） | 3.0 × 10⁻⁵ | 6.0 × 10⁻⁵ | 0.30（与 Axis A 同模型） |
| 以小胶质集为条件（敏感性） | 全部 1,000 | 0.331（0.086） | 5.9 × 10⁻⁵ | 1.2 × 10⁻⁴ | 0.041（与 Axis A 同模型） |

- Axis B 在各模型中均未检出（P = 0.12–0.13）。
- 描述性：Axis A 中不属于参照集的成员（前 250 时 83 个、1,000 时 64 个）单独检验仍富集（P = 0.0013、0.0014）。

**结论（按设计中预先写定的规则）：** 两种参照集下，Axis A 的条件 Holm P 均 < 0.05。Axis A 的 AD 遗传富集不能由“属于一般小胶质基因集”完全解释。反过来，前 250 个小胶质标志基因本身的 AD 富集较弱（P = 0.042），与 Axis A 同时放入模型后不再显著（P = 0.30）。

## 局限

- 参照集来自表达特异性而非遗传学；前 N 个的截断是事先固定的人为选择。Patir 2019 未能获取，可在作者提供文件后作为补充检验重跑。
- 结论为关联层面；与 R7 的另一半（APOE 区域排除，见 `AxisA_MAGMA_APOE_exclusion_20261016/`）合起来，Axis A 的 AD 富集既不依赖 APOE 区域，也不依赖一般小胶质基因背景。

## 复现

```sh
LC_ALL=en_US.UTF-8 R_LIBS=_workspace/r_libs_deconv_20261002:_workspace/r_libs Rscript 03_candidate_analysis/AxisA_MAGMA_conditional_microglia_20261016/code/run_conditional_magma.R
```

首次运行时清单中“进入基因分析”的计数因 `genes.raw` 行长不一而解析错误；已改为逐行取第一字段并断言 18,224 个基因后重跑。MAGMA 统计量不受影响（两次运行逐值相同）。
