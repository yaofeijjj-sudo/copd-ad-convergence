# T04 AD 年龄异常样本固定模块敏感性：运行前分析合同

状态：**已锁定，尚未查看敏感性模型结果**

锁定日期：2026-07-24

## 要回答的问题

GSE33000 current AD 网络在 AD/对照组内按 1.5×IQR 排除了 7 个年龄异常样本。T04 检验：在不重新定义网络模块的前提下，将这 7 例重新纳入，是否会实质改变 AD blue、AD turquoise 及其他 fixed module 的年龄/性别调整疾病效应。

## 固定不变的分析对象

- 继续排除 Huntington's disease 样本。
- current 参考网络为 `AD_WGCNA_no_residual_after_age_outlier`。
- 保留 current 的 8,807 个基因、17 个颜色标签（16 个非灰模块加 grey）和每个基因的模块归属。
- 保留 current 的探针选择、表达尺度、模块 PC1 载荷及方向。
- 保留 current 已独立识别的连接度异常样本 `GSM1424097` 的排除，不把它与年龄敏感性混在一起。
- 疾病效应模型保持为 `module eigengene ~ Disease_State + age_z + sex_male`；`age_z` 在各自分析队列内标准化。由于带截距模型中年龄的线性重标度不改变疾病系数，该疾病系数也等价于使用年龄（年）调整。
- current 分析为 459 例；T04 固定模块敏感性预计为 466 例，即恢复 7 个年龄异常样本但仍排除 1 个连接度异常样本。

## 固定模块投影

1. 从 current 年龄筛选后、连接度筛选前的 460 例学习分位数标准化目标；不因重新加入 7 例而改变参考样本的表达值。
2. 使用 current 已冻结的 selected probe 和 8,807 个基因。新增 7 例只映射到该固定特征空间。
3. 对每个模块，仅用 current 459 例学习基因中心、标准差、PC1 载荷、奇异值及方向；新增样本应用相同参数作样本外投影。
4. 在进入敏感性模型前，投影得到的 current 459 例 eigengene 必须复现保存的 current eigengene：每个模块 Pearson `r >= 0.999999` 且最大绝对误差 `<= 1e-10`。
5. 不重新切树、不重新合并模块、不重新选择软阈值，也不利用敏感性结果改变 Axis A/B 基因成员。

## 预先锁定的升级判据

以下任一条件满足，则 T04 判为实质不稳定，需要暂停 T05 及以后任务，并升级为同参数完整 WGCNA sensitivity；否则 current 网络保持主分析，完整重建不启动。

### Axis A/B 主判据

对 AD blue（Axis A 脑侧）或 AD turquoise（Axis B 脑侧）：

1. 疾病系数方向反转；
2. 95% 置信区间跨 0，或跨 17 个固定模块校正后的 BH FDR 不再小于 0.05；
3. 相对 current 的绝对疾病效应衰减超过 25%；
4. 在 16 个非灰模块中按绝对疾病效应排序，由 current 前 3 降至前 5 之外。

### 全模块辅助判据

对 16 个非灰模块的疾病系数向量：

1. current 与敏感性结果的 Pearson 相关低于 0.90；或
2. current 中 `|beta| >= 0.01` 的模块方向一致率低于 80%。

全模块判据用于发现广泛失稳；不能用它掩盖任一 Axis 主判据失败。

## 预设报告内容

- 7 个恢复样本的 ID、组别、年龄、性别和原始 IQR 排除原因；
- current 459 例与 fixed-module 466 例的全部模块效应、标准误、95% CI、P 值、BH FDR 和非灰模块排名；
- blue/turquoise 的方向、绝对及相对效应变化；
- 投影复现验收、样本与基因计数验收；
- 预设触发器逐项 PASS/FAIL；
- go/no-go 决策及其对 T05 的影响；
- Methods、Results 和 limitation 的中英文同步表述；
- 可追溯的输入、代码和输出哈希。

## 图形合同

- 核心结论：恢复 7 个年龄异常样本后，固定 AD 模块的协变量调整疾病效应是否保持。
- 图形类型：双面板 quantitative robustness figure。
- Panel a：blue/turquoise 在 current 与敏感性分析中的疾病系数及 95% CI。
- Panel b：16 个非灰模块 current 与敏感性疾病系数的对应关系，突出 blue/turquoise，显示恒等线。
- 图内不放总标题；Axis 标签直接标注；中性模块使用低饱和灰色。
- 输出 R 生成的 SVG、PDF 和 600 dpi TIFF，并提供逐点 source data。

