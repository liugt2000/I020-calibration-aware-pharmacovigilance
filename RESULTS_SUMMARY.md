# I020 第一批真实结果摘要

- **数据源**：openFDA drug/event.json（N=21,054,159 报告），2026-11-21 查询
- **药物字段**：`openfda.generic_name`（锁定）
- **状态**：四算法 + 经验校准（N1/N2 bracket）+ 活性对照已跑通；网络证据/notoriety 未跑

## 1. 药物报告量（决定性约束）

| 药物 | n 报告 | 有效 N1 惰性对照 |
|---|---|---|
| osimertinib | 33,148 | 17/20 |
| amivantamab | 3,858 | 7/20 |
| lazertinib | 865 | 5/20 |
| tepotinib | 691 | 4/20 ← 不足以拟合稳健 null |

**方法学发现 M1**：经验校准的可用性随药物报告量急剧下降。只有 osimertinib 能满足 ≥10 有效阴性对照的方法学门槛；tepotinib 的 N1 null 无法拟合（其 AKI 校准值为 NaN）。这本身是"FAERS 灌水"批评的一个具体体现：小样本新药的 ROR 看似惊人（tepotinib-AKI ROR=8.68），但连最基本的经验校准都做不到。

## 2. Index pairs 主结果

| Pair | 朴素 ROR (95%CI) | 四算法一致 | 校准 ROR (ML, N1) | q (N1) | 活性对照 ROR (95%CI) | 判定 |
|---|---|---|---|---|---|---|
| osimertinib–ILD | 7.17 (6.75–7.62) | 是 | 16.53 (LB 10.83) | <0.001 | **0.86 (0.79–0.93)** | 校准后衰减→类效应 |
| osimertinib–QT | 3.80 (3.33–4.33) | 是 | 8.75 (LB 5.65) | <0.001 | **3.51 (2.62–4.69)** | 保留→药物特异 ✅ |
| tepotinib–AKI | 8.68 (6.91–10.91) | 是 | NaN (N1 不足) | — | **2.92 (2.17–3.94)** | 活性对照保留，但校准不可评估 ⚠️ |
| lazertinib–ILD | 12.80 (9.32–17.58) | 是 | 22.50 (LB 16.38) | <0.001 | **1.27 (0.92–1.75)** | 跨 1→衰减为类效应 |
| amivantamab–IRR | 80.83 (74.69–87.48) | 是 | 150.2 (LB 116.5) | <0.001 | **1541.8 (96.3–24683)** | 对照 0 报告→极端值，需注明 |
| amivantamab–ILD | 7.83 (6.48–9.45) | 是 | 14.54 (LB 10.70) | <0.001 | **1.36 (1.00–1.86)** | 边缘跨 1→多为类效应 |

## 3. 关键方法学发现（按预设命题核对）

- **C1 ✅**：6 个 index pairs 全部四算法一致——但阳性对照外的阴性对照也被大量标记（见 null params），印证"四算法一致 ≠ 真实性"。
- **C2 ⚠️ 需修正表述**：校准后 ROR **点估计上升**（如 osimertinib-ILD 7.17→16.53），因为 μ 为负（惰性对照 under-reported）。衰减体现在**活性对照后**和**置信区间展宽**，而非校准本身降低点估计。这印证了我之前的判断，必须写清楚。
- **C3 ✅**：类效应事件（ILD）活性对照后衰减到 ~0.86–1.36；药物特异事件（QT、IRR）保留（3.51、1541）。预设命中。
- **C4 ❌ 预设未通过（重要阴性发现）**：网络邻近在 robust 与 attenuated 组间**无判别力**（mean z 差 −0.258，exact label-permutation p=0.40，两个网络阈值一致）。唯一表观显著的 tepotinib–AKI（z=−2.54, p=0.021）经敏感性分析证实为**循环论证**——去掉从已发表文章借来的 6 个 hub 后 z 翻转到 +0.67、p=1.0。详见新 §7。这不是失败，而是强化"多源网络整合需防循环偏倚"的方法学证据，须改写 Results §5.6 与 Discussion ¶5。
- **C5 ❌ 预设未通过（阴性发现）**：notoriety 时间分析未见任何 notoriety 模式——所有 change point 的 index−control 报告占比增量为负或可忽略（最大正增量仅 +0.008，无监管锚点吻合）。与经验校准结论一致：时序证据也不支持"通告驱动的报道膨胀"。注意 tepotinib 的 N1 对照序列全为 0、amivantamab 仅 11 例，时序分析对这两药信息量有限。详见 §8。
- **额外 M2**：活性对照揭示 osimertinib-ILD 的 ROR 从 7.17 跌到 0.86——即该信号在 EGFR-TKI 类内**根本不特异**，先前发表的高 ROR 主要是"EGFR-TKI 整体 vs 全库"的类效应。这是对 FAERS 灌水批评的最强实证回应。

## 4. null 参数（经验校准核心）

| 药物 | 面板 | n 对照 | μ (ML) | σ (ML) | μ (sample-SD) | σ (sample-SD) |
|---|---|---|---|---|---|---|
| osimertinib | N1 | 17 | -0.836 | 0.214 | -0.943 | 0.631 |
| osimertinib | N2 | 15 | 0.538 | 1.435 | 0.538 | 1.241 |
| lazertinib | N1 | 5 | -0.564 | 0.000* | -0.414 | 0.441 |
| amivantamab | N1 | 7 | -0.620 | 0.123 | -0.600 | 0.559 |

\* σ=0 是边界解，样本太小，不可靠。

## 5. 阳性对照保留率（校准未过度收缩的证据）

- osimertinib 5/5、lazertinib 5/5、amivantamab 4/5 → 校准保留了 genuine 事件，方法未过度收缩。这是重要的正面验证。

## 6. 写作影响

- 主 claim 从"多数信号校准后衰减"调整为更精确："**经活性对照重估，类效应信号衰减，药物特异信号保留；经验校准在低报告量药物上不可行，需报告适用边界**"。
- 骨架 Methods §4.5 需要补充"最低有效对照数门槛"和"邻界 σ 处理"规则。
- 图 3 瀑布图数据已就绪（results_active.csv）。
- **C4 阴性结果改变了文章的一处核心叙事**：不能再宣称"存活信号富集于网络机制证据"（Abstract Conclusions 与 Discussion ¶5 原话须删改）。改为："网络邻近分析未能在去偏后提供独立判别力；未控制 hub/循环偏倚时网络证据会产生假阳性——这本身就是对多源 FAERS 整合的警示。"主 claim 收缩为统计去偏（校准+活性对照），网络分析降级为"负面方法学对照"，仍扣住专刊 Network Medicine 主题（讨论网络医学方法的陷阱）。
- 分组标签仍由校准结果固定，符合冻结纪律；C4 属结果而非 prespecified success。

## 7. C4 网络邻近：完整阴性结果与循环论证诊断（2026-11-21）

### 7.1 主分析（STRING 诱导子网络，degree-matched null 10,000 次）

| Pair | 组 | d_c (score≥0.4) | z | p_emp | d_c (score≥0.7) | z |
|---|---|---|---|---|---|---|
| osimertinib–ILD | attenuated | 0.750 | −0.08 | 0.650 | 0.750 | −0.23 |
| osimertinib–QT | robust | 0.500 | −1.00 | 0.306 | 0.500 | −0.89 |
| tepotinib–AKI | robust | 0.143 | **−2.54** | **0.021** | 0.143 | **−2.94** |
| lazertinib–ILD | attenuated | 0.500 | −0.96 | 0.402 | 0.500 | −1.03 |
| amivantamab–IRR | robust | 1.000 | +0.80 | 1.000 | 1.000 | +0.84 |
| amivantamab–ILD | attenuated | 0.500 | −0.93 | 0.399 | 0.500 | −1.04 |

判别性检验：robust−attenuated 平均 z 差 = **−0.258**，精确 label-permutation（C(6,3)=20）单侧 p = **0.40**。→ C4 预设失败。

### 7.2 循环论证敏感性分析（决定性证据）

| 场景 | 靶点集 | 毒性基因集 | tepotinib–AKI z | p |
|---|---|---|---|---|
| S1 基线 | union（含发表 hub） | union（含 hub） | −2.54 | 0.020 |
| S2 去靶点 hub | primary（仅 MET） | union | **+0.67** | **1.000** |
| S3 去毒性 hub | union | clean（器官特异） | +0.62 | 0.891 |
| S4 双去 | primary | clean | +1.32（z 失效，null 方差→0） | 1.000 |

**解读**： tepotinib–AKI 的表观网络邻近完全依赖靶点侧的 6 个已发表 hub（HSP90AA1/AKT1/EGFR/CASP3/TNF/SRC）——这些 hub 正是 PMID 41582009 那篇文章的结果，用它来"验证"该信号属于循环。去掉后信号消失。S4 的 z=40.8（osimertinib–QT）是单基因靶点集导致 null 方差趋零的统计假象，不可用。

### 7.3 方法学结论（写入 Discussion）

1. 网络邻近度在本批 NSCLC 靶向药上**无独立判别力**：多数 index pair 的靶点集是单基因（EGFR 或 MET），network proximity 的分辨率与 null 稳定性不足。
2. 网络证据的表观显著性由**hub 基因选择**和**靶点集规模**驱动，而非毒性特异机制；这是网络医学整合的已知陷阱（degree-matching 只能部分缓解）。
3. amivantamab–IRR 的 z 为正，反映网络邻近对**免疫介导/输液反应**机制不适用——距离度量无法捕捉 Fc 受体/肥大细胞通路。
4. 正确定位：**网络分析作为负面方法学对照**，展示"未去偏的多源网络整合如何制造假证据"，而非作为阳性机制佐证。这比假装成功更有说服力，且直接服务专刊的计算建模/网络医学主题。

数据文件：`network_data/gene_sets.json`、`network_results.json`、`sensitivity_c4.json`。脚本：`build_gene_sets.py`、`network_proximity.py`、`sensitivity_c4.py`（固定种子 20261121）。

## 8. C5 Notoriety 时间分析（阴性结果，2026-11-21）

脚本 `notoriety.py` → `notoriety_results.json`。方法：每季度 index 事件报告占比减去同药 N1 对照占比（`(a−c)/n`），data-driven change-point（每侧≥4季度），对照 FDA 监管锚点季度。

| Pair | Change point | Δ (index−control fraction) | 监管锚点 | 判读 |
|---|---|---|---|---|
| osimertinib–ILD | 2019Q4 | −0.0533 | 2018Q2 | 负向，无 notoriety |
| osimertinib–QT | 2024Q3 | −0.0100 | 2018Q2 | 可忽略 |
| tepotinib–AKI | 2023Q1 | −0.1800 | 2021Q2 | 负向；N1 对照序列全 0，信息有限 |
| lazertinib–ILD | 2025Q3 | +0.0075 | 无 | 微小正增量，无锚点吻合 |
| amivantamab–IRR | 2022Q3 | −0.1728 | 2021Q2 | 负向；对照仅 11 例 |
| amivantamab–ILD | 2022Q4 | −0.0383 | 2021Q2 | 负向 |

**结论**：无一例出现"与监管锚点吻合的正向 notoriety 阶跃"。最大正增量 +0.008 且无锚点。与经验校准方向一致——时序证据同样不支持这些信号是通告驱动的报道膨胀。锚点说明：osimertinib ILD/QT 警示在 2015 初次标签即存在，无单一 sentinel 事件；tepotinib/amivantamab 2021Q2 加速批准。

## 9. 引用核对（致命错误已修正）

**重大更正**：原骨架把 osimertinib–ILD 网络毒理文章标为 PMC13343941——经 NCBI eutils 核实，该 PMCID 实为 1956 年德文《男性生育力测定》文章，完全无关。正确引用：**PMID 42106642, BMC Cancer 2026;26, doi:10.1186/s12885-026-15993-5**。若不查证直接投稿会被审稿人一击致命。

其余修正：PMID 42653711 卷期 19(1)→**19(8):1213**；PMID 41582009 期刊确认为 *Ren Fail* 2026;48(1):2610895、作者 Zheng R 等；PMID 41563236 期刊确认为 *Int J Surg*；PMID 34396179 卷期确认 1(3):347–356。全部引用已写入 v1.0 全文 §7。

## 10. 交付状态（投稿级完整稿件）

- **主交付物**：`../I020_manuscript_v1.0.md`（完整英文全文，所有占位符已用真实数据填充，Methods 写实际执行的数据源，C4/C5 阴性结果如实呈现，引用全部修正）。
- **图表**：`figures/Fig1_flow.png`、`Fig2_calibration.png`、`Fig3_waterfall.png`、`Fig4_notoriety.png`、`Fig5_network.png`、`Fig6_readuspv.png`（均 300 dpi，已目视质检）。
- **数据映射**：`fill_values.json`（全部占位符计算结果）+ `summarize_for_manuscript.py`（可复算）。
- **v0.1 骨架**保留为 prespecified 历史记录，已加 v1.0 指针与引用更正。
- **作者仍需完成**（列为 v1.0 §8 checklist）：外部时间戳存档拿公开 URL、作者/机构/基金等 metadata、Table 1–6 格式化、Int J Surg 卷期终校。

## 11. 审稿视角三项提升（2026-11-21，v1.1）

专家审稿视角识别出三项可提升点，已按顺序全部完成：

**① suspect-only 敏感性分析（数据重跑）**
- 问题：暴露字段未区分 suspect/concomitant，且聚合计数无法病例级去重（Fig6 原 note 自曝）。
- 执行：`fetch_suspect.py` 对所有 ~359 查询加 `drugcharacterization:"1"` 约束重跑；`pv_core.py` 复用核心逻辑；`compare_suspect.py` 逐项对比。
- 结果（决定性）：暴露计数最大变化 0.23%，index ROR 最大变化 2.38%，**全部 robustness call、null 参数、活性对照衰减判定、阳性对照保留（14/15）零变化**。根因：openFDA `generic_name` 字段主要为 suspect 药物填充，两种暴露定义近乎共线。
- 结论：concomitant 污染不是本结果的实质驱动因素——从"潜在硬伤"降级为"已量化并化解的疑虑"。但 case-level 去重仍不可行（openFDA 每报告唯一 ID，更新作为独立报告），已写入 Limitations 并注明 ASCII 重分析为对症补救。
- 新增稿件内容：§2.2（方法）、§3.2（Suspect-role sensitivity 小节）、Limitations、Fig6 note 更新、S2b。

**② 补充材料与表格实体化**
- 问题：全文引用 Supplementary S1–S10 与 Table 1–6，但文件全部不存在（虚构引用，形式审查即退回）。
- 执行：`supplementary/I020_supplementary_S1-S10_v1.0.md`（S1 冻结计划/S2 PT与对照/S2b suspect/S3 监管锚点/S4 人口学/S5 四算法/S6 时序/S7 网络/S8 null拟合/S9 代码地图/S10 READUS-PV checklist）+ `supplementary/I020_tables_v1.0.md`（Table 1–6 全部用真实数据填充）。S4 人口学数据为实时抓取的性别/年龄分布。

**③ Abstract 压缩 + Tier 逻辑自洽**
- 问题：Abstract 513 词超 JTM 上限；Tier A=QT+IRR 与正文"IRR unstable / AKI retained"矛盾。
- 执行：Abstract 压至 286 词（保留全部关键数字与阴性发现）；tier 规则加入"估计稳定性"判据（comparator 事件≥1 且 ROR 上限 CI<100），IRR（c=0, UB=24,683）降入 Tier C。
- 新 tier：**A=QT(1) / B=三个ILD(3) / C=AKI(校准不可估)+IRR(估计不稳定)(2)**。已同步 fill_values.json、§3.4/§3.7、Discussion ¶1、Table 6、Translation、Abstract。
- 稿件版本升至 v1.1。

**新增产物**：`fetch_suspect.py`、`pv_core.py`、`compare_suspect.py`、`counts_suspect.csv`、`counts_active_suspect.csv`、`suspect_*.csv`(4)、`suspect_comparison.json`、`supplementary/` 两份文件；Fig6 重新生成。**0 占位符、0 VERIFY、Abstract 286 词。**

## 12. 校准仿真：补强阳性对照论证（2026-11-21，v1.2）

**问题**：原"14/15 标签阳性对照保留"是空论证——强效应几乎必然通过，证明不了"校准无过度收缩"。

**执行**：`pc_simulation.py`（20,000 replicates，种子 20261121）→ `simulation_results.json`；`make_fig_sim.py` → `figures/Fig7_simulation.png`。DGP：`ln(ROR)_obs ~ N(μ+δ, σ²+SE²)`，SE 从 6 个真实 index pair 经验重采样。两套 null：N1 primary（μ=−0.836, σ=0.214）、N2 upper（μ=+0.538, σ=1.435）+ ideal 参考。三实验：A true-null type-I error、B power curve（ROR 1.2–6）、C 15 个临界对照（ROR 1.4–2.5）。

**结果（决定性）**：
- A：naive FPR **99.4%/91.6%**，校准 **5.1%/5.0%**（BH-FDR naive 99.4%/91.5%，cal 0/0）。→ 朴素检验在真实报告过程下结构性失效，这是"FAERS 灌水"的统计机制。
- B：主 N1 null 下校准效能 ROR2=**84%**、ROR3=99%、ROR≥4≥99%（vs ideal 63%/94%/100%）；N2 下 ROR6 仅 24%。
- C：主 null 校准保留 **9/15** vs naive 8/15（无系统性收缩）；N2 下校准 0/15 vs naive 14/15——14 个全是假阳性，暴露 naive 正对照论证的虚假性。

**结论**：校正是"选择性去偏"而非"破坏性收缩"——条件是阴性对照面板质量好（μ/σ 小）；N2 高离散面板下校准正确保持保守。这条把 14/15 从断言升级为量化属性，并给出适用边界。

**稿件变更**：§2.5（仿真设计）、§3.3（结果段）、Discussion 新增段、Abstract +1 句（318 词）、Fig. 7 图注、S11（完整表与复现命令）。版本 v1.2；哈希刷新；0 占位符/0 VERIFY。

## 13. 第二档三项收尾（2026-11-21，v1.3）

**① 版本锁定段（§2.8 + S12）**：查证并锁定全部外部源版本——openFDA 快照（2004→2026-11-21，N=21,054,159，API 不暴露 MedDRA 版本及固定截止季度，如实列为 reproducibility limitation）、ChEMBL_37（2026-05-01）、STRING 12.5（v12.0 仅为早期弃用下载文件名，未用于任何分析）、Open Targets API 26.9.0 / data 2026-09、Python 3.13.14 + pandas 2.3.3/numpy 2.5.4/scipy 1.18.1/networkx 3.7/matplotlib 3.11.2/requests 2.34.2、固定种子 20261121（仿真+网络 permutation）。S12 含逐项验证命令。

**② 跨库 replication 路线（Discussion + S13）**：基于核实的真实访问方式——EudraVigilance 通过 adrreports.eu 公开聚合报告（按物质×反应 PT，含年龄/性别/严重度/地区，每月更新）可重建 2×2 并用同一 N1/N2 PT 面板重拟合 null；VigiBase 公开层 VigiAccess（聚合计数，无年龄/国家明细）+ 向 UMC 申请 case-level（access level 2A limited / 2B full，no free text，签 Caveat，非成员收费），后者可直接实现我们缺失的 case-level 去重；第三条为 FAERS 季度 ASCII 的 DEMO/DRUG/REAC case + FDA_DT 去重重分析。S13 含 6 条最小复现规范（同一 PT 清单、本地标签重新 adjudication、固定 estimator、报告 μ/σ、同类活性对照、固定查询日期）+ 6 个核实过的源 URL。预期定性模式可复现，绝对 ROR 不要求一致——检验的是方法而非数值。

**③ 决策后果量化（§3.8 + Fig8 + S14）**：`decision_impact.py` → `decision_impact.json`；`make_fig_decision.py` → `figures/Fig8_decision.png`。决策漏斗（可触发 agent-specific 行动的 pair 数）：naive 四算法 6 → 校准可估 5（tepotinib-AKI 出局）→ 活性对照 2（三 ILD 降类效应）→ 稳定性 Tier A 1（osimertinib-QT；IRR 零分母出局）= **83% 减少**。结合 S11 仿真 FPR：每筛 6 对期望假行动 naive 5.73（FPR 95.5%）vs 校准 0.30（FPR 5.05%）≈ **19 倍降低**。明确标注为 decision-model 计算，非实际监管事件计数；类效应重分类仍保留 class-level 监测。已写入 Abstract/Translation/§3.8/Fig8 图注/S14。

**交付状态**：v1.3，全文 6810 词，Abstract 340 词（<350），0 占位符/0 VERIFY，8 张 300dpi 图，补充材料 S1–S14 + Table 1–6，哈希 manifest 一致。审稿视角第一档三项（v1.1）+ 第二档仿真（v1.2）+ 第二档三项（v1.3）全部清零。
