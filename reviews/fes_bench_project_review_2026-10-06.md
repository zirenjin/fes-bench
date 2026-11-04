# fes-bench 项目联合审查

日期 2026-10-06

## Review setup

- 输入范围包括项目叙述、下载目录中的 Colab 导出脚本、仓库内 notebook、GitHub 对应仓库、结果表、原始指标与 provenance 元数据。
- 评估边界是研究计划与部分证据包，不是完整投稿稿件。未提供完整引言、系统性的相关工作、冻结的已执行 notebook、完整图注或提交版讨论。
- 共享主张是准谐基线能够使学习到的自由能面跨温度、相和体系外推。
- 可见证据包括 4 个体系、9 个相、6 个相对自由能对，三类 OOD 切分，QH 可靠性审计，温度外推的学习模型结果，以及尚不完整的 LOPO 和 LOSO 结果。
- 三份相互盲审报告分别保存在 `reviewer_1.md`、`reviewer_2.md` 和 `reviewer_3.md`。审查重点在生成前固定，报告冻结后才进行本综合和数值一致性审计。

## 结论先行

当前项目最有价值的结果不是题目所写的“QH 基线让模型完成跨域外推”，而是一个更可信、也更有科学张力的负结果。

> QH 先验只有在适用域和任务目标同时匹配时才可能简化自由能学习。对当前数据，QH-only 明显失败，QH-residual 没有形成稳健优势，受约束的 thermo-form 反而显著更强。

现有证据不能支持原标题。温度外推中，QH-residual 的 ΔG MAE 为 31.239 meV/atom，representation tlog 为 31.557 meV/atom，zero baseline 为 4.521 meV/atom，thermo-form 为 0.401 meV/atom。QH-residual 的 ranking accuracy 仅 0.1045，并在 15 个 seed-pair crossing 事件中漏掉 12 个。LOPO 和 LOSO 还没有对应的学习模型结果。

因此当前状态应判定为“内部研究草稿，必须重构中心论点”，不宜继续以原标题推进。

## nature-proposal-writer QA

| 维度 | 分数 | 主要依据 |
| --- | ---: | --- |
| 研究问题清晰度 | 5.5 | 题目给出方向，但没有把可证伪的主假设、主要终点和成功阈值写成独立契约 |
| 科学张力 | 8.0 | 小 ΔG 对 Tc 极敏感，物理先验可能失效，学习模型与 thermo-form 形成真实冲突 |
| 证据匹配 | 3.5 | 标题覆盖三类 OOD，但学习证据只有温度 OOD，而且结果不支持正向主张 |
| 逻辑链 | 4.5 | 数据审计较强，但目标、消融、适用域、结果和结论没有闭环 |
| 方法可行性 | 5.5 | 仓库结构较成熟，但关键训练、PES 验收、TDEP、TI 上界和执行环境未冻结 |
| 创新性 | 6.0 | benchmark 与失败域分析有潜力，然而缺少 claim-matched prior-work 对比 |
| 风险边界 | 8.0 | 对泄漏、QH 失效、Hf 温度混杂和未完成任务披露充分 |
| 语言质量 | 5.0 | 项目叙述是工作笔记形态，notebook 版本不一致，术语密度高且中心命题漂移 |
| 平均分 | **5.75** | 低于 proposal 阈值 7.0 |

状态为内部草稿。建议先重构论证，再补实验。当前不应进入语言润色或投稿包装阶段。

## Claim to evidence 审计

| ID | 关键主张 | 当前证据状态 | 处理建议 |
| --- | --- | --- | --- |
| C01 | QH baseline 使模型跨温度、相和体系外推 | 不支持 | 删除或改写原标题 |
| C02 | QH-residual 在温度外推优于 representation-only | 证据混合 | MAE 仅小幅改善，ranking 和 crossing 更差，不能作为总体胜出 |
| C03 | QH-only 是有效物理基线 | 反证 | 三个 split 的 skill 均大幅为负且 crossing 全部或多数失败 |
| C04 | 1 meV/atom 可导致 SiO2 Tc 偏移约 130 K | 支持 | 由 7.449 至 7.551 micro-eV atom^-1 K^-1 的斜率可复核 |
| C05 | 高温相普遍落在 QH 不可靠域 | 支持于当前样本 | 5 个相超过 1% imaginary-mode 阈值 |
| C06 | QH 使 Hf、Ti、Zr 残差更容易学 | 依赖指标 | no-learning ratio 接近 1，但 learned anchored-shape ratio 显著降低，必须区分目标和拟合器 |
| C07 | 单侧 QH 可处理 mixed-regime phase pair | 不支持 | pair-level ΔG MAE 为 644.868 ± 465.591 meV/atom |
| C08 | thermo-form 是强基线 | 数值支持 | 但训练信息、自由度和泄漏约束仍需与学习模型完全对齐 |
| C09 | LOPO 已验证 learned phase extrapolation | 不支持 | learned rows 为 `n/a:not_trained_for_split` |
| C10 | LOSO 已验证 learned system extrapolation | 不支持 | learned rows 为 `n/a:not_trained_for_split`，Hf 还混入 temperature OOD |
| C11 | 三个 seed 足以建立稳定性 | 未建立 | aggregate 方差近乎为零，但 checkpoint 不在本地且 commit 为 unknown |
| C12 | repository 已形成可复现包 | 部分支持 | 结构与内部检查较好，但工作树 dirty、notebook 未执行、输入 hash 不一致 |
| C13 | 当前四体系覆盖主要 entropy mechanisms | 不充分 | 叙述中的 strong、partial、absent 分类没有独立分解证据 |
| C14 | frozen DPA 可作为足够精确的静态基线 | 未建立 | pairwise static error 为 16.3 至 38.3 meV/atom，PES fine-tune 验收缺失 |
| C15 | 结果可归因于 QH residual learning | 未建立 | 上游 PES error、reference construction、functional sensitivity 和 calibration 均可能混杂 |

## Cross-review synthesis (post-review; not shown to reviewers)

### Consensus strengths

三位审查者一致认可下列优点。

- 三类 OOD 切分是合理的研究框架，优于把所有 holdout 混成一个泛化指标。
- crossing、sign、ranking、missed crossing 和 coverage 比单独 MAE 更贴近相稳定性任务。
- QH reliability audit 和 single-side negative control 是项目中最有洞察力的部分。
- 项目对负结果、泄漏风险和未完成任务披露充分。
- thermo-form 的强结果值得上升为中心科学发现或主要对照。

### Consensus blocking concerns

1. **中心主张与结果冲突**
   关联 R1-M1、R2-M1、R3-M1。现有温度外推结果没有证明 QH-residual 的稳健优势，phase 和 system OOD 更没有学习结果。

2. **三类 OOD 的证据链不完整**
   关联 R1-M2、R2-M3、R3-M1。必须在完全冻结且匹配的 protocol 下完成 QH-only、QH-residual、representation-only 和 thermo-form。

3. **QH 适用域是核心限制**
   关联 R1-M3、R2-M4、R3-M2。5 个相不通过 1% 阈值，mixed-regime 结果极差。最终论断必须按 both-reliable、one-side-reliable 和 neither-reliable 分层。

4. **thermo-form 比较尚未解释**
   关联 R1-M6、R2-M1、R3-M3。若 thermo-form 在同等无泄漏条件下仍达到 0.401 meV/atom 和 2.796 K，它将改变论文主线。

5. **上游 PES 与 reference error budget 缺失**
   关联 R1-M4、R2-M9、R3-M6。当前 DPA 静态能误差与目标 ΔG 尺度相当或更大，不能直接把失败归因于 residual learner。

6. **复现包尚未冻结**
   关联 R1-M5、R2-M7、R3-M4。需要一个 clean commit、完整 checkpoint 或可持续引用、已执行 notebook、producer、环境和全输入 hash。

7. **crossing 与 uncertainty 汇总不透明**
   关联 R1-m1、R2-M6、R3-M7。Tc error 必须与成功率、missed/false crossing 和 denominator 同表报告。

8. **跨体系推断范围过宽**
   关联 R1-M2、R2-M8、R3-M5。LOSO 只有三个相关金属体系，系统是推断单元，不能用 frame 数量放大证据强度。

### Other consensus major concerns

- 原创性仍不可判断。三位审查者均要求 claim-specific prior-work comparison。
- 泄漏审计和 Hf temperature confounding 必须关闭后才能解释 OOD 结果。
- functional 和 truth sensitivity 缺失，使 metals 与 SiO2 的聚合比较不稳固。

### Where emphasis differs across reviewers

- Reviewer 1 最强调方法归因、PES error budget、训练公平性和复现闭环。
- Reviewer 2 最强调原标题与证据冲突、原创性边界以及 thermo-form 可能带来的研究主线转向。
- Reviewer 3 最强调系统级推断、crossing coverage 对非专业读者的误导风险，以及模型层级和术语对跨学科读者不够透明。

### Forensic consistency findings

#### AUD-01

- 状态 `confirmed_internal_error`
- 位置 README lines 108-110，QH-residual per-seed metrics
- 检查 raw seed 11、23、37 的五个 reference crossing
- 结果 每个 seed 都命中 SiO2 quartz-cristobalite，漏掉另外 4 个。总计 12 个 missed events。README 的 “all five crossing pairs have missed predicted crossings” 不成立。
- 修正 将表述改为“每个 seed 均漏掉 4/5 个 crossing，15 个 seed-pair 事件中共漏掉 12 个”。

#### AUD-02

- 状态 `provenance_gap`
- 位置 `crossing_errors_temp_extrap.csv` 与 Table 6
- 检查 crossing table 的 meta 生成时间和输入清单
- 结果 crossing table 生成于 learned QH-residual 和 tlog 汇入主表之前，Table 6 因而不含当前 learned predictors，虽然主表已经报告它们的 aggregate Tc 指标。
- 修正 重新生成 crossing table，并把 learned rows 的 per-seed 和 per-pair记录纳入同一 provenance chain。

#### AUD-03

- 状态 `confirmed_internal_error`
- 位置 notebook export lines 1-13 与 878-927
- 检查 “只展示 CSV，不重算实验”的声明
- 结果 Section C 直接从 reference_G 和 fqh 计算新的 entropy 和 Tc-shift 指标，不是读取 producer CSV。
- 修正 将该计算移入 versioned producer，或修正文档声明并将所有输入纳入 manifest。

#### AUD-04

- 状态 `confirmed_internal_error`
- 位置 notebook export lines 680-704 与 709-929
- 检查 export 顺序
- 结果 workbook、Markdown 和 manifest 在 Sections A-C 之前写出。`EXTRA_SHEETS` 在文件末尾定义后没有再次 export，因此新增诊断不在交付包中。
- 修正 把 export 放到最后，建立单一 `SHEETS` 集合并对所有输出做存在性和 hash 检查。

#### AUD-05

- 状态 `provenance_gap`
- 位置 current working tree 的 `result/tables/*.meta.json`
- 检查 meta 中声明的 SHA-256 与当前输入文件
- 结果 crossing error 的 4 个 meta、metric definitions、phase inventory、split definitions 和 system inventory 均存在 input hash mismatch。主 predictor comparison meta 未发现同类 mismatch。
- 修正 在 clean tree 上运行完整 reproduction，禁止手工改动 producer inputs，并让 CI 验证所有 meta hash。

#### AUD-06

- 状态 `aggregation_ambiguity`
- 位置 predictor comparison 的 `Tc_error_K`、`missed_crossings`、`pairs_covered`
- 检查 QH-residual 的 aggregation
- 结果 102.379 K 只来自每个 seed 唯一成功的 SiO2 crossing，共 3 个数字，而 missed crossing 为 12。`pairs_covered=5` 是被评估的 pair，不是成功预测的 pair。
- 修正 输出 `numeric_crossings / reference_crossings`，对 conditional Tc error 和 coverage 明确分栏。

#### AUD-07

- 状态 `unresolved_input_needed`
- 位置 QH-residual seed dispersion
- 检查三个 seed 的 aggregate MAE 和 checkpoint hash
- 结果 checkpoint hash 不同，但 aggregate MAE 几乎相同，seed std 约 5.4e-9 eV/atom。没有 checkpoint 文件和完整 log，不能判定是确定性收敛、强约束、错误复用或汇总问题。
- 修正 保存 checkpoint、训练日志、初始化摘要和 pair-level seed dispersion，再做结论。

#### AUD-08

- 状态 `passed`
- 检查 skill arithmetic 与 SiO2 slope conversion
- 结果 QH-only 和 QH-residual skill 与 `1-MAE/MAE_zero` 一致。7.449 至 7.551 micro-eV atom^-1 K^-1 对应 1 meV 约 132 至 134 K。

#### AUD-09

- 状态 `reproducibility_blocker`
- 位置 tracked notebook、下载的 export、repository environment
- 检查 notebook execution、依赖、测试与 checkpoint
- 结果 tracked notebook 没有保存输出且与 929-line export 不同步。本地环境缺少 pandas 与 pytest，checkpoint 不在本地，本次只能做语法检查和静态审计。
- 修正 发布锁定环境和独立 rerun log，在解释新增实验前先通过 headline-number reproduction gate。

### Minor revision checklist

- 为 QH-only、QH-residual、representation-only 和 thermo-form 增加统一模型表，列出输入、可训练项、参数量、calibration、训练标签和适用域。
- 统一 energy unit，并在表中直接给出 slope-to-K conversion。
- 区分 frame、phase、pair、system、fold、crossing 和 seed 这些统计单位。
- 解释 1% imaginary-mode threshold 的来源并做敏感性分析。
- 使用明确的状态词。建议只用 frozen、provisional、not trained、structural N/A、artifact unavailable。
- 修复 README Quick start 的 code fence。当前说明性 prose 位于 bash block 内。

### Broad-interest and significance readout

正向方法论文的 broad-interest case 目前很弱。负结果与 benchmark 论文的潜力更大。

若项目能证明“物理先验并不天然改善 OOD 泛化，适用域错配会使结果恶化，而受约束的 thermodynamic form 在小数据跨温度任务中更可靠”，它会形成机器学习、统计力学与材料模拟之间更普遍的结论。这个方向比继续维护原标题更有说服力。

## Risk / unsupported claims

- 原标题为不支持主张。
- 跨 phase 和 system 的 learned extrapolation 尚未评估。
- QH residual 的因果优势未建立。
- thermo-form 的公平性与泄漏边界尚未完全审计。
- 对 entropy mechanism coverage 的分类没有完整分解证据。
- 对 seed robustness、functional robustness 和 truth robustness 的结论均不可作出。
- 当前 repository 状态不能支撑“已形成完整复现包”的强表述。

## 建议的重构标题

首选方向

> Beyond interpolation: a benchmark reveals when quasi-harmonic priors help and fail in free-energy extrapolation

更保守方向

> Stress-testing physics-informed free-energy models across temperature, phases and material systems

如果 thermo-form 公平复核后仍保持显著优势

> Constrained thermodynamic forms outperform learned residuals for out-of-domain phase free energies

## 下一步

只做一个最关键的闭环。

1. 冻结 clean commit 和 environment，重新生成所有 table meta 与 crossing table。
2. 完成 temperature split 上严格匹配的四模型消融，并预先写明主要终点是 per-pair ΔG、crossing coverage 和 conditional Tc error 的联合判据。
3. 先复核 thermo-form 无泄漏和同等自由度边界。
4. 若 QH-residual 仍无稳健优势，立即转向“when physics priors help and fail”的负结果主线，再开展 LOPO 和 LOSO。
5. 只有 temperature 轴闭环后，才扩展 learned LOPO、family LOSO 和 physics-regime stress tests。

当前状态是 `proposal_score 5.75/10`，未通过 7.0 的 supervisor-facing gate。剩余最大风险是继续围绕一个已被当前结果削弱的正向标题补实验，从而让研究范围扩大但中心论点更弱。
