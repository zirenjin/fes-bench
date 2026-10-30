# S3–S5 证据与资源审计

审计日期：2026-10-27。结果只复用冻结的既有 C readout，或记录明确的阻断原因。没有把历史温度外推划分重命名为 LOPO。

## S3：统一 cristobalite-beta LOPO

C-simple 与 C-DPA 的同一 phase-LOPO readout 来自冻结 R 归档；其折叠、选择、表示和逐点预测哈希可在 `INPUT_HASHES.txt` 与原归档 provenance 中核验。S2 独立复算并与 V2 scorer 对照，验证了二者在目标相 LOPO 上的点差值和完整支持集指标。

T2 的历史定义重复加入静态项：`F_QH` 已包含 `E_DPA`，原式却为 `E_DPA + F_QH`。旧输出不能代表修正后的 T2。修正定义应为 `F_QH + r_theta + c_system`，其中 `r_theta` 为残差。T3 的现有 SiO2 运行属于温度尾部诊断；旧结构/元素类型映射运行已判无效，均不满足同一 LOPO。当前本地没有 DeepMD runtime，也没有可复用的目标标签隔离 LOPO checkpoint；到 Expanse 的 SSH 连接失败，因此没有提交训练任务。Delta AI 未用于计算。

R 归档的 source-only NPZ 对应 Ti LOSO，不是 SiO2 LOPO。SiO2 目标对没有合法 source-only readout 文件；本轮标记为不可用，避免用目标参考曲线补建系数。

## S4：静态项与热项分解

同源四相表可以用于其声明范围内的有限温热力学基线，但目前没有逐行映射的同相结构/静态能、共同零点/零点振动约定和固定相身份分支。故不能把 `G(T)-H0` 解释为受控的纯热自由能贡献，也不报告这类数值分解。

## S5：振动先验与消融

现有可靠性审计仅将 quartz 标为可用；beta-cristobalite 和 P63/mmc tridymite 不通过。C2221 修正版是固定晶胞谐振自由能，含不稳定 Γ 模且缺少相同参考层静态项，不能作为完整 QHA。由于不存在覆盖统一同一 LOPO 目标相的共同适用子集，本轮不做 M0–M3 训练比较。

详见 `data/EXPERIMENT_STATUS_DETAILED.csv`、`data/S1_RANK_SVD_AND_C2221_ADMISSION.json` 与原始 P/R 证据清单。
