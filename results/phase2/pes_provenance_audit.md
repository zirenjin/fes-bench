# PES provenance audit (P0-1)

审计日期：2026-09-20。审计对象是远端
`/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/`（用户约定的 `/share/jzr/`
工作区在该机器上解析到同一实验树）。本审计逐个读取实际 JSON/config 和训练日志；不以
方案文档中的描述替代配置证据。

## 结论

实验树中确实存在一条独立的 SiO₂ PES 微调 run，但它不是 FES head 训练：
`pes_dpa31_fes_full_singlehead_a100` 的 descriptor 和 energy fitting net 都可训练，loss
是 `ener`，使用 DFT energy/force/virial，训练 5000 步。相反，`fes_*`、
`results/continuous_sio2/*` 和 `hf_*_fes_*` 配置是冻结 DPA descriptor 加可训练 FES
head，loss 是 `fes`，不能拿来作为 PES 微调证据。

## 逐类配置证据

| 实验/配置 | 实际模型与 loss | 数据/步数 | 势来源判定 |
|---|---|---|---|
| `pes_dpa31_fes_full_singlehead_a100/input_v2_compat.json`（同内容的根目录 `input_pes_dpa31_fes_full_singlehead_a100.json`） | `descriptor.type=dpa3`, `descriptor.trainable=true`; `fitting_net.type=ener`, `fitting_net.trainable=true`; `loss.type=ener`, `pref_e=0.02→1`, `pref_f=1000→1`, `pref_v=0` | 56 个 DeepMD system；`training.numb_steps=5000`；数据根为 `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/pes/...` | **微调 PES** |
| `fes_dpa31_fullpes_singlehead_big_a100/input_v2_compat.json`（以及同前缀的 30+ 个变体） | `descriptor.trainable=false`; `fitting_net.type=fes`, `property_name=free_energy`, head trainable; `loss.type=fes` | `/fes/train/quartz_beta` + `/fes/train/cristobalite_beta`; `numb_steps=3000` | **FES head；不是 PES** |
| `results/continuous_sio2/{polynomial,tlog_polynomial}/seed*/config.json` | descriptor `trainable=false`; FES fitting net；`loss.type=fes` | 每个 seed `numb_steps=3000`；训练输出写入各自 FES run 目录 | **冻结 pretrained DPA + FES head；不是 PES** |
| `hf_fes_fewshot*`, `hf_tlog_fes_fewshot*`, `hf_property_fewshot*` | 实际配置均为冻结 descriptor、FES/property head；没有 energy/force PES loss | 各自 few-shot 数据和 3000-step FES/property 训练 | **冻结 pretrained DPA + 表征/FES head；不是 PES** |
| `fes_static_only_4phase*`、`fes_static_only_hf_hcp_bcc_pbe*` | `static_protocol.json` 只有固定结构、温度/压力输入和 source；没有训练配置或 PES 标签 | 1649 个 SiO₂ 温度点；Hf 1252 个温度点 | **静态推理 provenance；不能据此声称微调 PES** |

## SiO₂ PES 微调的实际收敛证据

读取 `pes_dpa31_fes_full_singlehead_a100/train.log`：

* step 2200 的 validation `rmse_f=1.76e-02 eV/Å` 是该日志中最低值；
* step 5000 的 validation `rmse_f=1.95e-02 eV/Å`，`rmse_e=7.52e-04 eV/atom`；
* checkpoint 实际保存为 `model.ckpt-3000.pt`, `model.ckpt-4000.pt`,
  `model.ckpt-5000.pt`，训练日志明确写出 5000 步完成。

因此本轮声子验证应使用这个单头 PES checkpoint（`head=None`），不能使用
`Domains_Alloy`/`Domains_SSE_PBE` FES 表征 head，也不能使用任何 `fes_*` checkpoint。

## LOPO 记录的核对

LOPO/continuous SiO₂ 配置逐个显示 frozen descriptor，且 FES loss/head 可训练；该记录与
实际配置一致。`pes_dpa31_fes_full_singlehead_a100` 是另一个 energy/force 微调实验，不能
被并入 LOPO 的 frozen-descriptor 表格；若旧表把它当作 frozen representation，应更正为
“独立 PES fine-tune，不属于 FES LOPO”。

## 严格排除项

本审计及后续 P0-3 训练只允许使用 DFT energy/force/virial PES 标签。任何 TI/free-energy
表、`fes/train` 数据或 FES head 都不进入 PES 训练；若发现配置混入其中，run 立即作废。

## Deviations from design

* 方案中的工作区别名 `/shared/jzr/fes_experiment_pbe_d3bj/` 在 thu-GenSi 上不存在，实际
  路径为 `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/`；只做路径解析，没有复制或
  改写别人的 run。
* 旧 run 的 PES loss 将 `pref_v` 设为 0；因此本审计如实记录为 energy/force 主导的
  PES 微调，不把它夸大为使用有效 virial 权重的训练。
