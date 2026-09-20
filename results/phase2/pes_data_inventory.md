# PES data inventory (P0-2)

审计日期：2026-09-20；只盘点可用于 PES 微调的 DFT energy/force/virial 类标签，不把
TI/FEP/free-energy 表当作训练数据。

## 总表

| 体系 | 盘点路径 | 可用 PES 数据 | 规模/覆盖 | 泛函一致性 | `pes_finetunable` |
|---|---|---|---|---|---|
| SiO₂ | `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/pes/` | **有** | 56 个 DeepMD systems、502 frames；quartz_beta 24 systems，cristobalite_beta 24 systems，tridymite_beta 8 systems（全部 `ood_tridymite`）；每个 raw set 有 `energy.npy`, `force.npy`, `virial.npy`, `coord.npy`, `box.npy` | 该实验树命名和 reference manifest 为 PBE-D3(BJ)，但 raw DeepMD 目录本身没有 functional 字段；因此标记为“来源口径 PBE-D3(BJ)，raw file 内未独立再验证” | `true` |
| Hf | `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/external/darus_ti_zr_hf/Hf_hcp_PBE.tar.gz`、`Hf_bcc_PBE.tar.gz` | **有，但标签不完整** | 每个压缩包共 725 个 OUTCAR，其中 `c.training_set_by_high_DFT` 可解析构型各 721 个（hcp 721、bcc 721）；README 明确分为 low/high DFT training sets；sample OUTCAR 含 `TOTEN` 和 `TOTAL-FORCE (eV/Angst)`；当前 OUTCAR `ISIF=0`，未发现可直接解析的 stress/virial 行 | 归档名和 README 明确为 PBE；与 Hf reference 的 PBE 口径一致 | `true`（energy/force；virial 缺口需在微调记录中明确） |
| Ti | `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/external/darus_ti_zr_hf/Ti_hcp.tab`, `Ti_bcc.tab` | **无** | 仅热力学表；列为 `T, G, H, S, C_p, C_V, V, α, B_T, B_S`，没有构型、energy/force/virial 张量 | 表头/文件没有 DFT label provenance；不能把 PBE 热力学表当 PES 标签 | `false` |
| Zr | `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/external/darus_ti_zr_hf/Zr_hcp.tab`, `Zr_bcc.tab` | **无** | 仅热力学表；列同 Ti，没有构型、energy/force/virial 张量 | 表头/文件没有 DFT label provenance；不能把 PBE 热力学表当 PES 标签 | `false` |

## SiO₂ 分割与帧计数

按 `energy.npy` 目录计数得到：`train=16`, `valid=16`, `test=16`,
`ood_tridymite=8`；合计 56 systems。每套 raw set 的 energy/force/virial/coord 数组均为
3 帧，所以总帧数为 502。训练配置实际列出同一批 56 systems，未使用 FES/TI 表。

## Hf 标签边界

Hf 压缩包的 725 个 OUTCAR 计数来自 `tar -tzf ... | grep 'OUTCAR$' | wc -l`；按
`c.training_set_by_high_DFT` 过滤并用 ASE 读取后实际构建 hcp/bcc 各 721 个构型。README
说明这些是 low/high DFT 训练集；抽样 OUTCAR 的 `TOTEN` 和 `TOTAL-FORCE` 可见，但其
`ISIF = 0` 表明这批归档不能直接声称包含 stress/virial 标签。Hf 的 PES 微调若继续，
应采用 energy+force loss（`pref_v=0`）或先补充独立 virial 数据；不得把 Hf `.tab` 的
热力学量伪装成 virial。

## 体系状态字段

已在以下文件写入 `pes_finetunable`：

* `data/sio2/system.json`: `true`
* `data/hf/system.json`: `true`（energy/force 可微调，virial 缺口见上）
* `data/ti/system.json`: `false`
* `data/zr/system.json`: `false`

## Deviations from design

* 设计表将 Hf/Ti/Zr 统一指向同一 DaRUS 数据集，但实际副本中 Ti/Zr 只有 thermo `.tab`，
  Hf 才有 OUTCAR 归档；按实际文件内容拆开记录，没有用表头推断不存在的 PES 数据。
* SiO₂ raw DeepMD 数据没有在目录内重复写 functional 元数据，因此没有把实验目录名当作
  独立验证；functional 一致性在 inventory 中明确标为来源口径而非 raw-file 自证。
