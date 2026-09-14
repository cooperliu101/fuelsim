# 四项未通过比较的独立诊断

2026-09-14 后续更新：用户已授权调整两个正式 B14 例题的体导热系数，
两程序采用相同物性并重新通过比较，见 [B14 当前记录](b14_conductivity_update/README.md)。
下文原设置下的失败状态和探测结果保留为历史证据。

2026-09-14 用户已授权实施 B5.23 和 B11.4 的分项容差，见
[当前两项容差约定](two_case_tolerance_qualification.md)。两个 B14 保持原设置。

2026-09-14 按用户要求先探测两个 B14，再分解 B5.23，最后单独核查 B11.4。
本次没有修改生产算法、正式输入卡、正式原生参考或验收门槛。四项正式比较仍未通过。
新增 17 个 B14 原生探测与两个成功的 B11.4 均匀温度探测，均使用 Windows
Abaqus 2025 RELr427，每个任务一个计算核心，并发合计不超过八核。

## B14：体导热量级相关的接触热导截断

原输入的体导热系数为 1e-12 W/(m K)，接触热导为 1000 W/(m² K)，
所有节点的温度和位移均被规定。用相同规定状态改变体导热系数、接触热导、
热导表格式、接触迭代设置，得到以下结果。表中有效接触热导取同一个非零热流
样本，并由两个较大体导热系数的原生运行分离线性体导热贡献后计算。

| 体导热系数 W/(m K) | 指定接触热导 W/(m² K) | 测得有效接触热导 W/(m² K) |
|---|---|---|
| 1e-12 | 1000 | 0.00158870255958 |
| 1e-12 | 500 | 0.00158870255958 |
| 1e-9 | 1000 | 1.58870255958 |
| 1e-12 | 0.0001 | 0.0001 |
| 1e-6 | 1000 | 1000 |
| 1e-6 | 500 | 500 |
| 1e-3 | 1000 | 1000 |

这些结果在本组模型上符合 `h_effective=min(h_requested, C*k_bulk)`，
其中测得 `C≈1.58870256e9 /m`。这个 C 只是当前网格和设置的识别结果，
没有进入 Fuelsim，也不声称适用于任意网格或 Abaqus 接触类型。
在原输入量级下，指定热导 1000 被截到约 0.0015887，因此热反力相差约 63 万倍。
体导热系数提高到 1e-6 后，CAX4T/CAX8T 最大原生热反力分别恢复为
12.8955918738 W 和 12.8955901076 W，与 Fuelsim 原结果处于相同量级。
这两个探测改变了体导热物性，因此没有替换正式参考。

以下排除性检查也已完成：

- 增加完整压力相关热导表，交换压力表与间隙表顺序，或只保留压力表，均未消除原量级截断。
- 强制不将严重不连续迭代直接转为平衡迭代，未改变结果。
- 释放一个非界面温度自由度并收紧热平衡控制，未改变异常接触热反力的量级。
- 在较高体导热量级，指定接触热导减半，分离出的接触热流也减半。
- 较高体导热量级下，保留原间隙表与使用压力表得到相同的接触热反力。

[官方热接触说明](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEITNRefMap/simaitn-c-thermalinteraction.htm)
区分开口与闭合热导定义，因此最初检查了热导表格式；上述实测排除了它作为本次
差异的单独解释。实际截断是否有可公开设置的解除入口、其内部 C 的计算规则，
尚未鉴定。不能把实测截断规则作为物理定律加入 Fuelsim。

完整字面输入、原生节点/材料点/接触 CSV、完成记录和脚本在
[contact_thermal_diagnosis](contact_thermal_diagnosis/)。
[conductance_probe.tsv](contact_thermal_diagnosis/conductance_probe.tsv) 保存定量结果，
[analyze.py](contact_thermal_diagnosis/analyze.py) 仅对保存的数据作代数分解，不求解模型。

## B5.23：固定初始质量与原生热容离散不等同

从两个程序实际输出的节点温度、位移及时间历史，独立重构八个外边界节点、
20 个时间增量的热容和导热贡献，共 160 个样本。这些节点不在接触界面，
重构不需要拟合接触热流。密度为初始值 100 kg/m³，比热使用当前节点温度。

Fuelsim 容量权重为 `rho_initial*cp(T_i)*w_reference_i`。
Abaqus 2025 对应权重为上述值再乘 `J_element/J_i`，其中 J_element 是整体
当前/参考体积比，J_i 是配对 Gauss 点的当前/参考体积比。
在非均匀变形下，这两个比值一般不同。

| 检查 | 最大绝对误差 |
|---|---|
| 用识别的原生规则重构 Abaqus 热反力 | 2.5921e-11 W |
| 用固定初始质量规则重构 Fuelsim 热反力 | 1.0005e-11 W |

最严重样本为 t=0.4 s、节点 24：

| 分量 | 数值 |
|---|---|
| Abaqus 原生热反力 | 4828.88685030 W |
| Fuelsim 热反力 | 4793.70648378 W |
| 总差值 | -35.18036652 W |
| 在同一原生状态下，仅切换热容规则的差值 | -35.28257278 W |
| 温度、位移及历史状态差异的贡献 | +0.10220626 W |

两个分量之和重现总差值；结果支持维持用户指定的固定初始质量实现。
本次没有把原生非均匀热容权重移入生产程序，也没有降低 0.5% 的热反力门槛。
完整分解见 [b523_decomposition.tsv](thermal_residual_diagnosis/b523_decomposition.tsv)，
可用 [analyze_b523.py](thermal_residual_diagnosis/analyze_b523.py) 从保存的参考与生产输出复现。

## B11.4：均匀温度场中的极小热反力

正式输入所有温度节点始终为 600 K；没有温度梯度或温度变化，理论热反力为零。
另外只改变均匀温度基准为 300 K 和 1200 K，保留机械加载，重新计算原生响应。
每组比较覆盖 14 个时刻、294 个节点样本。

| 程序和均匀温度 | 最大节点热反力绝对值 W | 最大时刻总热反力绝对值 W |
|---|---|---|
| Abaqus，300 K | 7.15305008e-12 | 5.34e-28 |
| Abaqus，600 K | 1.43061002e-11 | 1.07e-27 |
| Abaqus，1200 K | 2.86122003e-11 | 2.14e-27 |
| Fuelsim，600 K | 5.82596616e-15 | 1.97e-31 |

原生极小热反力严格随温度基准缩放，且整体相互抵消；这些证据支持舍入误差解释，
而不是物理传热或应力误差。正式 600 K 比较差值仍为 1.43067424e-11 W，
所以仍不满足原来的 1e-11 W 绝对门槛。
结果见 [b114_uniform_temperature.tsv](thermal_residual_diagnosis/b114_uniform_temperature.tsv)。
第一次探测因未替换 ALL 节点集温度边界且提取器不匹配而无效；错误输入保存在
`thermal_residual_diagnosis/failed_attempt_inputs/`，对应失败日志在 `roundoff_logs/`。
成功重算使用正确的 B11.4 提取器，其记录在 `roundoff_confirmed_logs/`。

## 后续建议，尚未实施

1. B14 可以将两个程序的体导热系数同时设为 1e-6 W/(m K)，避免原生热导截断。
   这会改变例题物性；非界面微小体热流应改为与原生值比较绝对误差，保留
   1e-8 W 的绝对误差门槛，而不再断言热流自身小于该值。界面及机械比较门槛不变。
2. B5.23 保留固定初始质量，并可单独讨论将该例题节点热反力最大逐点相对门槛由 0.5%
   调整为 1%；0.05 W 绝对条件、整体相对二范数与峰值门槛、其他正式量门槛和所有样本覆盖保持不变。
3. B11.4 可单独将本例热反力绝对门槛由 1e-11 W 调整为 2e-11 W；不扩展到
   其他例题或其他物理量。

以上是基于诊断提出的新调整，不属于此前已授权的八项容差修改，尚未采用。

## 复现

在 Windows PowerShell 中分别运行 `contact_thermal_diagnosis/run_all.ps1`，
用 `-TaskFile` 指定同目录的 `tasks.json`、`pressure_tasks.json`、`scale_tasks.json`、
`iteration_tasks.json`、`validation_tasks.json` 或 `roundoff_tasks.json`。
每批最多八个单核任务，各批应顺序执行。使用不同 `-LogDirectory` 保留各批状态。
所有输入均为完整保存的 `.inp` 文件，运行程序不生成或改写问题定义。

在仓库根目录，三个分析入口为：

```bash
python verification/abaqus/2025/contact_thermal_diagnosis/analyze.py
python verification/abaqus/2025/thermal_residual_diagnosis/analyze_b523.py
python verification/abaqus/2025/thermal_residual_diagnosis/analyze_b114.py
```

后两项还读取 `build/blackbox` 中实际生产程序生成的原始 Exodus 输出。
相关正式测试仍使用 `ctest --test-dir build -j8 --output-on-failure`；
本次诊断没有将原生探测加入日常 CTest。

本次按八个并发额度重新运行四项正式比较及两个元数据检查。四项比较仍失败，
参考 SHA256 和验证矩阵检查通过，耗时 5.37 秒；结果见
[相关测试复查日志](thermal_residual_diagnosis/ctest.log)。三个独立分析脚本均成功运行。
