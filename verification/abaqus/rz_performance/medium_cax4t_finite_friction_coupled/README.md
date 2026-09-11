# CAX4T 中等规模：瞬态热、有限应变、摩擦与蠕变塑性耦合

本例在已有 CAX4T 中等规模有限应变摩擦模型中，为包壳加入 Norton 蠕变与
J2 线性各向同性硬化塑性。生产可执行文件未修改，验证资料没有加入 CTest。

## 物理定义和加载口径

保留 7424 个单元、7670 个节点、23010 个求解自由度，以及原有弹性、热膨胀、
导热、法向罚刚度 `1e14 Pa/m`、有限滑移节点到面接触、摩擦系数 `0.2` 和
滑移容许比例 `0.001`。芯块仍为热弹性材料。

包壳非弹性参数沿用 B13 综合验证：

- Norton 等效蠕变率为 `1e-5 * (q / 5e6)^3 s^-1`，`q` 为 Pa 制等效应力。
  Abaqus 对应 `*Creep, law=TIME` 的系数为 `8e-26`、应力指数为 3、时间指数为 0。
- 初始屈服应力为 `5e6 Pa`，线性各向同性硬化模量为 `2e9 Pa`。
- 不计摩擦、塑性和蠕变耗散生热；这些耗散与外加体热源分别处理。

必须同步调整参考的热分析步骤：原 Abaqus `STEADY STATE` 步骤明确报告
`CREEP AND SWELLING EFFECTS ARE OMITTED IN THIS STEP`，无法用于蠕变对比。
双方因此启用真正的瞬态热容量，使用原输入中的密度和比热：芯块
`10970 kg/m^3`、`300 J/(kg K)`，包壳 `6500 kg/m^3`、`330 J/(kg K)`。
文件名前缀 `quasistatic` 指机械上忽略惯性，本例热过程为瞬态。

功率曲线显式定义为 0～20 秒从零线性上升至 `2e8 W/m^3`，双方每个增量使用
相同的区间平均体热源。Fuelsim 使用已有的 `TimeFunctions` 与
`heat_source_time_evaluation=interval_average`；Abaqus 使用 `HEAT_RAMP` 幅值。
第 1 秒的平均体热源为 `5e6 W/m^3`，第 20 秒为 `1.95e8 W/m^3`，
最终曲线端点为 `2e8 W/m^3`。20 个固定 1 秒增量、边界约束及外侧 600 K 温度不变。

Abaqus 的步骤类型和未指定幅值时的加载规则见官方
[耦合温度—位移步骤](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-coupledtemperature-displacement.htm)和
[分布热流幅值](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-dflux.htm)。
原先稳态步骤遗漏蠕变以及转为瞬态后未显式指定功率幅值的尝试，作为失败参考诊断
保存在 `diagnostics/`，不作为生产实现误差或正式计时证据。

## 精度与同时激活

末态比较覆盖全部节点、全部材料积分点和 65 个接触节点。检查原来的温度、
位移、应力和摩擦场，并新增弹性、塑性、蠕变张量，以及两个等效非弹性标量。
末态所有指标均通过相对 L2、相对绝对峰值和最大逐点相对误差 `0.01%` 门槛；
零参考量独立检查绝对误差，不使用分母下限。

Abaqus 额外导出每个增量、每个包壳积分点的等效塑性应变、等效蠕变应变及
等效应力。用双方输出直接核实同一点、同一增量内的两种非弹性增量均大于
`1e-14`，该值只用于排除数值噪声的激活计数，不用于放宽精度分母。
从第 14 秒开始有共同激活的积分点，到第 19、20 秒，包壳全部 4096 个积分点
均在双方同时产生塑性和蠕变增量。

本算例的逐步蠕变增量与各自步末等效应力代入 Norton 后向 Euler 公式一致，
详见 `comparison.json` 中每步的 `fuelsim_integration` 和 `abaqus_integration`。
这一结论来自当前输出，不能推广成 Abaqus 在所有加载问题中都使用隐式蠕变积分。

额外的逐增量等效应变检查保留了一项严格门槛差异：第 18 秒，单元 7145、
Abaqus 积分点 3 的等效塑性应变，Fuelsim 为 `8.475406985e-9`，Abaqus 为
`8.470959280e-9`，绝对差 `4.447705047e-12`，逐点相对误差 `0.0525053%`。
该时刻塑性全场相对 L2 误差为 `0.0000126443%`。全部逐增量等效应变仍通过
项目 `0.5%` 门槛，但并非所有历史点都通过额外的 `0.01%` 目标。
检查器保留这项 `passed=false`，没有调整分母、材料参数或误差阈值。

外部比较范围为末态完整场，以及所有 20 个时刻的等效塑性、等效蠕变历史。
没有将两个标量历史的比较扩大称为全部张量的逐时刻比较。

## 性能和验证方法

只有最终物理、热源和材料历史口径对齐后才进行正式性能测量。
双方固定逻辑 CPU 0、单线程、直接求解，Fuelsim 使用 MUMPS。
每套程序一次预热、两次正式运行，关闭结果文件输出，各次串行执行，
不与编译或回归同时运行。外部耗时与内部求解时间分别记录。

| 指标 | Fuelsim | Abaqus |
|---|---:|---:|
| 正式外部耗时，第 1 次 | 46.4613 s | 57.8471 s |
| 正式外部耗时，第 2 次 | 46.8574 s | 59.7397 s |
| 正式外部耗时平均值 | 46.6594 s | 58.7934 s |
| 每次非线性迭代数 | 58 | 71 |
| 每次完成增量数 | 20 | 20 |
| 末态单侧切向摩擦力合计 | 36.605970933 N | 36.605970829 N |

Fuelsim 本次外部耗时减少 `20.64%`。其正式内部求解时间分别为
45.9578、46.3311 秒；Abaqus `JOB TIME SUMMARY` 分析时间分别为 54、55 秒。
内部时间统计范围不同，不替代上表的外部计时。末态最大逐点相对误差为
`0.0000156323%`，出现在节点切向力。

同一主机上 Fuelsim 在 WSL2 运行，Abaqus 在 Windows 运行；逻辑 CPU 绑定
不能保证跨虚拟机边界落在相同物理核心。两次正式重复仅支持本算例、本次环境的结论。
此次热过程也发生了变化，不能将与此前稳态热弹性摩擦例题的时间差完全归因于材料耦合。

相关 CAX4T 局部切线、有限应变蠕变塑性、蠕变积分、摩擦和历史契约共 5 项测试通过。
本次只增加完整输入与只读验证资料，没有修改生产代码。

## 复现

Fuelsim 精度计算：

```bash
build/fuelsim -i verification/fuelsim/quasistatic_rz_performance_medium_cax4t_finite_friction_coupled.fsi
```

Abaqus 在上级目录的 Windows PowerShell 中运行：

```powershell
.\run.ps1 -SourceDirectory <source> -Size medium -Element cax4t -Strain finite -Friction -Coupled -Runs 1 -ResultsDirectory <accuracy>
.\run.ps1 -SourceDirectory <source> -Size medium -Element cax4t -Strain finite -Friction -Coupled -Timing -Runs 3 -ResultsDirectory <timing>
```

检查末态和材料历史：

```bash
env NETCDF_LIBRARY=/home/cooper/miniforge/envs/moose/lib/libnetcdf.so \
  /home/cooper/miniforge/envs/moose/bin/python \
  verification/abaqus/rz_performance/medium_cax4t_finite_friction_coupled/check_coupled.py \
  <fuelsim-result.e> \
  verification/abaqus/rz_performance/medium_cax4t_finite_friction_coupled/rz_performance_medium_cax4t_finite_friction_coupled \
  --report /tmp/cax4t-coupled-check.json
```

Fuelsim 正式计时：

```bash
/home/cooper/miniforge/envs/moose/bin/python \
  verification/abaqus/rz_performance/medium_cax4t_finite_friction/run_timing.py \
  --root "$PWD" --executable "$PWD/build/fuelsim" --results /tmp/cax4t-coupled-timing \
  --element cax4t --coupled
```

`comparison.json` 保留全部末态指标、逐步等效应变指标、同时激活计数和积分诊断。
`summary.json` 汇总精度边界与计时，`provenance.json` 保存输入、可执行文件、
网格、参考和完整结果摘要。压缩 CSV、日志与失败参考诊断用于复查和复现。
