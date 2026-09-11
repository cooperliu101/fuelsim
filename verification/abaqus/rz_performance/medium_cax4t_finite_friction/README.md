# CAX4T 中等规模有限应变摩擦对比

在现有中等规模有限应变 CAX4T 算例上增加机械摩擦。生产程序仍为 `22df663`，
没有修改求解器或单元代码，没有改变网格、材料、热源、时间步和收敛容差。
本例作为手动验证资料保存，没有加入 CTest。

## 物理与数值设置

- 保留原来的 7424 个 CAX4T 单元、7670 个节点和 23010 个求解自由度。
- 保留有限应变、节点到面的有限滑移接触、法向罚刚度 `1e14 Pa/m`，
  间隙热导率 `0.4 W/(m K)` 和最小热间隙 `1e-6 m`。
- 增加库仑摩擦系数 `mu=0.2`、弹性滑移容许比例 `slip_tolerance=0.001`。
  参数沿用已有 B8.1 有限滑移摩擦验证。Abaqus 使用相同参数，并显式设置
  `*Gap Heat Generation` 为零；Fuelsim 本例不将摩擦耗散转化为热源。
- 保留 20 个单位增量，将芯块体热源线性加载至 `2e8 W/m^3`，包壳外侧
  温度为 600 K，芯块和包壳底部轴向位移为零。每个增量求稳态热平衡，
  不含热容项；时间步之间提交有限应变和摩擦历史。
- 精度输入和计时输入的物理定义完全相同，计时输入关闭结果文件输出。
  两个程序均使用单处理器直接求解，Fuelsim 使用 MUMPS。

## 精度结果

比较所有末态节点、材料积分点和接触节点，保持原有相对 L2、相对绝对峰值和
最大逐点相对误差均小于 `0.01%` 的门槛。零参考量使用独立绝对误差约束，
没有增加相对误差分母下限。新增的摩擦牵引、切向力和滑移按符号比较。

| 字段 | 最大逐点相对误差 |
|---|---:|
| 温度 | 0.000000267% |
| 自由位移向量 | 0.000003052% |
| 应力张量 | 0.000004321% |
| 接触压力 | 0.000003794% |
| 接触间隙 | 0.000003803% |
| 法向接触力 | 0.000003795% |
| 摩擦牵引 | 0.000042464% |
| 累计切向滑移 | 0.000045481% |
| 切向接触力 | 0.000042462% |

全部指标通过，原始指标保存在 `accuracy.json`。Fuelsim 输出包含初始状态及
20 个增量，Abaqus 本次只保存末态场；因此上述为末态全场对比，不声称完成了
两个程序的逐增量全场比较。

Fuelsim 末态 65 个接触节点全部接触，61 个滑动、4 个粘着，带符号切向力合计
`138.145184 N`，最大摩擦牵引为 `1.824022 MPa`，最大累计切向滑移为
`2.946060 um`。库仑限制的最大超出量仅为双精度舍入量 `2.33e-10 Pa`。
`friction_diagnostics.json` 还保存每个增量的接触、滑动数量和切向力。
这些数量是 Fuelsim 的历史诊断，末态摩擦场另外通过了 Abaqus 精度检查。

## 正式耗时

固定逻辑 CPU 0，所有线程环境变量为 1，各次求解串行执行。
每个程序先预热一次，再正式运行两次；计时阶段没有同时编译或运行测试。
外部时间包含程序启动和求解，内部时间另列，不混用两种计时口径。

| 程序 | 第一次正式外部耗时 | 第二次正式外部耗时 | 外部均值 | 非线性迭代次数 |
|---|---:|---:|---:|---:|
| Fuelsim | 36.3522 s | 36.5685 s | **36.4603 s** | 47 |
| Abaqus | 45.6105 s | 43.6428 s | **44.6266 s** | 54 |

Fuelsim 本次外部用时少 **18.30%**。两个程序均完成 20 个增量，Fuelsim
没有时间步重试。Fuelsim 内部时间分别为 35.7978 和 35.8432 秒，Abaqus
`JOB TIME SUMMARY` 的分析时间分别为 41 和 40 秒。

环境仍为同一主机上的 WSL2 Fuelsim 与 Windows Abaqus；跨虚拟机边界固定逻辑
处理器不能证明落在同一物理核心。两次正式重复只支持本算例、本次环境下的结论。
原无摩擦计时没有在本次重新运行，不能据不同批次结果精确归因摩擦带来的耗时增量。

## 验证与复现

本次仅增加输入和扩展只读对比、运行脚本，生产可执行文件未变。
CAX4T 单元测试、B8.1 有限滑移摩擦对标、摩擦历史契约测试共 3 项通过；
原无摩擦精度报告重新计算后各项指标完全不变。

Fuelsim 精度入口：

```bash
build/fuelsim -i verification/fuelsim/quasistatic_rz_performance_medium_cax4t_finite_friction.fsi
```

Abaqus 在 Windows PowerShell 中运行（源目录填写本仓库对应的 Windows/WSL 路径）：

```powershell
.\run.ps1 -SourceDirectory <source> -Size medium -Element cax4t -Strain finite -Friction -Runs 1 -ResultsDirectory <accuracy>
.\run.ps1 -SourceDirectory <source> -Size medium -Element cax4t -Strain finite -Friction -Timing -Runs 3 -ResultsDirectory <timing>
```

Fuelsim 正式计时：

```bash
/home/cooper/miniforge/envs/moose/bin/python \
  verification/abaqus/rz_performance/medium_cax4t_finite_friction/run_timing.py \
  --root "$PWD" --executable "$PWD/build/fuelsim" --results /tmp/cax4t-friction-repeat
```

精度检查：

```bash
env NETCDF_LIBRARY=/home/cooper/miniforge/envs/moose/lib/libnetcdf.so \
  /home/cooper/miniforge/envs/moose/bin/python verification/abaqus/rz_performance/compare.py \
  <fuelsim-result.e> \
  verification/abaqus/rz_performance/medium_cax4t_finite_friction/rz_performance_medium_cax4t_finite_friction \
  --element cax4t --incremental --friction --report /tmp/cax4t-friction-accuracy.json
```

`provenance.json` 保存可执行文件、输入、网格、脚本和结果摘要；三个压缩 CSV
保存 Abaqus 参考，压缩日志保存原始计算和计时证据。`summary.json` 保存正式时间、
迭代次数、精度结论和适用边界。
