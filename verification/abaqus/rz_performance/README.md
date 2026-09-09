# 轴对称生产输入性能与 Abaqus 验证

原 `fuelsim_m1_single_core_benchmark` 已删除，现使用 `fuelsim -i` 运行完整输入卡。
本模型检查稳态热弹性接触，不包含塑性、蠕变或有限应变。

| 输入 | 芯块径向×轴向 | 包壳径向×轴向 | 单元数 | 节点数 | 自由度数 |
|---|---:|---:|---:|---:|---:|
| medium | 100×64 | 16×64 | 7,424 | 7,670 | 23,010 |
| large | 200×64 | 32×64 | 14,848 | 15,210 | 45,630 |

芯块半径 4.120 mm、高 10 mm；包壳内外半径 4.122/4.692 mm、高 10.020 mm。
初始气隙 2 μm，两个区域独立节点；网格在每个方向均匀划分。
芯块轴线径向位移为零，芯块和包壳底面轴向位移为零，包壳外表面温度为 600 K。
两区域初始温度均为 600 K。芯块体热源经 20 个等幅加载增量从零增至 2e8 W/m³。
其余未施加热边界条件的表面绝热；机械接触无摩擦，不施加外部压力。

保留原程序的几何、导热、弹性、热膨胀和接触参数。为与 Abaqus 公式匹配，
单元由默认 Quad4 改为 CAX4T，热接触显式采用节点到面，机械接触同样为节点到面。
两套程序使用相同的节点编号、连接关系、边集、20 个固定加载增量和直接线性求解。
Fuelsim 使用 MUMPS；Abaqus 使用其原生直接求解器，不能声称两者使用同一个求解器实现。

`../generate_rz_performance_meshes.py` 只生成网格和 Abaqus 网格包含文件。
四张完整生产输入卡位于 `../../fuelsim/steady_rz_performance_{medium,large}{,_timing}.fsi`。
`_timing` 卡仅移除结果输出配置；验证卡写出最终结果。
Abaqus 的 `.inp`、材料和接触表均已显式保存；温度相关导热率使用 0.2 K 间隔的
解析函数采样，实际温度范围必须处于表内。间隙导热表在初始间隙附近使用 1e-8 m 间隔。
表格不按比较误差拟合。模型密度和比热不参与两套稳态方程。

## 精度入口

```bash
build/fuelsim -i verification/fuelsim/steady_rz_performance_medium.fsi
NETCDF_LIBRARY=/home/cooper/miniforge/envs/moose/lib/libnetcdf.so python \
  verification/abaqus/rz_performance/compare.py \
  verification/fuelsim/steady_rz_performance_medium_results.e \
  verification/abaqus/rz_performance/rz_performance_medium \
  --report verification/abaqus/rz_performance/medium_accuracy.json
```

large 使用对应同名文件。检查程序只读取生产输出和 Abaqus 参考，不调用求解器。
它验证全部节点、全部四个材料积分点和全部接触节点的编号与数量。
温度、自由位移向量、应力张量、接触压力、间隙、法向力和总法向力分别检查相对
二范数、相对绝对峰值和最大逐点相对误差，当前统一要求小于 0.01%，不设置分母下限。
张量范数包含两倍剪切分量平方。两套程序的规定零位移分量分别按 1e-12 m 检查，
不把 Abaqus 规定零位移处的舍入残值用作相对误差分母。
该精度结论只覆盖最终加载状态；不代表网格收敛，也不比较未输出的中间增量场。
稳态生产输出不包含支承反力和边界热反力，这两项不包含在本次比较范围内。

## 计时入口和范围

```bash
python benchmarks/run_rz_performance.py \
  --windows-source '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance'
```

计时前停止其他编译和求解任务。每档每个程序预热一次，随后记录两次；顺序运行，
使用单进程、单线程，Fuelsim 固定到 WSL 逻辑 CPU 0，Abaqus 固定 Windows CPU 0。
WSL 与 Windows 的逻辑 CPU 编号不保证对应同一物理核心，因此这是同一主机上的
跨平台运行观察，不能解读为同一操作系统下的纯求解器性能差异。
Fuelsim 使用外部单调时钟计时，Abaqus 使用 PowerShell Stopwatch；两者都覆盖
程序启动、输入读取和完整求解。Abaqus 计时包含启动及许可证等开销。
另行保存 Fuelsim 内部求解时间和 Abaqus JOB TIME SUMMARY 分析时间，不混用。
两套程序均关闭场结果、历史结果和重启动文件输出，仍保留运行日志及 Abaqus 必需文件。

大网格继续作为手动验证资料，不加入轻量自动回归。

## 初始配置结果（2026-09-09，优化前）

两次正式测量的均值如下，预热运行不参与均值。

| 规模 | Fuelsim 外部总时间 | Abaqus 外部总时间 | Fuelsim / Abaqus | Fuelsim 内部求解时间 | Abaqus 分析时间 | 非线性迭代数（Fuelsim / Abaqus） |
|---|---:|---:|---:|---:|---:|---:|
| medium | 70.34 s | 41.80 s | 1.683 | 69.84 s | 37.00 s | 145 / 52 |
| large | 163.90 s | 83.13 s | 1.971 | 163.37 s | 78.50 s | 158 / 52 |

初始回溯线搜索配置下，两档 Fuelsim 都慢于 Abaqus。非线性迭代数也更多；本次记录这一事实，
没有修改生产求解器或通过放宽精度门槛取得更好的速度。启动与分析时间的定义及跨平台
边界见上文。原始样本见 `timing.json`，每次运行的 `.log`、`.dat.gz`、`.msg.gz` 和 `.sta.gz`
保存了对应的工作量与完成状态。

| 比较量 | medium 最大逐点相对误差 | large 最大逐点相对误差 |
|---|---:|---:|
| 温度 | 2.69319322e-07% | 2.66266004e-07% |
| 自由位移向量 | 3.6555955e-06% | 3.57305734e-06% |
| 应力张量 | 4.7116045e-06% | 4.63624859e-06% |
| 接触压力 | 4.06128122e-06% | 4.00893005e-06% |
| 接触间隙 | 4.06366757e-06% | 4.00405276e-06% |
| 逐节点法向力 | 4.06232467e-06% | 4.00994981e-06% |
| 总法向力 | 3.46598801e-06% | 3.40510667e-06% |

全部检查均通过 0.5% 门槛，规定零位移在两套程序中均通过 1e-12 m 绝对门槛。
温度、位移和应力检查覆盖 7,670/15,210 个节点与 29,696/59,392 个材料积分点；
两档均检查全部 65 个接触节点。完整三项相对指标和绝对误差见两个 `*_accuracy.json`。

移除旧性能目标后的 Release 统一回归：271/271 通过，164.48 秒。
比较器另用 1% 非零误差和超限零参考误差验证拒绝行为。大网格保持手动运行。

最终压缩提取脚本已重新读取两套 Abaqus 结果数据库，六张参考表解压后的内容与首次提取逐字节一致。
`SHA256SUMS` 校验本目录保存的输入、脚本、参考与测量证据；`provenance.json` 另记录网格、生产输入、可执行文件和构建环境。


## 中等规模求解优化（2026-09-09）

当前中等规模两张完整输入卡使用 `line_search = basic`，即优先使用完整牛顿更新；
`backtracking_fallback = true` 表示该次求解失败时，从同一个初值改用回溯线搜索重试。
本次全部 20 个加载增量直接收敛，没有触发回溯重试。大规模输入卡保持原设置。
网格、CAX4T 公式、材料、接触、20 个加载增量、MUMPS 和全部收敛容差均保持不变。
这是对已有求解功能的输入配置优化，没有修改生产数值核或求解器实现。

逐增量对照表明，原回溯线搜索强制残量下降，在第 7、13、14 个增量分别需要
23、37、18 次迭代；完整更新分别需要 3、5、4 次。第一增量的完整更新允许残量
从 0.07638 短暂上升到 0.07740，随后下降到 9.389e-6 和 2.306e-11。
最终总迭代数从 145 降到 65，残量计算从 262 次降到 85 次。
详细逐次残量保存在 `medium_optimized/iteration_diagnosis.json`。

PETSc 分阶段诊断还显示，矩阵压缩后重新进行 MUMPS 结构分析占有较大开销：
原配置累计约 24.5 秒，完整更新后约 10.9 秒。现有代码压缩掉精确为零的条目，
因此接触状态变化可能改变用于分解的矩阵结构，不能直接假定每次结构相同。
本次通过减少迭代次数降低这部分开销，没有改变结构复用规则。
分阶段诊断开启了监视输出，只用于解释耗时，不作为正式速度比较样本。

当前比较程序将三项相对误差门槛从 0.5% 收紧到 0.01%，保持零参考值的独立
绝对误差规则。优化后的最终场完整比较见 `medium_optimized/accuracy.json`。

只重复中等规模并独立保存测量证据：

```bash
python benchmarks/run_rz_performance.py --size medium \
  --results-directory verification/abaqus/rz_performance/medium_optimized \
  --windows-source '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance' \
  --windows-results '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance\medium_optimized'
```

该目录与优化前的原始样本分开保存。`provenance.json` 保存实际输入、运行脚本和
可执行文件的散列值；`timing.json` 保存一次预热和两次正式测量。目录外的旧
`provenance.json`、`timing.json`、`summary.json` 仍对应优化前配置。


正式计时结果如下，均值不包含预热，单位为秒。

| 程序 | 第一次正式测量 | 第二次正式测量 | 外部总时间均值 | 内部求解或分析时间均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim | 31.5170 | 32.0542 | **31.7856** | 31.2888 | 65 |
| Abaqus | 41.8601 | 41.7760 | **41.8180** | 37.0000 | 52 |

Fuelsim 外部运行时间为 Abaqus 的 0.7601 倍，耗时减少 **23.99%**，速度为
1.316 倍。两次 Fuelsim 正式运行均快于两次 Abaqus 正式运行。
比较边界仍是前述同一主机、WSL2 与 Windows、各自 CPU 0、单线程、20 个固定增量；
不是同一操作系统下的纯内核速度比较。

| 比较量 | 优化后最大逐点相对误差（%） |
|---|---:|
| 温度 | 2.69319307e-7 |
| 位移向量 | 3.65559911e-6 |
| 应力张量 | **4.71157019e-6** |
| 接触压力 | 4.06128122e-6 |
| 接触间隙 | 4.06366757e-6 |
| 逐节点法向力 | 4.06232464e-6 |
| 总法向力 | 3.46593512e-6 |

全部相对二范数、绝对峰值相对误差和最大逐点相对误差均小于 **0.01%**。
规定零位移仍单独满足 1e-12 m 的绝对门槛。覆盖范围保持为全部 7,670 个节点、
29,696 个材料积分点和 65 个接触节点的最终状态，位移和应力按向量与张量范数比较。

此次未修改 CMake、生产求解层或物理实现。相关输入解析、单进程和多进程求解器、
验证矩阵回归 **8/8 通过**（5.05 秒），命令为：

```bash
ctest --test-dir build -j4 --output-on-failure \
  -R '^fuelsim_(input|solver|hex8_solver|hex20_solver|mpi_solver|mpi_hex8_solver|mpi_hex20_solver|verification_matrix)_tests$'
```

检查程序另外验证 0.005% 误差通过、0.02% 误差失败，零参考值的绝对门槛独立有效，
空数据和非有限数被拒绝。新脚本的单规模运行和独立输出目录已通过上述完整计时运行验证。
