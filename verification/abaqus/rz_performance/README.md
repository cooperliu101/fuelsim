# 轴对称生产输入性能与 Abaqus 验证

原 `fuelsim_m1_single_core_benchmark` 已删除，现使用 `fuelsim -i` 运行完整输入卡。
基础版本检查小应变稳态热弹性接触，不包含塑性或蠕变；有限应变 CAX4T 版本及其未通过的精度结果见文末。

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


## 中等规模完整牛顿更新优化（2026-09-09，载荷预测前）

中等规模两张完整输入卡使用 `line_search = basic`，即优先使用完整牛顿更新；
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


## 中等规模载荷预测优化（2026-09-09）

当前中等规模验证和计时输入卡在 `[Executioner]` 中显式启用
`use_linear_load_predictor = true`。执行器使用最近两个已接受平衡解及其载荷因子，
线性外推下一次节点场初值。前两次成功求解保持原来的初值策略，之后按实际载荷间距
外推，指定的温度和位移边界值重新施加。材料和接触历史不做预测，仍从完整的已提交
状态求解。预测失败时恢复该状态，从上一个平衡解重试同一载荷；两次均失败才缩小
载荷增量，失败尝试的迭代和耗时也进入总计。

该功能默认关闭。本次仅中等规模两张输入卡启用，不修改大规模输入、载荷路径、
MUMPS 设置、20 个固定增量、收敛容差或 0.01% 精度门槛。

| 加载增量 | 完整牛顿更新 | 加上载荷预测 |
|---|---:|---:|
| 1—2，各增量 | 3 | 3 |
| 3—10，各增量 | 3 | 2 |
| 11 | 5 | 4 |
| 12 | 3 | 2 |
| 13 | 5 | 3 |
| 14 | 4 | 3 |
| 15—20，各增量 | 3 | 2 |
| 合计 | **65** | **46** |

非线性迭代和 Jacobian 计算次数减少 29.23%，残量计算次数从 85 降到 66。
18 次预测全部成功，载荷增量没有缩小。Abaqus 使用相同的 20 个加载增量，记录为
52 次迭代。逐次残量和最终场验证保存在 `medium_load_predictor/`。
最终温度、自由位移向量、应力张量、接触压力、间隙、逐节点法向力和总法向力的
三项相对误差全部通过 0.01% 门槛；最大逐点相对误差为 4.71168438e-6%，出现在应力张量。
该结论仍只覆盖全部节点、材料积分点和接触节点的最终状态，不扩大为逐应力分量或
中间加载状态的相对误差保证。

复现正式计时的命令与上一节相同，将两处结果目录的末级名称改为
`medium_load_predictor`；运行器仍只选择 `--size medium`，先预热再测量两次。
`medium_optimized/` 的输入散列和测量值属于前一阶段，不以当前输入重新解释旧记录。


验证方面，Release 构建及统一回归 271/271 通过（167.29 秒）。新增内部检查确认
默认关闭预测，线性问题的预测保持各节点场不变且减少迭代；非等间距载荷检查随后
强化为每次最多一次牛顿迭代，实际产生 220 次预测、180 次预测失败恢复，最终温度、
径向位移和轴向位移分别与基准平衡解一致。强化后的测试单独复查通过，耗时 4.42 秒。
这些检查继续使用原有内部测试程序和共享输入卡，没有增加用户测试例题。


正式外部计时如下。每个程序先预热一次，再测量两次，顺序运行；计时期间没有编译
或其他求解任务。单核单线程、MUMPS、结果文件输出关闭以及跨操作系统边界与前述一致。

| 程序 | 第一次正式测量 | 第二次正式测量 | 外部总时间均值 | 内部求解或分析时间均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim，载荷预测 | 22.5491 s | 22.4468 s | **22.4980 s** | 21.9932 s | **46** |
| Abaqus，本次重测 | 39.8064 s | 39.6734 s | **39.7399 s** | 36.0000 s | **52** |

Fuelsim 耗时为本次 Abaqus 的 0.5661 倍，减少 **43.39%**，速度为 1.766 倍。
相对前一阶段记录的 Fuelsim 31.7856 秒，耗时进一步减少 **29.22%**。
前一阶段的 Abaqus 均值为 41.8180 秒，本节使用本次重新测得的 39.7399 秒比较，
不将两个批次的 Abaqus 时间混用。两次 Fuelsim 正式运行均完成 20 个加载增量、
46 次迭代和 18 次成功预测，预测失败次数及载荷缩小次数均为零。
原始日志、三项精度指标、逐增量迭代、可执行文件与源文件散列值均保存在
`medium_load_predictor/`，其 `summary.json` 汇总本次结果。


## ADlite 0.2.3 更新后重测（2026-09-09）

当前依赖升级到 ADlite 0.2.3，固定提交为
`fd319e00234e18280319d141f17d9fa015c2501b`，替代 0.2.2 的
`6b8af513aa5abb956a1b246d4d1a5c4f5fc7d8a2`。生产数值实现、求解器和中等规模两张
输入卡均未改变。两版本使用相同 Conda 编译器、Release 链接时优化、CPU 0 和单线程。
升级前保留的可执行文件重新预热一次并测量两次，作为本次直接比较基线；不只引用
上一节的历史运行时间。安装新版本时保留旧依赖目录，配置时显式更新缓存中的
`adlite_DIR`，并检查最终链接命令使用 `adlite-0.2.3/lib/libadlite.a`。

全部证据保存在 `medium_adlite023/`。`adlite022_timing.json` 和三个对应日志记录
升级前重测，`timing.json` 记录升级后的程序及本次 Abaqus 重测；`upgrade.json`
记录两个版本、库文件和升级前可执行文件的散列，`provenance.json` 记录升级后程序。
复现命令沿用 `--size medium`，将两处结果目录的末级名称指定为 `medium_adlite023`。


升级后的中等规模仍完成 20 个加载增量、46 次迭代和 18 次成功预测。完整精度报告
与 0.2.2 版本的 `medium_load_predictor/accuracy.json` 逐字节一致，最大逐点相对误差
为 4.71168438e-6%，全部规定指标继续满足 0.01%。这说明本算例的误差指标和求解
工作量保持一致，正式时间差不会混入加载步或非线性迭代次数变化。


新版本安装脚本自检通过，Fuelsim Release 完整回归 **271/271 通过**（166.37 秒）。
依赖版本及固定提交已同步到 CMake、安装脚本、持续集成配置、可复现构建检查和文档。
旧版程序来自 `7320aa2`，其散列与前一阶段已归档的程序一致；对应源码和旧版 ADlite
均可用于重建基线。升级前后输入卡保持逐字节不变。


本次正式运行结果如下，均不包含预热：

| 程序与依赖 | 第一次正式测量 | 第二次正式测量 | 外部总时间均值 | 内部求解或分析时间均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim，ADlite 0.2.2，本次重测 | 22.7744 s | 22.8797 s | **22.8271 s** | 22.3187 s | 46 |
| Fuelsim，ADlite 0.2.3 | 22.3808 s | 22.5081 s | **22.4445 s** | 21.9005 s | 46 |
| Abaqus，本次重测 | 41.8360 s | 41.8069 s | **41.8215 s** | 38.0000 s | 52 |

升级后的 Fuelsim 比本次旧库基线少用 **0.3826 秒（1.68%）**，改善幅度较小。
与上一节历史均值 22.4980 秒也很接近，因此不据此声称普遍或大幅加速。
相对本次 Abaqus，Fuelsim 耗时减少 **46.33%**。该比例使用当前批次的 Abaqus
41.8215 秒；前一批的 Abaqus 为 39.7399 秒，两批比例的差异不能全部归因于 ADlite。
两版本 Fuelsim 的全部正式运行均为 20 个增量、46 次迭代、18 次成功预测，
精度门槛和最终场覆盖范围与前述一致，完整记录见 `medium_adlite023/summary.json`。


## ADlite 0.2.3 关闭显式 SIMD 后端对照（2026-09-09）

本次固定 ADlite 提交 `fd319e00234e18280319d141f17d9fa015c2501b`，使用
`ADLITE_ENABLE_SIMD=OFF` 构建独立的 `adlite-0.2.3-nosimd` 安装目录。
原源码目录出现的未提交修改未纳入此对照；标量版本从固定提交的干净副本构建。
原有启用版本的库及可执行文件保留，输入卡、求解器、编译器、Release 链接时优化、
单核和单线程设置均保持一致。这里关闭的是 ADlite 手写 SIMD 后端，不额外禁止
编译器自动向量化，也不改变 PETSc 或 MUMPS 的构建方式。

构建记录确认 `ADLITE_HAS_AVX2_BACKEND=0`。链接实际安装库的探测还确认宽度 5 和 6
的缩放、累加和组合运算都选择标量计算函数，见 `medium_simd_off/dispatch.txt`。
本次对照测试期间，本地 `build/fuelsim` 使用关闭版本。对照完成后按用户要求恢复
默认启用 SIMD 的 `adlite-0.2.3` 安装库；仓库默认安装配置始终维持启用 SIMD。
恢复后重新完成 Release 构建及统一回归，271/271 项测试通过；运行时探测确认
宽度 5 和 6 的缩放、累加和组合运算选择 AVX2 后端。

复现时先按常规方式构建启用版本，并保存为 `/tmp/fuelsim-adlite023-simd-on`。
关闭配置可通过先创建安装脚本的构建缓存来复现：

```bash
git clone /home/cooper/ai_project/ADlite /tmp/fuelsim-adlite023-fixed-source
git -C /tmp/fuelsim-adlite023-fixed-source checkout fd319e00234e18280319d141f17d9fa015c2501b
cmake -S /tmp/fuelsim-adlite023-fixed-source \
  -B /tmp/fuelsim-adlite023-fixed-nosimd-build \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_BUILD_TYPE=Release -DADLITE_ENABLE_SIMD=OFF
scripts/install_adlite.sh /tmp/fuelsim-adlite023-fixed-source \
  /home/cooper/ai_project/fuelsim-dependencies/adlite-0.2.3-nosimd \
  /tmp/fuelsim-adlite023-fixed-nosimd-build /home/cooper/miniforge/envs/moose
```

上述源码目录须为该固定提交的干净 Git 工作树。Fuelsim 使用常规 Release 配置，
同时将 `CMAKE_PREFIX_PATH` 和 `adlite_DIR` 指向该版本化安装目录及其 `lib/cmake/adlite`。
全部构建仍使用 `--parallel 4`，统一回归使用 `ctest --test-dir build -j4 --output-on-failure`。

正式计时使用同一张完整输入卡，不生成或修改问题定义：

```bash
python benchmarks/run_rz_simd_comparison.py \
  --simd-on /tmp/fuelsim-adlite023-simd-on --simd-off build/fuelsim
```

运行器先分别预热，再按“启用、关闭”和“关闭、启用”的顺序测量两轮，外部单调时钟
覆盖程序启动和完整求解。`medium_simd_off/` 保存两版本程序散列、逐次记录及精度证据。
本次不重测 Abaqus，避免将先前批次的 Abaqus 时间当成本轮 SIMD 对照的直接基线。


关闭版本的 Release 构建及 ADlite 自检通过，Fuelsim 统一回归 **271/271 通过**
（167.94 秒）。精度报告与启用版本逐字节一致，最大逐点相对误差仍为
4.71168438e-6%，满足 0.01% 门槛。

| ADlite 显式 SIMD 后端 | 第一轮正式测量 | 第二轮正式测量 | 外部总时间均值 | 内部求解时间均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| 启用 | 22.5567 s | 22.4650 s | **22.5109 s** | 22.0014 s | 46 |
| 关闭 | 22.2189 s | 22.1642 s | **22.1916 s** | 21.6555 s | 46 |

本次关闭版本平均少用 **0.3193 秒（1.42%）**，两轮均略快，但差异较小。
该观察只对应本中等规模算例，不据此推广其他单元、导数宽度或模型的性能。
全部六次运行都完成 20 个增量、46 次迭代，且没有预测失败。
原始样本及汇总见 `medium_simd_off/timing.json` 和 `medium_simd_off/summary.json`。


## elements 重构后中等规模复测（2026-09-10）

测量程序为 `da5ad75`，ADlite 0.2.3 开启 SIMD，实际链接库在当前 CPU 上的
导数宽度 5、6、7、10 均选择 AVX2。沿用完整生产输入、MUMPS、单进程单线程、
CPU 0 绑定及关闭结果输出的计时条件。两套程序各预热一次，再正式测量两次。

| 程序 | 第一次外部耗时 | 第二次外部耗时 | 外部耗时均值 | 内部求解或分析均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim | 22.6176 s | 22.2837 s | **22.4507 s** | 21.9486 s | 46 |
| Abaqus | 42.4289 s | 48.7008 s | **45.5648 s** | 37.5000 s | 52 |

按本批正式测量均值，Fuelsim 耗时减少 **50.73%**，
Abaqus/Fuelsim 耗时比为 **2.03**。Abaqus 两次分析时间分别为 37、38 秒，
外部总时间的波动更大；不能把外部总时间与分析时间之间的差额完全归因于许可证。
这是同一主机的 WSL/Windows 跨平台观察，CPU 编号不保证对应同一物理核心，
也不代表多次统计基准或所有算例的加速比例。

当前生产程序另外运行逐字复制的完整验证输入卡，最终场与受追踪的 Abaqus 参考
比较，所有指标均满足 0.01% 门槛；最大逐点相对误差为
**4.71168438e-06%**。本次重新运行的是 Abaqus 计时作业，
精度检查使用已有最终场参考，未重新生成 Abaqus 场输出。

原始计时、压缩运行日志、精度报告及输入和程序散列位于
[`medium_elements_refactor/`](medium_elements_refactor/summary.json)。未修改源码或物理输入。


## 中等规模 CAX4RT 版本（2026-09-10）

新增完整输入 `verification/fuelsim/steady_rz_performance_medium_cax4rt.fsi` 及
对应 `_timing.fsi`。与 CAX4T 相比仅改变单元型号、注释和结果文件名，网格仍为
7,424 个单元、7,670 个节点、23,010 个自由度；材料、边界、接触、20 个固定加载
增量、MUMPS 和载荷预测设置不变。Abaqus 新增对应 CAX4RT 输入与同节点连接的网格。
CAX4RT 使用一个活跃材料积分点及其固有沙漏控制，不能用 CAX4T 四点应力作参考。

运行完整输入：

```bash
build/fuelsim -i verification/fuelsim/steady_rz_performance_medium_cax4rt.fsi
python benchmarks/run_rz_performance.py --size medium --element cax4rt \
  --results-directory verification/abaqus/rz_performance/medium_cax4rt \
  --windows-source '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance' \
  --windows-results '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance\medium_cax4rt'
```

重新测量时应指定新的结果目录以保留本次证据。`run.ps1 -Size medium -Element cax4rt`
在不带 `-Timing` 时运行并提取 CAX4RT 最终场。`compare.py` 使用 `--element cax4rt`
选择一个材料积分点，默认 CAX4T 仍检查四点；对 CAX4RT 输出错误选择 CAX4T 会被拒绝。
本次另行核对了不活跃的 12 个应力字段全部为 NaN，原 CAX4T 精度检查仍通过。

单进程单线程、绑定逻辑 CPU 0、ADlite 0.2.3 SIMD 开启、关闭结果输出，预热一次后
正式测量两次，外部计时包括启动、输入读取和完整求解：

| 程序 | 第一次外部耗时 | 第二次外部耗时 | 外部耗时均值 | 内部求解或分析均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim CAX4RT | 10.9232 s | 11.5372 s | **11.2302 s** | 10.4939 s | 46 |
| Abaqus CAX4RT | 31.7888 s | 29.7110 s | **30.7499 s** | 26.5000 s | 53 |

Fuelsim 比本批 Abaqus 耗时少 **63.48%**，耗时比为
**2.74**；相较紧邻前一批 Fuelsim CAX4T 的 22.4507 秒，减少
**49.98%**。该跨型号比较属于不同积分算法的
同模型运行观察，并非同一离散算子的代码优化。两套程序分别运行于 WSL 和 Windows，
逻辑 CPU 编号不保证同一物理核心，不将此比例推广到其他模型。

本次重新运行 Abaqus CAX4RT 输出作业生成对应参考，同时单独运行 Fuelsim 验证输入。
温度、位移、应力张量、接触压力、间隙和接触反力的全部比较指标通过原有 0.01% 门槛，
最大逐点相对误差为 **4.53764245e-06%**。仅比较最终加载状态。
结果、原始压缩日志和散列见 [`medium_cax4rt/summary.json`](medium_cax4rt/summary.json)。
未修改单元数值实现；此算例为手动性能验证，没有新增重复的自动回归。


## 中等规模 CAX8T 版本（2026-09-10）

新增完整生产输入 `steady_rz_performance_medium_cax8t.fsi` 和对应 `_timing.fsi`。
沿用中等规模的 7,424 个单元分区，在每条边增加共享的位移中间节点，保持原有
角点坐标、角点连接、材料、边界、接触参数和 20 个固定加载增量。单元数相同，
但自由度数不同，不能将跨型号时间差称为同规模的算法性能差异。

| 量 | CAX4T/CAX4RT | CAX8T |
|---|---:|---:|
| 单元数 | 7,424 | 7,424 |
| 几何与位移节点数 | 7,670 | 22,762 |
| 温度自由度数 | 7,670 | 7,670 |
| 总自由度数 | 23,010 | 53,194 |
| 每单元材料积分点数 | 4 / 1 | 9 |
| 接触从属节点数 | 65 | 129 |

网格生成器仅写几何，使用
`python verification/abaqus/generate_rz_performance_meshes.py --element cax8t`。
Fuelsim 和 Abaqus 使用同一组角点、中间节点和单元连接；本次检查了所有中间节点
恰好位于对应边的中点。物理输入卡已完整保存，运行器不会生成或修改输入定义。

```bash
build/fuelsim -i verification/fuelsim/steady_rz_performance_medium_cax8t.fsi
python benchmarks/run_rz_performance.py --size medium --element cax8t \
  --results-directory verification/abaqus/rz_performance/medium_cax8t \
  --windows-source '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance' \
  --windows-results '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance\medium_cax8t'
```

重新测量时使用新的结果目录以保留原有证据。`run.ps1 -Size medium -Element cax8t`
不带 `-Timing` 时运行并提取 Abaqus CAX8T 的最终场。

沿用单进程、单线程和逻辑 CPU 0 绑定，Fuelsim 使用 MUMPS 和开启 SIMD 的
ADlite 0.2.3；两套程序关闭结果输出，各预热一次，再正式测量两次：

| 程序 | 第一次外部耗时 | 第二次外部耗时 | 外部耗时均值 | 内部求解或分析均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim CAX8T | 55.4502 s | 55.6876 s | **55.5689 s** | 54.9603 s | 48 |
| Abaqus CAX8T | 104.1527 s | 106.0978 s | **105.1253 s** | 100.5000 s | 52 |

Fuelsim 外部耗时少 **47.14%**，Abaqus/Fuelsim 耗时比
为 **1.89**。这是同一主机的 WSL/Windows 跨平台观察，外部时间包含
程序启动等开销，逻辑 CPU 0 不保证同一物理核心，不推广到其他模型或硬件。

精度使用 `compare.py --element cax8t`，核对所有 7,670 个温度角点、22,762 个
位移节点、66,816 个材料积分点以及 129 个接触节点。温度自由度掩码由网格角点
独立核对，边中节点温度另检查为角点线性插值。Abaqus 没有独立热自由度的节点
允许缺少热反力或温度字段，不把 NaN 当成数值参考；所有角点温度仍严格检查。
九点应力使用与现有 CAX8T 验证相同的积分点顺序。节点压力使用
`contact_recovered_pressure_fuel_cladding` 对比 Abaqus CPRESS，遵循
[既有 CAX8T 恢复规则](../b114_cax8t_recovery_validation.md)，不改变单元或接触计算。

全部最终场指标通过原有 0.01% 门槛，最大逐点相对误差为
**4.58855909e-06%**。Fuelsim 的 129 个接触节点全部活跃，
法向反力合计约 664.04 N。本次重新运行了 Abaqus CAX8T 输出作业，未借用其他
型号的应力参考；不声称逐增量场均已对比。原 CAX4T/CAX4RT 比较脚本检查仍通过。

完整证据见 [`medium_cax8t/summary.json`](medium_cax8t/summary.json)。未修改生产
数值实现；新增算例作为手动性能验证，不加入重复的自动回归。


## 中等规模 CAX8RT 版本（2026-09-10）

完整生产输入为 `steady_rz_performance_medium_cax8rt.fsi` 和对应 `_timing.fsi`。
与 CAX8T 使用同一 `rz_performance_medium_cax8t.e` 二次网格，只有型号、注释和
输出文件名不同。7,424 个单元、22,762 个节点、53,194 个自由度、材料、接触、
小应变设置及 20 个固定加载增量保持一致。体积分从每单元九点改为四点。

```bash
build/fuelsim -i verification/fuelsim/steady_rz_performance_medium_cax8rt.fsi
python benchmarks/run_rz_performance.py --size medium --element cax8rt \
  --results-directory verification/abaqus/rz_performance/medium_cax8rt \
  --windows-source '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance' \
  --windows-results '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance\medium_cax8rt'
```

重新测量时使用新的结果目录。Abaqus 输出作业使用 `run.ps1 -Size medium -Element cax8rt`
且不带 `-Timing`，精度检查使用 `compare.py --element cax8rt`。

沿用单进程、单线程、逻辑 CPU 0、MUMPS 和 ADlite 0.2.3 SIMD 设置，关闭结果输出，
各预热一次，随后正式测量两次：

| 程序 | 第一次外部耗时 | 第二次外部耗时 | 外部耗时均值 | 内部求解或分析均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim CAX8RT | 44.2350 s | 43.4817 s | **43.8583 s** | 43.2325 s | 48 |
| Abaqus CAX8RT | 86.1217 s | 86.0712 s | **86.0964 s** | 81.5000 s | 53 |

Fuelsim 外部耗时少 **49.06%**，耗时比为 **1.96**。
相对上一批同网格 Fuelsim CAX8T 的 55.5689 秒，减少
**21.07%**；这属于不同积分规则的同模型比较，
不声称是相同数值算子的代码优化。仍是 WSL/Windows 跨平台观察，包含启动开销，
逻辑 CPU 编号不保证同一物理核心，不推广到其他模型或硬件。

本次重新生成 Abaqus CAX8RT 最终场，检查所有 7,670 个角点温度、22,762 个节点
位移、29,696 个活跃积分点应力和 129 个接触节点。四点顺序按现有轴对称规则映射为
`[0,1,3,2]`；二次接触压力沿用已验证的恢复压力。全部指标通过原有 0.01% 门槛，
最大逐点相对误差 **4.63968479e-06%**。

129 个接触节点全部活跃。预留 q4—q8 的 30 个应力及坐标字段全部为 NaN，
不参与统计。原 CAX4T、CAX4RT、CAX8T 精度检查复查仍通过。只比较最终状态，
未修改生产数值代码，未增加重复自动回归。原始证据见
[`medium_cax8rt/summary.json`](medium_cax8rt/summary.json)。


## 中等规模 CAX4T 有限应变对比（2026-09-10）

**以下为原 steady 执行方式的未通过记录。** 后续已定位增量历史差异，
匹配 Abaqus 的增量执行结果见下节。原记录保留，不覆盖其误差或性能数据。
两套程序均完成 20 个固定加载增量，
但最终场不满足既有 0.01% 门槛。没有修改单元数值实现、物理参数或误差门槛。

新增 `steady_rz_performance_medium_cax4t_finite.fsi` 和对应 `_timing.fsi`。
沿用原 CAX4T 中等规模网格：7,424 个单元、7,670 个节点、23,010 个自由度。
两个区域均改为 `strain = finite`，Abaqus 设置 `nlgeom=YES`；材料、接触、
载荷、20 个增量、MUMPS、载荷预测和求解器容差保持不变。

```bash
build/fuelsim -i verification/fuelsim/steady_rz_performance_medium_cax4t_finite.fsi
python benchmarks/run_rz_performance.py --size medium --element cax4t --strain finite \
  --results-directory verification/abaqus/rz_performance/medium_cax4t_finite \
  --windows-source '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance' \
  --windows-results '\\wsl.localhost\Ubuntu\home\cooper\ai_project\fuelsim\verification\abaqus\rz_performance\medium_cax4t_finite'
```

重复测量时选择新的结果目录。Abaqus 参考使用 `run.ps1 -Size medium -Element cax4t
-Strain finite` 且不带 `-Timing` 重新生成。比较脚本不变，使用 `--element cax4t`。
单进程、单线程、逻辑 CPU 0、ADlite 0.2.3 SIMD 开启、关闭结果输出，各预热一次后
正式测量两次，得到以下实际耗时：

| 程序 | 第一次外部耗时 | 第二次外部耗时 | 外部耗时均值 | 内部求解或分析均值 | 非线性迭代数 |
|---|---:|---:|---:|---:|---:|
| Fuelsim | 28.4315 s | 28.5623 s | **28.4969 s** | 28.0025 s | 46 |
| Abaqus | 43.7750 s | 45.8381 s | **44.8066 s** | 39.5000 s | 53 |

本批外部耗时比为 1.57，Fuelsim 实际耗时少
36.40%。由于精度未通过，这不是精度达标后的性能结论。
WSL/Windows 跨平台、启动开销以及逻辑 CPU 编号的限制仍适用。

最终场检查覆盖全部节点、29,696 个材料积分点和 65 个接触节点；最大逐点误差如下：

| 比较量 | 最大逐点相对误差 | 原 0.01% 门槛 |
|---|---:|---|
| 温度 | 0.011720% | 未通过 |
| 自由位移向量 | 0.079766% | 未通过 |
| 应力张量 | 0.773119% | 未通过 |
| 接触压力 | 0.159140% | 未通过 |
| 接触间隙 | 0.159140% | 未通过 |
| 节点法向反力 | 0.159157% | 未通过 |
| 接触总反力 | 0.049146% | 未通过 |

规定零位移的绝对误差检查通过。65 个接触节点全部活跃，Fuelsim 法向反力合计约
667.46 N。最大应力误差位于芯块第 6001 单元的 Abaqus 第 2 积分点，参考张量
范数约 9.44 MPa、差值约 0.0730 MPa；不是零参考值舍入造成的误差。该 0.7731%
最大应力误差也超过轴对称既有 0.5% 逐点门槛。该批测量时尚未定位差异原因。

原始日志、独立 Abaqus 参考、完整误差指标和最大应力点信息保存在
[`medium_cax4t_finite/`](medium_cax4t_finite/summary.json)。本次增加的是手动对比输入与
未通过结果记录，没有把未通过算例加入自动回归，也没有放宽任何验收条件。

## CAX4T 有限应变增量历史诊断与匹配结果

**全部最终场指标通过，最大逐点相对误差为 0.00000430214%。**
这低于本次要求的 0.1%，也通过比较脚本原有的 0.01% 门槛。
全部 7,670 个节点、29,696 个材料积分点和 65 个接触节点均参与比较。
没有剔除原最大应力误差点，也没有改变材料、网格、载荷、接触参数或误差算法。

原因是两种执行方式的有限应变历史不同。`SteadyProblem::compute_contribution`
调用 `compute_cax4_thermoelastic`，后者向单元传入零的已接受节点状态和空材料历史；
`SteadyProblem::commit_internal_state` 只提交接触历史。因此原输入中的 20 个
载荷步用于逐步求解，但体单元每次仍从初始构形计算有限应变。
Abaqus 的 20 个增量则逐次累积材料历史，并使用相邻增量的中间构形计算热传导。

两个独立计算验证了这个判断：

1. 新增完整 Abaqus 输入
   `rz_performance_medium_cax4t_finite_single_increment.inp`，只把原 20 个增量
   改为一次达到同一最终载荷。原 Fuelsim steady 结果与它的所有最终场指标
   通过，最大逐点相对误差为 0.00000445988%。
2. 新增完整 Fuelsim 输入
   `verification/fuelsim/quasistatic_rz_performance_medium_cax4t_finite.fsi`，
   使用已有 transient 执行器提交每个接受增量的节点状态和材料历史。
   `include_thermal_time_term = false` 关闭热容项，因此每个增量仍求解稳态
   热平衡；本算例也不包含位移惯性、塑性或蠕变。固定增量为 1，结束时间和
   载荷上升时间均为 20，与原 Abaqus 参考的 20 个等分增量一致。

这是算例执行配置的修正，没有修改生产单元或求解器实现，也没有改变
`SteadyProblem` 的现有行为。需要累积有限应变历史的本算例使用新输入。
该路径完成 20 个增量，无失败重试，共 46 次非线性迭代，65 个接触节点全部活跃。

| 比较量 | 原 steady 最大逐点误差 | 匹配增量历史后的最大逐点误差 |
|---|---:|---:|
| 温度 | 0.011720% | 0.000000265783% |
| 自由位移向量 | 0.079766% | 0.00000336718% |
| 应力张量 | 0.773119% | 0.00000430214% |
| 接触压力 | 0.159140% | 0.00000393020% |
| 接触间隙 | 0.159140% | 0.00000392508% |
| 节点法向反力 | 0.159157% | 0.00000393122% |
| 接触总反力 | 0.049146% | 0.00000310324% |

规定零位移继续按绝对误差检查并通过。相对 L2、相对绝对峰值和最大逐点
相对误差全部保存在
[`medium_cax4t_finite_incremental/accuracy.json`](medium_cax4t_finite_incremental/accuracy.json)。
参考仍为原 `medium_cax4t_finite/` 中的 20 增量 Abaqus 数据，没有重新拟合或替换。
检查脚本的 `--incremental` 选项要求初始状态及 20 个单位增量均存在，
最终载荷因子为 1；跨程序精度结论只覆盖最终状态，不声称已逐增量比较全部历史场。

```bash
env -u PETSC_OPTIONS OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 taskset -c 0 build/fuelsim \
  -i verification/fuelsim/quasistatic_rz_performance_medium_cax4t_finite.fsi
env NETCDF_LIBRARY=/home/cooper/miniforge/envs/moose/lib/libnetcdf.so \
  python verification/abaqus/rz_performance/compare.py \
  verification/fuelsim/quasistatic_rz_performance_medium_cax4t_finite_results.e \
  verification/abaqus/rz_performance/medium_cax4t_finite/rz_performance_medium_cax4t_finite \
  --element cax4t --incremental --report /tmp/cax4t-finite-accuracy.json
```

运行时可把完整输入逐字复制到隔离目录，并保持 `../meshes/` 网格路径有效。
单增量 Abaqus 诊断在包含原网格和材料 `.inc` 文件的目录中直接运行新 `.inp`，
再使用 `extract_results.py` 提取结果。归档诊断的作业名沿用
`rz_performance_medium_cax4t_finite`，其实际单增量输入由
[`provenance.json`](medium_cax4t_finite_incremental/provenance.json) 中的 SHA256 标识。

新输入的生产日志、可执行文件和输入摘要、单增量 Abaqus 原始日志及参考保存在
[`medium_cax4t_finite_incremental/`](medium_cax4t_finite_incremental/summary.json)。
原四种小应变中等规模算例的全部比较指标复查保持不变。
**新执行路径尚未做关闭输出后的重复计时，原 1.57 倍耗时比不能用于此处。**
