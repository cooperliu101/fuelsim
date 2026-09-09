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
二范数、相对绝对峰值和最大逐点相对误差，统一要求小于 0.5%，不设置分母下限。
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

## 本次结果（2026-09-09）

两次正式测量的均值如下，预热运行不参与均值。

| 规模 | Fuelsim 外部总时间 | Abaqus 外部总时间 | Fuelsim / Abaqus | Fuelsim 内部求解时间 | Abaqus 分析时间 | 非线性迭代数（Fuelsim / Abaqus） |
|---|---:|---:|---:|---:|---:|---:|
| medium | 70.34 s | 41.80 s | 1.683 | 69.84 s | 37.00 s | 145 / 52 |
| large | 163.90 s | 83.13 s | 1.971 | 163.37 s | 78.50 s | 158 / 52 |

当前配置下，两档 Fuelsim 都慢于 Abaqus。非线性迭代数也更多；本次记录这一事实，
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
