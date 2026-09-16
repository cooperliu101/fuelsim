# C3D20RT 中等规模例题的误差定位与修正

本次以提交 `149dea569b59f506f6890db6fbf3c28048a9a756` 为基线，继续比较
1,152 个体单元、19,524 个自由度模型的全部二十个时间步。基线完整指标见
[before_comparison.csv](before_comparison.csv)。该轮修正使用原来的严格 0.1%
要求；零参考量保留原绝对容差，非零参考量不设置相对误差分母下限。
用户随后要求极小量与绝对容差比较；当前验收规则和完整结果见
上级目录的 [完整报告](../README.md)。下文保留物理修正时的原始比较记录。

## 先用同一节点历程隔离误差

手动诊断程序 `fuelsim_c3d20rt_state_replay` 读取 Abaqus 的全部节点温度和
位移历程，将同一历程逐步交给 Fuelsim 的材料更新和残量计算。它不求解平衡，
不用于替代 `fuelsim -i` 的生产程序验收。材料点只根据原生节点坐标重建的
积分点位置配对，并逐单元检查一一对应。

修正前，在这个指定节点历程上，燃料应力的最大逐点相对误差仅为
`2.86114e-7%`，包壳材料更新也接近原生结果；但 0.15 s 的界面总传热量
仍相差 `4.16602%`。前十一帧的燃料内部热残量最大值低于 `9e-11 W`。
这说明应优先检查热接触及随后升温阶段的热物性表示，不能把耦合计算中的
燃料应力误差直接归因于体单元应力更新。完整诊断见
[before_state_replay.csv](before_state_replay.csv)。

## 修正热接触的平均间隙

二十节点单元的温度只在角点求解，机械接触同时使用角点和边中点约束。
旧实现用整个温度积分区域的采样点平均间隙。原生节点历程中的传热量表明，
本例题应采用温度角点所对应的机械平均间隙：只有该角点的机械积分规则
参与间隙平均，而边中点规则仍参与温度平均和传热面积。

有限滑移的间隙还要投影到次接触面随变形更新的单位法向。主接触面用于
定位投影点和计算温度插值，不能用它的法向替代曲面上的次面法向。

代码因此分别累计间隙平均面积和温度平均面积，并保留两个商的完整几何
导数。全局装配预先保存每个采样点属于哪个机械约束，以及次面参考法向
的朝向。热接触可以单独启用，不依赖机械接触残量或历史是否启用。

采用相同 Abaqus 节点历程后，界面总传热量的最大相对差降至
`0.000484723%`。局部倾斜平面测试同时检查不同法向、不同平均面积、
解析间隙、两侧热量守恒及残量导数；中心差分方向导数相对误差为
`7.92247e-11`。全局曲面检查覆盖角点间隙与热量的关系、热接触独立启用、
装配后的几何导数，以及两个进程的结果一致性。装配后的中心差分方向
导数相对误差为 `5.33975e-9`，原始记录见
[derivative_checks.txt](derivative_checks.txt)。

这一修正作用于共用的二十节点单元有限滑移热接触路径。它没有修改体单元
积分、本构模型、固定初始质量规则、机械接触的间隙定义或力学参数。
小滑移仍使用原有的热积分规则。当前原生证据来自本圆柱例题的完整历程，
不能据此声称所有曲面和材料组合均已验证。

## 加密 Abaqus 的导热系数表

燃料导热系数的物理公式为 `k(T) = 3824/T + 0.61`。原 Abaqus 表在
590 至 620 K 每隔 0.1 K 取值，但 620 K 后直接跳到 1200 K。原生模型
实际升温至约 649.65 K，已超出原先加密的温度区间。该区间的线性插值
与物理公式最大相差约 `1.98746%`。

两张完整 Abaqus 输入卡均新增 620.1 至 660.0 K 的四百行数据，每行直接
由原公式求值，没有根据验证结果拟合参数。在原实际温度范围内，表格的
最大相对插值误差降至约 `6.34e-7%`。逐帧温度范围与插值检查见
[conductivity_table_check.csv](conductivity_table_check.csv)。该表的温度
上限来自修改前原生历程；更新参考的最高温度为 649.664845 K，同样在
加密范围内。

随后重新运行 Windows 原生 Abaqus 2025，完成全部二十个增量，并重新
提取全部节点、接触与材料点字段。新的运行目录、输入哈希与原生记录见
上级目录的 `provenance.txt`、`SHA256SUMS` 和
[abaqus_runner.txt](abaqus_runner.txt)。更新后的指定历程诊断见
[after_state_replay.csv](after_state_replay.csv)；全部二十帧燃料内部热残量
最大值为 `2.95210e-10 W`。

为分清代码与参考表的影响，保留了完整交叉比较：

- [code_only_comparison.csv](code_only_comparison.csv) 使用修正代码的生产结果，
  与旧 Abaqus 参考比较。
- [reference_only_comparison.csv](reference_only_comparison.csv) 使用旧 Fuelsim
  生产结果，与更新后的 Abaqus 参考比较。
- [before_solver_tightening_comparison.csv](before_solver_tightening_comparison.csv)
  同时采用修正代码和更新参考，仍使用原非线性收敛条件。

前两种组合仍存在非零有限量超过 0.1% 的比较结果，不能只修正其中一项。

## 求解收敛条件与最终生产比较

原收敛条件下，0.55 s 的 5843 号自由节点仍有约 `4.87159e-4 N` 的
未平衡力，超过零参考反力的 `1e-5 N` 绝对容差。Fuelsim 的
`reaction_force_*` 保存所有节点的未缩放残量：在约束节点上它是反力，
在自由节点上它是未平衡力。因此应收紧求解收敛条件，不能把自由节点输出
强制改成零来通过比较。

本例题的两张 Fuelsim 输入卡统一将全局绝对收敛容差从 `1e-7` 改为
`1e-9`，相对收敛容差从 `1e-8` 改为 `1e-10`，机械残量绝对收敛容差
从 `1e-4` 改为 `1e-6`。温度残量绝对收敛容差仍为 `1e-8`。
本输入未启用分字段收敛判断，实际停止使用全局残量和步长条件；机械
未平衡力是否达标，仍由生产结果的逐节点原始残量检查。
这些是求解停止条件，该轮物理修正没有改变比较程序的验收门槛。

最终结果和剩余边界见上级目录的 [完整报告](../README.md)、
[comparison.csv](../comparison.csv) 和
[comparison_by_frame.csv](../comparison_by_frame.csv)。

## 复现指定节点历程的诊断

在仓库根目录执行以下命令。它是独立诊断目标，不注册为生产例题测试。

```bash
cmake --build build --target fuelsim_c3d20rt_state_replay --parallel 4
gzip -dc verification/abaqus/2025/c3d20rt_medium/reference_nodes.csv.gz > /tmp/c3d20rt_nodes.csv
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  taskset -c 0 build/fuelsim_c3d20rt_state_replay \
  verification/fuelsim/transient_c3d20rt_medium_friction.fsi \
  /tmp/c3d20rt_nodes.csv /tmp/c3d20rt_replay
python3 verification/abaqus/2025/c3d20rt_medium/diagnose_replay.py \
  /tmp/c3d20rt_replay verification/moose/m58_integrated_hex20_mesh.e \
  /tmp/c3d20rt_replay_summary.csv
```

修改前的完整参考、输入和源码可从基线提交取得。复算基线时应使用该提交
的源码和参考；不要把更新后的参考误当作旧参考。
