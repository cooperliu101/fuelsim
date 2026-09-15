# 多条轴向线的截面分区试验

本试验检查沿板宽增加相互耦合的轴向线，能否用更少的计算成本替代端部丰富模态。
它通过实际 `fuelsim -i` 入口运行，沿用原有四类载荷和完整三维参考，保留全域
位移、轴向应力、完整应力张量、弯矩、应变能的严格小于 1% 门槛。低于或高于
门槛都保存实际数值；局部误差分解仅用于诊断。

板厚 2 mm、板宽 20 mm、长度 200 mm，材料为 E=70 GPa、nu=0.3 的均匀线弹性体。
按用户选择，根部 z=0 改用三个不共线节点上的六个标量位移约束：

| 支撑位置，单位 m | 置零的位移分量 |
| --- | --- |
| A=(0,0,0) | u_x、u_y、u_z |
| B=(0.001,0,0) | u_y、u_z |
| C=(0,0.01,0) | u_z |

刚体位移 u=t+omega cross X 代入后，六个方程的矩阵满秩。其余根部节点保持自由，
不额外约束幅值导数。完整三维参考施加同样的六个约束；三维网格只增加三个
支撑节点集，坐标、连接、原有节点集及边集均逐项保持原样。
这是一种三点支撑；反力在离散节点传递，支撑点附近的局部应力仍计入全域比较。
四类载荷及完整三维网格见[生产比较说明](../verification/section_modal/README.md#同条件三维比较)。

## 离散与实现

几何、材料、轴向节点、截面积分和载荷与原端部模型相同。新网格 `plate_width.e`
保持 `plate_local.e` 的全部坐标与单元连接，增加五个宽度位置及三个支撑节点集。
位置全部来自 Exodus；生产程序不生成分区坐标。

选择三个节点集表示两个相连的宽度分区，选择五个表示四个分区。两侧边界及
内部交界分别对应 3 条、5 条轴向线。因此分区数和轴向线数相差一，不能混用。
另用只选择两条边界线的一个宽度分区作为对照。

每条线使用连续的分段线性权重 h_r(y)。相邻分区共享线上的位移表示，权重满足
sum_r h_r(y)=1。先调用已有截面求解器生成 12 个种子模式，补齐其 Phi1、Phi2 中
尚未独立表示的位移场，再采用已有的独立幅值运动学。这个选择保留经典拉伸、
两方向弯曲、泊松场和扭转翘曲；它与原来的导数耦合多项式离散不是简单的代数换元。

在已有 QUAD8 节点上定义分区场：

```text
Psi_ri = I_h[h_r(y) Phi_i(x,y)]
u(x,y,z) = sum_ri Psi_ri(x,y) a_ri(z)
```

I_h 表示使用原截面有限元的共享节点插值。乘积以该插值后的场为准；没有使用
超出 QUAD8 空间的解析乘积再与不同积分规则混合。分区边界必须沿截面单元边界。
原经典位移场显式保留，随后对候选乘积做两遍截面 L2 正交化，去除线性相关项，
并保持已验证的物理反射空间。最终输出的幅值是这个等价空间中的稳定坐标，
不应把每个幅值直接解释成某条线的某个物理位移。
宽度分区自身也必须满足相应反射；否则关闭该轴的反射分解，避免生成输入中
没有选择的镜像分区场。不对称分区和穿过截面单元内部的非法分区另有内部检查。

以一个场 Psi 为例，其应变包含：

```text
epsilon_xx = Psi_x,x a
epsilon_yy = Psi_y,y a
epsilon_zz = Psi_z a'
gamma_xy   = (Psi_x,y + Psi_y,x) a
gamma_xz   = Psi_z,x a + Psi_x a'
gamma_yz   = Psi_z,y a + Psi_y a'
```

代码仍使用六分量张量剪应变，以上 gamma 在进入材料接口时乘 1/2。
这些场的横向梯度产生线间耦合；单元刚度来自完整的三维材料应变能。
没有设置经验弹簧、独立梁之间的拟合连接或人工剪切刚度。
每个轴向单元仍使用五次 Hermite 插值和六点积分，材料点保持完整三维应变及应力。

全长均采用同一个分区空间，本轮关闭端部凝聚。一个、两个、四个分区分别形成
25、38、64 个独立截面方向，73 个轴向节点对应 5475、8322、14016 个幅值及导数
未知量；六个物理位移约束的乘子另计。输入中的 `count=12` 是种子模式数，输出
`mode_count` 是实际总方向数，`seed_mode_count` 和 `width_line_count` 分别记录
种子模式数与线数。`width_enrichment_modes` 还计入补齐的独立导数字段。

实现集中在 `enrich_section_width`，复用原材料、截面运动学、轴向装配、物理载荷
投影和场恢复。`width_lines` 输入暂时与端部凝聚互斥，作为研究功能使用。

## 验证与复现

内部检查覆盖经典拉伸、两方向弯曲、扭转的能量和逐点应变，分区加密后的空间
包含关系，截面正交性以及材料点残量的方向导数与虚功。生产比较还核对截面
单元公共面的位移连续、全部材料点覆盖、应力积分内力、物理约束和外功。

```bash
# 仅生成带宽度节点集的研究网格，生产输入卡均已完整提交。
python - <<'PY'
from pathlib import Path
from verification.section_modal.mesh import write
write(Path('verification/section_modal/plate_width.e'), nx=4, ny=8, nz=128,
      graded=True, local_interfaces=True, coarsen_interior=True, width_lines=True,
      point_support=True)
write(Path('verification/section_modal/plate_support_reference.e'), nx=4, ny=8, nz=1024,
      point_support=True)
PY

ctest --test-dir build -j8 --output-on-failure
python verification/section_modal/width_study.py \
  build/fuelsim build/blackbox/section_width_support build/blackbox/section_width_support/reference \
  --reference-only --cpu 0
python verification/section_modal/width_study.py \
  build/fuelsim build/blackbox/section_width_support build/blackbox/section_width_support/reference --cpu 0
python verification/section_modal/width_study.py \
  build/fuelsim build/blackbox/section_width_support build/blackbox/section_width_support/reference \
  --mpiexec /home/cooper/miniforge/envs/moose/bin/mpiexec
```

研究脚本首先重算采用三点支撑的四个完整三维参考。它在实际生产进程外计时，固定
CPU、单线程及 MUMPS，计入全部截面处理、求解、材料点恢复和文件输出。
非均匀载荷采用两轮相反顺序比较 12 个全长模式及三个宽度分区方案，所有这些
方案均不采用端部凝聚；其他载荷运行一轮。三维参考每类载荷运行一轮。
参考文件和可执行文件的校验值随结果保存。比较脚本在计时区间之外运行。

任何一项研究方案不满足原有 1% 精度指标时，研究脚本在完成所有计算并写出
`study.csv` 后返回 1；每个算例的 `accuracy_status` 和 `failed_metrics` 明确记录
结果。求解失败和物理一致性检查失败则立即报错，不能作为离散误差继续统计。
双进程命令只检查两个进程与单进程的完整输出一致性，不代表通过三维精度门槛。

## 已暂停的宽度分条试验结果

原始 14 次生产比较均未通过 1% 精度门槛，记录保存在
[`results/width/study.csv`](../verification/section_modal/results/width/study.csv)。
这些记录对应可执行文件 `175e4538e52bf8bd9df590c2780dca7ea905bf022ba0da42b698523e59e63306`，
不是后续实体端部实现的计时结果。各行同时保存参考文件校验值。

下表为半宽非均匀载荷，两次单核外部计时的平均值；应力和能量误差为百分比。

| 宽度方向布置 | 独立截面方向 | 全局方程，含六个约束 | 时间，秒 | 完整应力 L2 | 能量 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原始模态 | 12 | 2,634 | 1.71 | 87.8506 | 74.5216 |
| 两条轴向线 | 25 | 5,481 | 2.04 | 86.9659 | 72.9319 |
| 三条轴向线 | 38 | 8,328 | 2.38 | 82.6380 | 63.2870 |
| 五条轴向线 | 64 | 14,022 | 3.20 | 80.3291 | 59.5912 |

分条增加了截面空间，完整应力误差有所下降，但这些计算不足以准确表示三点
支撑附近的局部三维变形。五条轴向线时，均布横向载荷超过 99.99% 的应力误差
平方积分位于根部前 6 毫米。内部残量、虚功和并行一致性通过，不等于通过
完整三维精度比较。该分支按用户要求暂停，代码和失败结果保留作为研究记录。
后续[真实三维端部试验](section-solid-ends.md)使用不同的局部运动学空间。

## 学术关系

截面展开与一维轴向有限元组合已有相关研究，例如
[Carrera 等的任意截面位移展开](https://iris.polito.it/handle/11583/3014693)以及
[节点相关截面运动学](https://iris.polito.it/handle/11583/2692893)。这些工作提供了
研究背景，不能证明本例的性能或 1% 精度。本试验评价的是当前 Fuelsim 离散、
种子模式和边界条件下的具体选择。
