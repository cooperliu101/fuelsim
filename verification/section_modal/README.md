# 截面模态生产求解与三维板比较

本目录所有结构计算均运行 `fuelsim -i <case.fsi>`。完整输入卡显式保存；测试只
逐字复制输入卡及网格。Python 仅生成网格、调度生产程序和读取输出，不装配或
重新求解力学问题。三维参考采用 Fuelsim C3D20T，本次比较不属于新增 Abaqus 鉴定。

## 生产输入和适用范围

在稳态 Cartesian 三维输入中增加以下节即可选择截面模态求解；省略时运行实体模型：

```text
[SectionModes]
 count = 12
 [ends]
  count = 288
  lower_interface = lower_fine
  upper_interface = upper_fine
 []
 [transition]
  count = 192
  lower_interface = lower_interface
  upper_interface = upper_interface
 []
[]
```

`count` 为 3 至 768，且不能超过截面节点数的三倍或可生成的独立候选方向。
自动畸变求解的 512 个横向自由度上限仍保留。生产入口读取一个
沿正 z 方向固定挤出的 HEX20 Exodus 网格，提取共享 QUAD8 截面和轴向端节点。
轴向单元允许非均匀长度，截面不要求矩形规则网格，命名材料块必须沿 z 保持不变。
缺失、重叠、连接不一致或无法构成同一固定截面的单元会明确报错。

本阶段只接受固定温度、无本征应变的小应变线弹性，`Executioner` 使用
`type=steady` 和 `load_steps=1`。求解器使用 PETSc/MUMPS 直接分解；Solver 字段
为 `linear_solver`、`preconditioner`、`direct_factorization`、`absolute_tolerance`、
`relative_tolerance`、`maximum_iterations`。节点集或边集可以约束物理位移，
压力和牵引必须采用参考构形。所有温度约束等于对应区域初温。接触、体力、热载荷、
时间推进、重启动和材料历史更新尚未接入这条生产路径。

截面场依次包括三个经典整体模式、可选扭转、两个独立转角、经典泊松松弛、
独立轴向翘曲，以及由二维平衡方程和畸变特征问题生成的方向。有限谱、非刚体
剪切零空间、轴向平衡修正和横向平衡修正各有独立的数学定义；推导见
[设计说明](../../docs/section-modal-mechanics.md)。候选场经两遍截面 L2 正交化及
线性相关检测，生产 Ritz 方向与离线保存的原始特征向量明确区分。

经典耦合模式在离线接口中永久保留。生产模型在两个独立转角具备后使用混合弯曲，
`gamma_xz=wx'+rx`、`epsilon_zz=X*rx'`，Euler–Bernoulli 子空间为 `rx=-wx'`。
若某个泊松或翘曲导数场已经完全包含在独立位移场中，生产离散采用对应独立物理
幅值；否则保留其原有导数项。不能把高连续性多项式空间中的导数关系当成无害的
代数换元，也不能删除尚无独立表示的翘曲。轴向伸长、两方向纯弯曲、扭转的完整
应变和能量均有独立幅值形式的内部验证。

每个活动模式在轴向节点保存幅值及其一、二阶导数，采用五次 Hermite 插值和六点积分。
不设置子节时，机械未知量数为 `3 * mode_count * axial_nodes`。子节指定端部局部
修正：`lower_interface` 以下和 `upper_interface` 以上的节点使用相应模式数；
可省略其中一个选择器。启用局部凝聚时全长至少需要 10 个方向，以完整保留独立
转角、泊松场和扭转翘曲所构成的经典空间；本次验证固定使用 12 个全长模态。
选择器必须是 Exodus 中完整内部端截面的节点集。
多个区域重叠时取较大的模式数，因此可以向内部逐层缩小局部基。两端局部区域
不能连成覆盖全长的区域，生产代码不推断或硬编码这些区域的位置。

不同区域共用同一截面基的前若干方向。缺失的高阶方向在交界处 q、q'、q'' 均为零。
端部全部内部幅值及相应约束乘子通过局部平衡消去，最后的全局方程只保存内部
及交界节点的少量模态。求解后恢复全部局部场；二维截面节点没有成为全局未知量。
物理位移约束使用带主元的正交化去除重复行，独立约束的乘子数另外报告。
当前会明确拒绝同时跨越局部与保留未知量的单条位移约束。
数学推导及成本边界见 [少量全长模态与端部修正](../../docs/section-modal-efficiency.md)。

## 材料点与输出

生产输入必须提供 `Outputs.csv` 和 `Outputs.history`，模态路径暂不输出 Exodus。
汇总文件记录模式分类数量、幅值自由度、独立约束、能量、外功、物理残量以及
线性残量修正次数。矩阵乘法残量和材料点虚功残量分别保存。`mode_count` 是全长
模式数，`section_basis_size` 和各分类数量描述离线截面基。`dof_count` 是凝聚后
保留的幅值数，`recovered_local_dofs` 是恢复的局部幅值数，`active_axial_dofs`
是两者之和。`global_system_size` 与 `local_system_size` 还计入相应约束乘子。
截面预处理、轴向装配、凝聚、求解与残量修正、物理场恢复的秒数分别记录。

历史文件按 `kind` 区分记录：

- `node`：全部源网格节点的坐标和三个物理位移。
- `point`：全部轴向与截面材料点的坐标、位移、积分权重、六分量应变和应力。
- `section`：轴向积分位置的轴力、两弯矩、扭矩和单位轴长应变能。
- `amplitude`：轴向节点、模式编号、幅值及其一、二阶导数。
- `local_amplitude`：凝聚后恢复的局部节点及其全部活动模态幅值。

应变按 `[xx,yy,zz,xy,yz,xz]` 保存，剪应变是张量分量。相对于材料加权中心，
`N=integral(szz)`、`Mx=integral(Y*szz)`、`My=-integral(X*szz)`。
外功包括给定载荷及非零给定位移的约束反力做功。

固定线弹性切线允许先积分截面材料切线，再与轴向形函数导数乘积收缩；这一积分
顺序与完整材料点切线逐项核对。它不替代材料点响应：每个轴向位置仍重构三维应变，
调用材料接口，独立计算应力、残量和能量。固定的截面运动学系数可以复用，复用
与逐次计算的应变、应力、残量、切线及能量有逐项完全一致性检查。
细长板求解使用这些物理残量修正直接
解，复用原有分解，不添加人工刚度。未来非线性切线应由已有的逐材料点入口更新。

## 同条件三维比较

几何为 H=0.002 m、W=0.02 m、L=0.2 m，E=70 GPa、nu=0.3、T=300 K。
根部 z=0 的三个物理位移均为零，模态模型使用相同的完整夹持约束。

| 算例 | 载荷 |
| --- | --- |
| axial | 端面轴向牵引 1 MPa，合力 40 N |
| bending | 端面 x 正负两半的轴向牵引为正负 10 kPa，纯力偶 My=-0.0002 N m |
| transverse | x 正面均布向内压力 10 Pa，线载荷 -0.2 N/m |
| nonuniform | x 正面、y 非负半宽施加 x 向牵引 -20 Pa，同时引起弯曲和扭转 |

弯曲端面牵引为两段常数，其局部应力也参与比较，没有替换成解析线性应力。

保留的粗网格 `plate.e` 为 2/4/8 个 HEX20，用于低模态生产路径和约束回归。
严格精度检查使用以下网格：

| 模型 | x/y/z 单元数 | 机械未知量或积分点 |
| --- | --- | --- |
| 模态输入 `plate_local.e` | 4/8/72 | 121 个截面节点，73 个轴向端节点；124416 个材料点 |
| 三维参考 `plate_reference.e` | 4/8/1024 | 510315 个位移自由度；884736 个材料点 |
| 三维加密检查 `plate_reference_check.e` | 4/8/512 | 255339 个位移自由度；442368 个材料点 |

模态轴向网格在两端加密：`z(s)=L*(s-0.85*sin(2*pi*s)/(2*pi))`，s 均匀划分。
HEX20 轴向边中点取两个端点的算术平均，因此比较单元仍为仿射挤出。
所有几何和离散规模均来自提交的 Exodus 网格；生产代码不生成这些网格。
手动全长扫描的 `plate_refined.e` 保留 128 层；`plate_local.e` 保留该网格两端各八层，
在其余平滑区间每隔一个截面平面取一个，形成 72 层。它保持几何和二维截面连接，
并增加距两端约 6 mm、30 mm 处的四个截面节点集。12 个模态覆盖全长，局部使用
288、192 个方向；这些高阶方向不进入内部全长方程。精度与成本须同时查看，不能
仅用凝聚后方程大小推断加速。

比较覆盖三维参考的全部位移节点和全部原生材料点，包括夹持端和受载端。
模态点输出通过 3×3×6 多项式插值重构到参考位置，位移还检查单元公共面连续性，
并与独立的源节点输出核对。该插值仅用于这里的仿射、固定均匀线弹性验证网格。

严格小于 1.0% 的指标为：全节点位移向量相对 L2、按体积加权的轴向应力相对 L2、
完整应力张量相对 L2、两方向弯矩相对 L2，以及总应变能相对误差。应力张量内积
使用 `[1,1,1,2,2,2]`；弯矩在全部三维轴向积分平面计算。这个门槛不是最大逐点
误差门槛。轴向拉伸的零弯矩采用 1e-10 N m 的绝对误差检查，不设相对分母下限。

额外检查物理平衡残量不超过 1e-7、物理约束误差不超过 1e-9、外功与两倍能量
的相对差不超过 1e-8，并重新积分点应力核对生产截面内力。误差 CSV 在判定前写出，
失败结果不会被丢弃。原始记录见 `results/`。

12 个全长模态、32256 个端部局部幅值和 1044 个最终全局幅值的四类生产比较全部
通过。完整应力张量相对 L2 误差依次为 0.248804%、0.400454%、0.835404%、0.872976%。
每项位移、轴向应力、弯矩和能量结果见 `results/*_local_comparison.csv`，规模、
物理残量和内部各阶段计时见 `results/*_local_summary.csv`。这些内部计时来自
回归运行；受控性能比较单独使用 `benchmark.py` 的外部计时。

三维参考从 512 层到 1024 层的独立变化为：共同节点位移 0.00533%，轴向应力
0.407%，完整应力 0.863%，弯矩 0.169%，能量 0.00587%。1024 层参考弯矩相对
载荷平衡解析式 `My(z)=-0.2*(L-z)^2/2` 的误差为 0.0465%。这个检查只加密轴向，
不构成截面网格或夹持角点逐点应力已经收敛的证明。

## 复现

```bash
build/fuelsim -i verification/section_modal/transverse_local_12.fsi
ctest --test-dir build -j8 -R 'fuelsim_(section_modal_|modal_beam_tests|classic_section_tests|section_warping_tests|section_distortion_tests)' --output-on-failure
ctest --test-dir build -j8 --output-on-failure

# 固定单个 CPU、单线程、相同 MUMPS 的外部计时；包含预处理、求解和结果输出。
python verification/section_modal/benchmark.py build/fuelsim build/blackbox/section_modal_benchmark --cpu 0

# 手动三维参考检查。网格程序不生成或修改输入卡。
python verification/section_modal/mesh.py
build/fuelsim -i verification/section_modal/transverse_mesh_check.fsi
python verification/section_modal/refinement.py \
  verification/section_modal/transverse_mesh_check_solid.e \
  build/blackbox/section_modal_accuracy_transverse/transverse_refined_solid.e \
  verification/section_modal/reference_refinement.csv
```

验证读取需要 Python、NumPy 和 netCDF4，生产库没有新增数值依赖。
CMake 的 `FUELSIM_SECTION_PYTHON` 可以指定解释器。精度测试预留八个 CTest 调度
槽，使大型分解在 `-j8` 回归中逐个运行；其他测试仍并发，实际每个串行生产任务
使用一个线程。形成约束矩阵后立即释放原始刚度矩阵，避免分解期间保留两份。
运行日志直接写入文件，即使任务被系统终止也能保留已经输出的诊断。
双进程验证比较全部位移、材料点应力、积分能量及输出位置，要求相对差小于 1e-8。
重复网格和模态扫描保留为手动材料，不加入自动测试。
