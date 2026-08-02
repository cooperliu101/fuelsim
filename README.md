# fuelsim

`fuelsim` 是一个依赖精简的 C++17 核燃料性能有限元程序。统一可执行程序
通过输入卡选择 `SteadyProblem` 或 `TransientProblem`，并从一个 Exodus
文件自由组合任意数量的物理区域。稳态问题包含：

- 一个 `.e` 文件中的多个独立 2D 轴对称 RZ Quad4 区域，区域间节点不合并；
- 每个区域独立的稳态温度相关热传导、体积热源和小应变热弹性；
- 由 `primary`、`secondary` 边集定义的 STS 气隙导热；
- 与 MOOSE/JAX 实现一致的 secondary 节点到 primary 线段 NTS 无摩擦罚接触；
- ADlite 生成体单元和界面的局部 Jacobian；
- PETSc SNES、KSP 和 AIJ 稀疏矩阵完成串行 Newton 求解。

瞬态问题在同一接触离散上实现：

- M2.1：一致热容矩阵、Backward Euler 物理时间积分、热源斜坡、
  committed/trial/commit/rollback 和失败步 cutback；
- M2.2：通用小应变 J2 Norton 蠕变、J2 线性各向同性硬化塑性及同一
  材料点的全隐式耦合；
- M2.3：芯块热膨胀闭合初始间隙的 PCMI 回归，包壳同时累积塑性与
  蠕变历史，并与同网格、同时间步 MOOSE 算例比较；
- 燃料、包壳各自独立的积分点 `double` 历史，Newton 回调内只生成
  ADlite trial 状态；
- 瞬态热传导与准静态力学耦合，PETSc 工作区跨时间步复用。
- 版本化二进制检查点保存完整 committed 状态，并通过模型签名、字节序、
  长度和校验和执行严格重启动；
- Exodus 场结果包含节点温度/位移/接触量、四积分点应力与非弹性历史，以及
  守恒界面总量。

M2 材料参数仅用于算法和软件验证，不代表真实燃料或包壳经验模型。
具体状态契约、算法和边界见 [M2 设计说明](docs/m2.md)。

M0 单燃料圆柱仍作为解析解和 MOOSE 回归基线保留。项目不定义自己的
C++ 模板，不引入 Eigen、Boost、JSON/YAML、日志库或第三方测试框架。
除 ADlite 外，PETSc 是唯一直接外部数值依赖，Exodus 是唯一直接网格 I/O
依赖；MPI、BLAS 和直接求解器等只允许作为 PETSc 的传递依赖。

## 构建

先安装 ADlite 到临时前缀：

```bash
source /home/cooper/miniforge/etc/profile.d/conda.sh
conda activate moose

cmake -S /home/cooper/ai_project/ADlite \
  -B /tmp/adlite-fuelsim-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DADLITE_BUILD_EXAMPLES=OFF
cmake --build /tmp/adlite-fuelsim-build --parallel
cmake --install /tmp/adlite-fuelsim-build \
  --prefix /tmp/adlite-fuelsim-install
```

### 主开发入口：旧 PETSc + 直接 Exodus API

PETSc 只负责求解，不需要启用 Exodus。Exodus 单独构建为串行 I/O 库并复用
`moose` Conda 环境中的 NetCDF；fuelsim 不使用 DMPlex，也不使用
`PetscViewerExodusII`。先安装独立 Exodus：

```bash
git clone --branch v2024-06-27 --depth 1 \
  https://github.com/gsjaardema/seacas.git /tmp/seacas-exodus-src

./scripts/build_exodus.sh \
  /tmp/seacas-exodus-src \
  /home/cooper/.local/exodus-2024-06-27 \
  /home/cooper/miniforge/envs/moose
```

然后使用原 Conda PETSc 构建 fuelsim：

```bash
env \
  PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  PKG_CONFIG_PATH=/home/cooper/miniforge/envs/moose/lib/pkgconfig \
  cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
    -DCMAKE_PREFIX_PATH=/tmp/adlite-fuelsim-install \
    -DSEACASExodus_DIR=/home/cooper/.local/exodus-2024-06-27/lib/cmake/SEACASExodus

env PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  cmake --build build --parallel

ctest --test-dir build --output-on-failure
```

`fuelsim_exodus` 将一个 2D 非结构 Quad4 文件的节点、连接关系、元素块、
节点集和边集转换为 fuelsim 自有网格对象；专项 CTest 通过 Exodus API
写出并逐项回读。每个 `[Regions/<name>]` 选择一个元素块；接触不再重复声明
块名，而是由边集相邻单元自动确定所属区域。这样同一入口既能分析单独芯块
或包壳，也能组合芯块—包壳以及芯块—包壳1—包壳2。所有 MOOSE 对比都会
读取仓库内由对应 MOOSE 输入生成的 `*_mesh.e`。一般非结构 Quad4 可进入
I/O 层，生产问题会保留每个选中块的原始节点坐标和 Quad4 连接关系，不再
重建张量积 RZ 网格。单元仍须具有有效的正 Jacobian；当前接触实现另外要求
参与接触的边集是轴对称圆柱面。直接 Exodus I/O 为串行操作；当前求解器也
只允许一个 MPI rank。

运行稳态燃料—包壳工况：

```bash
./build/fuelsim -i verification/fuelsim/steady_fuel_cladding.fsi
```

运行带包壳塑性—蠕变耦合的瞬态 PCMI 工况：

```bash
./build/fuelsim -i \
  verification/fuelsim/transient_fuel_cladding_pcmi.fsi
```

输入卡采用严格、带版本号的 MOOSE 风格分段文本；未知段、未知键、重复键和
不适用于所选本构模型的参数都会立即报错。网格尺寸和几何只从 Exodus 文件
读取，不在输入卡中重复维护。完整字段说明见
[输入卡说明](docs/input-card.md)。

瞬态工程输出示例：

```text
[Executioner]
  type = transient
  # 其余时间推进字段省略
  restart = previous.checkpoint
[]

[Outputs]
  exodus = fields.e
  checkpoint = latest.checkpoint
  checkpoint_interval = 5
[]
```

检查点只包含最后一个成功提交态，不包含 Newton trial 或失败时间步。严格
状态与文件契约见 [M3 工程化说明](docs/m3.md)。

运行 M2.3 PCMI—MOOSE 验收：

```bash
./build/fuelsim_m2_pcmi_solver_tests \
  verification/moose/m23_pcmi_coupled_cladding_rz_mesh.e
```

运行 MOOSE Exodus 网格驱动的 M1 验收：

```bash
./build/fuelsim_m1_exodus_moose_tests \
  verification/moose/m1_fuel_cladding_gap_rz_mesh.e
```

运行内部节点畸变、不可转换为张量积网格的 M1 全场验收：

```bash
./build/fuelsim_m1_unstructured_moose_tests \
  verification/fuelsim/steady_fuel_cladding_unstructured.fsi \
  verification/moose/m1_fuel_cladding_unstructured_rz_all_nodes_final.csv \
  verification/moose/m1_fuel_cladding_unstructured_rz_fuel_surface_final.csv
```

稳态示例将热源分成 20 个线性载荷步，以稳定跨越接触活动集的切换。
PETSc 选项仍可在命令行覆盖，例如：

```bash
./build/fuelsim \
  -i verification/fuelsim/steady_fuel_cladding.fsi \
  -snes_monitor -ksp_error_if_not_converged
```

当前只支持一个 MPI rank。

程序会同时输出问题构造、PETSc 设置、非线性求解、残量回调和 Jacobian
回调的内部计时。20 个载荷步复用同一问题几何、SNES、Vec、Mat、矩阵非零
结构和回调缓冲区；热源只更新各区域的具体核参数。

M2 同样复用几何和 PETSc 工作区，但每个成功时间步会提交 nodal 温度/位移
以及积分点非弹性历史。失败尝试从最后一个 committed 状态重启，不会把
未收敛 Newton 状态作为下一次初值。

## 数值契约

全局状态采用统一的 field-major 排列：

```text
[T(:), ur(:), uz(:)]
```

每个体单元及界面贡献均固定为 12 个局部自由度，因此 PETSc 只保留一条
装配路径。体单元顺序为：

```text
[T0..T3, ur0..ur3, uz0..uz3]
```

热界面的四个节点按
`[secondary0, secondary1, primary0, primary1]` 排列，对应：

```text
[Ts0, Ts1, Tp0, Tp1,
 urs0, urs1, urp0, urp1,
 uzs0, uzs1, uzp0, uzp1]
```

当前间隙和界面定律为：

```text
g = (Rp + urp) - (Rs + urs)
h = gap_conductivity / max(g, minimum_gap)
q = h * (Ts - Tp)
p = contact_penalty * max(-g, 0)
```

`q>0` 表示热量由 secondary 流向 primary，`p>0` 表示压缩接触压力。
`contact_penalty` 的单位为 `Pa/m`。热接触在当前 secondary 表面
`2*pi*r*J` 上积分，并把相反热流投影到 primary 节点。每个 secondary
节点只有一个有效机械投影，节点反力按当前 secondary 半边面积集总后，通过
primary 线段形函数分配相反反力。两种界面残量均离散守恒，投影、面积和
界面定律都由 ADlite 线性化。

默认燃料高度为 `10.000 mm`，包壳高度为 `10.020 mm`。顶部额外的
`20 um` 轴向裕量用于防止燃料热膨胀后越过包壳接触面。

PETSc Dirichlet 行采用 `F_i=x_i-g_i`，Jacobian 只清约束行并置单位对角。
默认使用 Newton basic line search，用户传入的 PETSc 选项可覆盖默认值。

M2 热容残量和非弹性更新为：

```text
Rcap_i = integral(N_i*rho*cp*(Tnew-Told)/dt*2*pi*r*dA)
q_creep + 3*G*dt*A*(q_creep/q_ref)^n = q_trial
delta_ep = max((q_trial-yield_old)/(3*G+H), 0)

q_trial = q + 3*G*(delta_ep+delta_ec)
delta_ec = dt*A*(q/q_ref)^n
q = yield_old + H*delta_ep               coupled active branch
```

材料点可选择 `elastic`、`norton_creep`、`j2_plasticity` 或
`norton_creep_j2_plasticity`。耦合分支先判断蠕变松弛后是否仍然超出
屈服面，再同时求解塑性和蠕变增量。

## 当前验收

CTest 覆盖：

- 严格输入语法、物理问题调度、未知键和本构条件字段拒绝；
- 单区域、三层区域和两对 `primary/secondary` 接触的自由组合；
- 稳态与瞬态输入卡读取 MOOSE Exodus 网格的端到端求解；
- RZ 体积积分、形函数和梯度恒等式；
- 体单元 AD Jacobian 与中心差分方向导数；
- 开放、最小热隙饱和和闭合接触三种界面分支的 AD Jacobian；
- STS 热流和 NTS 节点反力守恒；
- 接触投影半开区间的唯一性、越界失活和加高包壳覆盖；
- 带热源实心圆柱解析温度；
- 无应力自由热膨胀；
- 开口端 Lamé 厚壁圆筒压力解；
- 开放气隙燃料—包壳圆柱的解析热阻解；
- 闭合气隙的 20 步端到端求解；
- 匹配物理与网格设置的 MOOSE 全场温度、径向/轴向位移、接触压力分布及
  总反力对比；非零参考场统一检查相对 L2、相对绝对峰值和最大逐点相对
  误差；
- M0、M1、M2.1、M2.2 四条材料路径和 M2.3 PCMI 全部直接读取各自 MOOSE
  输入生成的 Exodus 网格，核对块、节点集和边集后再执行对比；
- M2 瞬态体单元 AD Jacobian、一致热容矩阵和均匀绝热升温解析解；
- J2 塑性闭式径向返回、卸载和活跃分支 AD 切线；
- Norton `n=1` 解析根、非线性局部残量和活跃分支 AD 切线；
- 塑性—蠕变耦合应力平衡、屈服一致性、退化分支、极端尺度和 AD 切线；
- M2 committed 状态不受残量/Jacobian 回调影响，接受步提交且拒绝步回滚；
- 两步 M2 PETSc 工作区复用和 20 步同点非零蠕变/塑性燃料—包壳演示；
- MOOSE 瞬态热容、Norton、J2 及两套耦合载荷路径参考。
- 芯块热膨胀闭合 1 um 间隙、5 个节点接触、包壳塑性—蠕变耦合的
  20 步 PCMI—MOOSE 端到端回归。

独立 MOOSE 输入、结果快照和运行条件位于
`verification/moose/`。

默认 M1 最终步对全部 528 个节点以及 11 个燃料接触节点进行比较。温度、
径向位移、轴向位移和接触压力的相对 L2 误差分别为 `0.00989%`、
`0.06382%`、`0.01774%` 和 `0.22071%`；相对绝对峰值误差分别为
`0.00055%`、`0.12765%`、`0.02530%` 和 `0.33282%`；最大逐点相对误差
分别为 `0.08710%`、`0.38834%`、`0.75155%` 和 `0.33282%`，十二项均
低于 `1%`。总接触反力的三项单值误差为 `0.0443%`，11 个燃料表面节点
均成功投影且处于接触状态。

内部节点畸变的非张量 M1 网格保留全部原始 Exodus 连接关系。与 MOOSE
逐节点全场比较时，温度、径向位移、轴向位移和接触压力的相对 L2 误差分别
为 `0.00990%`、`0.06393%`、`0.01772%`、`0.21727%`；同一顺序下，相对
绝对峰值误差分别为 `0.00064%`、`0.12758%`、
`0.02529%`、`0.32887%`，最大逐点相对误差分别为 `0.08697%`、
`0.38810%`、`0.75970%`、`0.32887%`，十二项均低于 `1%`。逐点相对误差
不对精确零参考值做除法；径向和轴向位移分别有 11 和 48 个零参考点，其
最大绝对差均为 `0`。

M2 的 MOOSE 最小参考同样比较全部网格节点的三个场。均匀瞬态升温的温度
三项相对误差均为 `0`；径向和轴向位移参考场全为零，三项绝对误差也均为
`0`。Norton 的应力、
等效蠕变和位移误差分别为 `0.0019924%`、`0.0089635%` 和
`0.0036523%`；J2 应力和等效塑性应变与解析值在输出精度内一致。耦合
牵引路径的应力、等效塑性、等效蠕变和位移误差分别为
`0.00000394%`、`0.001188%`、`0.000472%` 和 `0.000408%`；位移控制路径
四项误差也均小于 `0.001%`。

M2.3 PCMI 算例采用弹性芯块和耦合 Norton—J2 包壳。20 s 末，5 个芯块
表面节点均投影并接触。全部 50 个节点的温度三项误差最大为 `0.04522%`；
轴向位移三项误差最大为 `0.07670%`。径向位移的相对 L2、相对绝对峰值和
最大逐点相对误差分别为 `0.08796%`、`0.11889%` 和 `0.25280%`；后两项
的局部差异来自顶部包壳内表面节点，最大绝对差约 `6.45 nm`，因此径向
位移采用全场 L2 `<0.1%`、另外两项 `<0.5%` 的回归门槛。包壳平均等效
应力误差为 `0.00431%`，平均等效塑性和蠕变应变误差分别为 `0.0489%` 和
`0.0359%`。接触压力三项误差分别为 `0.06243%`、`0.10375%` 和
`0.10375%`，总接触力误差为 `0.00027%`。

包壳 8 个单元的全部 32 个 Gauss 点也与 MOOSE 一一对齐。等效应力、
等效塑性应变和等效蠕变应变的逐点相对 L2 误差分别为 `0.0124%`、
`0.1389%` 和 `0.0769%`，相对绝对峰值误差分别为 `0.0334%`、`0.2756%`
和 `0.1884%`，最大单点相对误差分别为 `0.0334%`、`0.2920%` 和
`0.1884%`。逐点回归门槛为相对 L2 `<0.2%`，相对绝对峰值与最大逐点
相对误差均 `<0.5%`。

## 单核性能

固定 CPU 和所有数值库线程数为 1 后，默认 1,584 DOF 算例的五次运行中位数
为 fuelsim `2.00 s`、MOOSE `3.21 s`，fuelsim 约快 `1.61` 倍。与提交
`319e173` 的同机配对基线相比，对象复用使 fuelsim 从 `2.23 s` 降至
`1.95 s`，墙钟时间缩短 `12.6%`。

仓库另提供不纳入 CTest 的 23,010 DOF、20 载荷步手动基准：

```bash
./build/fuelsim_m1_single_core_benchmark
```

两次单核运行的均值为 fuelsim `40.45 s`、MOOSE `54.97 s`，fuelsim 约快
`1.36` 倍。完整网格、命令、观测范围和计时口径见
`benchmarks/README.md`。

## 当前边界

首版使用相同轴向单元数，但允许包壳略高于燃料；每个燃料热边的两个
Gauss 点仍必须投影到同一个包壳线段。不支持通用非匹配 mortar、摩擦、
位移惯性或有限应变。

M2.2 只提供与具体材料无关的等温 Norton 幂律和线性硬化 J2 算法。暂不
实现温度、辐照、燃耗、孔隙率或应力相关的真实燃料/包壳经验关联，也不
支持有限应变、非共轴多机制、裂变气体、燃料重定位或开裂。后续真实模型
必须在当前状态事务、局部 AD Jacobian 和独立 MOOSE 回归基础上逐项加入。
