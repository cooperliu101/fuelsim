# fuelsim

`fuelsim` 是一个依赖精简的 C++17 核燃料性能有限元程序。默认 M1
算例包含：

- 两个独立的 2D 轴对称 RZ Quad4 网格：燃料与包壳节点不合并；
- 稳态温度相关热传导、燃料体积热源和小应变热弹性；
- 燃料侧 STS 积分并投影到包壳线段的气隙导热；
- 与 MOOSE/JAX 实现一致的燃料节点到包壳线段 NTS 无摩擦罚接触；
- ADlite 生成体单元和界面的局部 Jacobian；
- PETSc SNES、KSP 和 AIJ 稀疏矩阵完成串行 Newton 求解。

M0 单燃料圆柱仍作为解析解和 MOOSE 回归基线保留。项目不定义自己的
C++ 模板，不引入 Eigen、Boost、JSON/YAML、日志库或第三方测试框架。
除 ADlite 外，唯一直接外部数值依赖是 PETSc；MPI、BLAS 和直接求解器等
只允许作为 PETSc 的传递依赖。

## 构建

本机 PETSc 位于 `moose` Conda 环境。先安装 ADlite 到临时前缀：

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

随后构建和测试 `fuelsim`：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/tmp/adlite-fuelsim-install
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

运行默认 M1 工况：

```bash
./build/fuelsim
```

默认工况将热源分成 20 个线性载荷步，以稳定跨越接触活动集的切换。
PETSc 选项仍可在命令行覆盖，例如：

```bash
./build/fuelsim -snes_monitor -ksp_error_if_not_converged
```

当前只支持一个 MPI rank。

默认程序会同时输出问题构造、PETSc 设置、非线性求解、残量回调和 Jacobian
回调的内部计时。20 个载荷步复用同一个 M1 几何、SNES、Vec、Mat、矩阵非零
结构和回调缓冲区；热源只更新具体的燃料核参数。

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
`[fuel0, fuel1, cladding0, cladding1]` 排列，对应：

```text
[Tf0, Tf1, Tc0, Tc1,
 urf0, urf1, urc0, urc1,
 uzf0, uzf1, uzc0, uzc1]
```

当前间隙和界面定律为：

```text
g = (Rc + urc) - (Rf + urf)
h = gap_conductivity / max(g, minimum_gap)
q = h * (Tf - Tc)
p = contact_penalty * max(-g, 0)
```

`q>0` 表示热量由燃料流向包壳，`p>0` 表示压缩接触压力。
`contact_penalty` 的单位为 `Pa/m`。热接触在当前燃料表面
`2*pi*r*J` 上积分，并把相反热流投影到包壳节点。机械接触以燃料表面节点
为 secondary、包壳线段为 primary；每个燃料节点只有一个有效投影，节点
反力按当前燃料半边面积集总后，通过包壳线段形函数分配相反反力。两种界面
残量均离散守恒，投影、面积和界面定律都由 ADlite 线性化。

默认燃料高度为 `10.000 mm`，包壳高度为 `10.020 mm`。顶部额外的
`20 um` 轴向裕量用于防止燃料热膨胀后越过包壳接触面。

PETSc Dirichlet 行采用 `F_i=x_i-g_i`，Jacobian 只清约束行并置单位对角。
默认使用 Newton basic line search，用户传入的 PETSc 选项可覆盖默认值。

## 当前验收

CTest 覆盖：

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
- 匹配物理与网格设置的 MOOSE 温度、位移、接触压力分布及总反力对比。

独立 MOOSE 输入、结果快照和运行条件位于
`verification/moose/`。

默认 M1 最终步与 MOOSE 的接触压力相对 L2 误差为 `0.2207%`，最大节点
相对误差为 `0.3328%`，总接触反力相对误差为 `0.0443%`；11 个燃料表面
节点均成功投影且处于接触状态。

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

## M1 边界

首版使用相同轴向单元数，但允许包壳略高于燃料；每个燃料热边的两个
Gauss 点仍必须投影到同一个包壳线段。不支持通用非匹配 mortar、摩擦、
蠕变、塑性、裂变气体、燃料重定位或历史变量。这些物理应在当前局部残量、
AD Jacobian 和回归验收稳定后逐项加入。
