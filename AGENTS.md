# fuelsim Agent Guide

本文件适用于整个仓库。

## 目标与当前范围

`fuelsim` 使用 C++17 开发核燃料性能有限元程序。当前 M1 是串行 2D
轴对称 RZ 燃料—包壳稳态热弹性求解器：

```text
体单元/界面固定 12 DOF -> ADlite 局部 Jacobian
                       -> PETSc 统一稀疏装配 -> SNES Newton
```

M1 使用两个独立的结构化 Quad4 网格、燃料侧 STS 热接触、燃料节点到
包壳线段的唯一 NTS 机械接触和 field-major 全局自由度：

```text
[T(:), ur(:), uz(:)]
```

M2.1 在相同空间离散上增加 Backward Euler 一致热容、物理时间步和
committed/trial/commit/rollback。M2.2 增加通用 J2 Norton 蠕变、J2
线性硬化塑性及两者在同一材料点的全隐式耦合；当前不包含真实燃料或包壳
经验模型。

## 依赖与 C++ 约束

- 仅使用 C++17。
- 不定义项目自己的 C++ 类模板、函数模板、表达式模板或标量泛型层。
- 允许使用 `std::vector`、`std::array` 等标准库模板。
- 自动微分只能使用用户的 ADlite 软件包和具体类型
  `adlite::Scalar`。
- 除 ADlite 外，唯一允许的外部数值依赖为 PETSc；Exodus 仅作为直接网格
  I/O 依赖。
- MPI、BLAS、LAPACK、MUMPS、Hypre 等只能作为 PETSc 的传递依赖，不得由
  fuelsim 单独发现或链接。
- 不引入 Eigen、Boost、fmt、JSON/YAML、CLI、日志或第三方测试框架。
- 测试使用普通 C++ 自检程序和 CTest。
- 公共头文件不得暴露 PETSc 类型；PETSc 保持在 `fuelsim_petsc` 实现层。
- 不使用 `FetchContent` 或构建时网络下载。
- 类成员变量统一采用 MOOSE 风格的前缀下划线，如 `_parameters`；不得使用
  `parameters_` 后缀。

## 数值契约

- 局部自由度顺序固定为
  `[T0..T3, ur0..ur3, uz0..uz3]`。
- ADlite 只按 12 个体单元或界面局部自由度播种，禁止按全局自由度播种。
- RZ 积分测度为完整的 `2*pi*r*detJ*w`。
- 应变和应力分量顺序为 `[rr, zz, hoop, rz]`，`rz` 是张量剪应变。
- 历史变量使用 `double` 保存；只有 trial state 使用 ADlite。
- M2 热容使用参考构形一致质量矩阵：
  `N_i*rho*cp*(T_new-T_old)/dt`，不包含位移惯性。
- M2 的所有 Newton、线搜索和失败重试必须从同一 committed 积分点状态
  重算 trial；只能在最终收敛解上重算一次并提交。
- 失败时间步必须同时回滚 nodal state、积分点历史、物理时间和热源；失败
  `SolveResult.state` 不得作为重试初值。
- 轴对称 J2 内积的 `rz` 项必须乘 2；塑性与蠕变应变增量必须无迹。
- Norton 首版为 `rate=A*(q/q_ref)^n` 的等温后向 Euler 更新；J2 首版为
  线性各向同性硬化径向返回。
- 耦合分支必须先求 creep-only 松弛应力；只有该应力仍超过当前屈服应力
  时才激活塑性。活跃时塑性和蠕变共用最终 J2 应力方向并同时满足应力平衡、
  Norton 后向 Euler 方程和塑性一致性条件，禁止一次性算子分裂。
- Norton 与耦合标量根使用对数域保护，等效应力使用 ADlite `hypot`；
  不得通过直接平方和或显式构造超范围幂律系数破坏极端尺度。
- PETSc Dirichlet 约束使用 `F_i=x_i-g_i` 和只清行的
  `MatZeroRows(..., diagonal=1)`。
- 当前只支持一个 MPI rank，不得把每个 rank 重复装配全模型称为并行。
- 不隐式夹持异常材料值或几何值；非法结构输入应明确报错。
- 界面间隙为 `g=(Rc+urc)-(Rf+urf)`；开放为正、穿透为负。
- 气隙导热为 `h=k_gap/max(g,g_min)`。
- 法向压力为 `p=penalty*max(-g,0)`，罚参数单位为 `Pa/m`。
- 热接触在燃料侧当前表面测度上积分，并将相反热流投影到包壳节点。
- 机械接触采用唯一 NTS 投影；燃料节点反力按当前燃料半边面积集总，并按
  包壳线段形函数分配相反反力。
- 热接触与机械接触都必须离散守恒。
- 一个 M1 载荷路径只能构造一次问题几何，并在所有载荷步复用同一组
  SNES、Vec、Mat、非零结构和回调缓冲区；载荷步只更新具体热源参数。
- 内部计时使用单调时钟，至少区分问题构造、求解器设置、非线性求解、残量
  回调和 Jacobian 回调，并记录回调次数及 PETSc 工作区构造次数。

## 架构边界

- `fuelsim_core`：网格、自由度、材料、Quad4 RZ 核和问题定义，仅依赖
  ADlite。
- `fuelsim_exodus`：使用 Exodus API 在 `.e` 文件和 fuelsim 自有非结构
  Quad4 网格之间转换，保留元素块、节点集和边集的 ID 与名称；不使用
  DMPlex，不暴露 Exodus 类型。
- `fuelsim_petsc`：PETSc 会话、稀疏装配和 SNES 求解。
- `NonlinearProblem` 只作为求解器端口；不得扩张成 MOOSE 式对象工厂。
- M1 燃料和包壳节点必须保持独立；默认包壳高度比芯块高 `20 um`，界面通过
  轴向投影耦合。
- 不复制 MOOSE 的对象工厂、继承层次或输入参数系统。
- 不复制 jax_fuel 的运行时声明式 Kernel 注册系统。
- 新物理先形成具体、可验证的局部残量，再考虑通用化。
- M2 使用具体的 `M2Problem`、`M2TimeStepper` 和
  `Quad4RzTransientKernel`；不得把时间状态职责塞入 PETSc 回调。
- 不增加材料工厂或标量泛型层；真实燃料/包壳关联应在通用状态事务稳定后
  作为单独里程碑实现。
- 除非用户明确要求，不增加旧 API 别名、适配器或兼容层。

## 必须执行的验收

后续开发使用原 `moose` Conda PETSc，并直接链接独立的串行 Exodus I/O
库。PETSc 不需要启用 Exodus；不得使用 DMPlex 或 PETSc Exodus viewer。
先安装 ADlite 和 Exodus，然后配置 fuelsim：

```bash
env \
  PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  PKG_CONFIG_PATH=/home/cooper/miniforge/envs/moose/lib/pkgconfig \
  cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH=/tmp/adlite-fuelsim-install \
  -DSEACASExodus_DIR=/home/cooper/.local/exodus-2024-06-27/lib/cmake/SEACASExodus
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

涉及 CMake、PETSc 求解层或依赖配置的修改必须完成上述入口的完整回归。

PETSc/MPICH 测试在受限沙盒内可能出现 `OFI EP enable failed`。遇到该错误应在
沙盒外重跑，不能归因于 fuelsim 数值实现。

每次物理修改至少检查：

1. 对应局部 AD Jacobian 与中心差分方向导数；
2. 相关解析解；
3. 默认端到端求解；
4. 匹配物理、罚参数、加载路径和网格设置的 MOOSE 对标量。

M1 的 Exodus 接入还必须运行 `fuelsim_m1_exodus_moose_tests`，读取
`verification/moose/m1_fuel_cladding_gap_rz_mesh.e`，不得在测试内用
fuelsim 重新生成等价网格。当前 M1 允许从一般非结构 Quad4 文件读取元数据，
但转换到求解网格时，每个 `fuel`/`clad` 块必须是完整张量积 RZ 网格。

M2 还必须检查：

1. residual/Jacobian 重复调用不修改 committed 历史；
2. 接受步只提交一次，拒绝步完整回滚；
3. 活跃蠕变、塑性和耦合分支的 AD 切线分别通过中心差分；
4. 多时间步只构造一次 PETSc 工作区；
5. MOOSE 瞬态温度、应力、位移、等效塑性应变和等效蠕变应变误差均小于
   `0.1%`。

性能修改还必须检查：

1. 默认 1,584 DOF 算例与修改前提交的同机配对计时；
2. 固定 CPU 且 OMP/OpenBLAS/MKL/NUMEXPR 均为 1 线程；
3. `benchmarks/` 中 23,010 DOF、20 步算例至少完成一次；
4. MOOSE 使用相同网格、物理、载荷步、直接求解器并关闭文件输出。

构建成功不等于数值验收通过。只有相关 CTest、解析指标和 MOOSE 指标全部
满足门槛后，才能声称功能完成。
