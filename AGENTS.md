# fuelsim Agent Guide

本文件适用于整个仓库。

## 目标与当前范围

`fuelsim` 使用 C++17 开发核燃料性能有限元程序。当前 M1 是串行 2D
轴对称 RZ 燃料—包壳稳态热弹性求解器：

```text
体单元/界面固定 12 DOF -> ADlite 局部 Jacobian
                       -> PETSc 统一稀疏装配 -> SNES Newton
```

M1 使用两个独立的结构化 Quad4 网格、匹配 Line2 界面、参考构形和
field-major 全局自由度：

```text
[T(:), ur(:), uz(:)]
```

## 依赖与 C++ 约束

- 仅使用 C++17。
- 不定义项目自己的 C++ 类模板、函数模板、表达式模板或标量泛型层。
- 允许使用 `std::vector`、`std::array` 等标准库模板。
- 自动微分只能使用用户的 ADlite 软件包和具体类型
  `adlite::Scalar`。
- 除 ADlite 外，唯一允许的外部数值依赖为 PETSc。
- MPI、BLAS、LAPACK、MUMPS、Hypre 等只能作为 PETSc 的传递依赖，不得由
  fuelsim 单独发现或链接。
- 不引入 Eigen、Boost、fmt、JSON/YAML、CLI、日志或第三方测试框架。
- 测试使用普通 C++ 自检程序和 CTest。
- 公共头文件不得暴露 PETSc 类型；PETSc 保持在 `fuelsim_petsc` 实现层。
- 不使用 `FetchContent` 或构建时网络下载。

## 数值契约

- 局部自由度顺序固定为
  `[T0..T3, ur0..ur3, uz0..uz3]`。
- ADlite 只按 12 个体单元或界面局部自由度播种，禁止按全局自由度播种。
- RZ 积分测度为完整的 `2*pi*r*detJ*w`。
- 应变和应力分量顺序为 `[rr, zz, hoop, rz]`，`rz` 是张量剪应变。
- 历史变量以后使用 `double` 保存；只有 trial state 使用 ADlite。
- PETSc Dirichlet 约束使用 `F_i=x_i-g_i` 和只清行的
  `MatZeroRows(..., diagonal=1)`。
- 当前只支持一个 MPI rank，不得把每个 rank 重复装配全模型称为并行。
- 不隐式夹持异常材料值或几何值；非法结构输入应明确报错。
- 界面间隙为 `g=(Rc+urc)-(Rf+urf)`；开放为正、穿透为负。
- 气隙导热为 `h=k_gap/max(g,g_min)`。
- 法向压力为 `p=penalty*max(-g,0)`，罚参数单位为 `Pa/m`。
- 界面两侧必须使用同一个参考燃料表面测度并保持热量、法向力守恒。

## 架构边界

- `fuelsim_core`：网格、自由度、材料、Quad4 RZ 核和问题定义，仅依赖
  ADlite。
- `fuelsim_petsc`：PETSc 会话、稀疏装配和 SNES 求解。
- `NonlinearProblem` 只作为求解器端口；不得扩张成 MOOSE 式对象工厂。
- M1 燃料和包壳节点必须保持独立，首版只支持匹配轴向界面。
- 不复制 MOOSE 的对象工厂、继承层次或输入参数系统。
- 不复制 jax_fuel 的运行时声明式 Kernel 注册系统。
- 新物理先形成具体、可验证的局部残量，再考虑通用化。
- 除非用户明确要求，不增加旧 API 别名、适配器或兼容层。

## 必须执行的验收

先在 `moose` Conda 环境中安装 ADlite，然后：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/tmp/adlite-fuelsim-install
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

PETSc/MPICH 测试在受限沙盒内可能出现 `OFI EP enable failed`。遇到该错误应在
沙盒外重跑，不能归因于 fuelsim 数值实现。

每次物理修改至少检查：

1. 对应局部 AD Jacobian 与中心差分方向导数；
2. 相关解析解；
3. 默认端到端求解；
4. 匹配物理、罚参数、加载路径和网格设置的 MOOSE 对标量。

构建成功不等于数值验收通过。只有相关 CTest、解析指标和 MOOSE 指标全部
满足门槛后，才能声称功能完成。
