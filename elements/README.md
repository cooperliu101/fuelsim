# fuelsim 单元计算库

本目录与 fuelsim 保持在同一个仓库，但能够独立配置、构建和测试。当前提供完整的
CAX4T（四节点轴对称温度—位移耦合单元）局部计算。库只依赖 C++17 标准库和
ADlite，不查找或链接 PETSc、MPI、Exodus，也不包含全局网格或求解器。

## 调用与状态所有权

入口为 `fuelsim::elements::evaluate_cax4t`，声明位于
`include/cax4t.hpp`。接口参考 Abaqus UEL 的职责边界，使用具体
C++ 结构体和普通函数，不实现 Abaqus 的二进制接口或参数兼容层。

`Cax4tInput` 借用材料对象、参考几何、当前节点状态、已接受节点状态及材料历史。
这些对象必须在调用期间有效，函数不会保留引用，也不会修改输入。参考几何可通过
`make_quad4_rz_geometry` 预计算并重复使用。每个单元的局部顺序固定为
`[T0..T3, ur0..ur3, uz0..uz3]`，应力和应变顺序为 `[rr, zz, hoop, rz]`，
其中 `rz` 是张量剪应变。

- `committed_history == nullptr` 选择稳态热弹性计算；通常传入零位移的已接受
  节点状态，关闭热容项。返回历史中的应力有效，其他历史字段不代表瞬态更新。
- 非空历史选择瞬态材料更新，要求时间增量有限且大于零。关闭
  `include_thermal_time_term` 仅关闭热容，不关闭蠕变或塑性积分。
- `time` 为当前计算时刻，`time_step` 为本次材料积分的时间增量，均使用 SI 单位。
- 返回的 `residual` 是局部方程不平衡量，`jacobian` 是按行存储的
  `d(residual)/d(state)`，没有额外负号。未请求矩阵时该数组为零。
- `history` 是试探结果。fuelsim 在整个状态验证成功后决定是否接受；丢弃一次
  调用的结果不影响任何后续调用。函数不执行时间推进或状态提交。
- `stored_heat_rate` 与 `generated_heat_rate` 使用本次残量的热积分测度，满足
  四个温度残量之和等于储热率减体热源发热率。它们是功率，不是时间积分后的能量。
  全局边界热流、接触热流、机械能量汇总和两个半步的诊断合并仍由 fuelsim 管理。
- 非法试探温度、位移或变形几何抛出 `std::domain_error`，交给调用方的线搜索或
  缩小时间步重试机制。非法时间增量、已接受节点数据或请求组合抛出
  `std::invalid_argument`。函数不夹持数值，不通过失败结果返回部分历史。

CAX4T 保留原有四点积分、轴对称选择性体积处理、Hughes-Winget 转动、平均温度
热膨胀和专用热算子。运动学自动微分只使用六个局部变量，本构关系只使用五个
局部变量，再以闭式节点链组装 12×12 矩阵。

## 与 fuelsim 的依赖关系

```text
fuelsim_io / fuelsim_solver
             ↓
        fuelsim_core
             ↓
     fuelsim::elements
             ↓
       adlite::adlite
```

CAX4T 实现在 `src/cax4t.cpp`。共用材料对象、材料函数注册表、坐标类型以及
四节点轴对称几何和运动学也由本库持有，其他 fuelsim 单元继续使用这些共用实现。
已有共用类型继续位于 `fuelsim` 命名空间；新单元入口位于 `fuelsim::elements`。
旧材料头文件已经迁移，未设置转发头文件或类型兼容别名。

本次不迁移其他体单元、普通边界积分、接触搜索或全局网格。其他单元不需要通过
CAX4T 入口计算。库没有 fuelsim_core 的反向链接，也不通过包含父目录源码完成构建。
独立仓库、安装包及版本发布流程留待接口成熟后实施。

## 构建与验证

在 fuelsim 仓库根目录执行以下命令，独立构建只需事先安装 ADlite：

```bash
cmake -S elements -B build-elements \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH=/home/cooper/ai_project/fuelsim-dependencies/adlite-0.2.3
cmake --build build-elements --parallel 4
ctest --test-dir build-elements -j4 --output-on-failure
```

`fuelsim_cax4t_element_tests` 直接链接单元库，检查解析受限热膨胀、热源与热容守恒、
小应变和有限应变下塑性与蠕变同时作用的中心差分方向导数，以及历史不可变和失败
试探后的重复性。同一个测试也注册在 fuelsim 的统一 CTest 中。

完整生产计算仍通过 `fuelsim -i <case.fsi>` 验证，保留 B7、B15 和其他现有算例的
参考数据及验收门槛。单元库的局部测试不能替代全局装配、接触和时间推进测试。

## 本次迁移验证记录

迁移前基线为 `7f760e0`。使用相同 GCC 14.3、Release 配置和 ADlite 0.2.3，
在一个畸变四边形上比较两种应变形式与四种材料分支（热弹性、蠕变、塑性、
蠕变与塑性耦合）。每组检查 12 个残量、144 个切线矩阵条目和 72 个历史数值，
共 1,824 个数值的十六进制浮点输出逐字一致。这是限定局部状态的迁移一致性检查，
不代表对所有输入的逐位一致性证明。

独立构建的单元测试通过。新增耦合方向导数检查同时使用 `1e-5` 和 `3e-6`
两种差分步长，小应变最大相对误差为 `4.95543e-8`，有限应变为 `4.9283e-8`，
均低于 `2e-7` 门槛。较小步长在小应变热残量相减时出现舍入放大，因此未通过
降低验收要求来接受该误差。

fuelsim Release 构建通过，统一 `ctest --test-dir build -j4 --output-on-failure`
为 **272/272 通过**，包含 B7、B15、综合接触、重启动和并行等价性检查。受沙盒
通信初始化限制，完整回归在沙盒外运行。外部程序对比使用已追踪的 Abaqus/MOOSE
参考结果，本次没有重新生成这些参考数据。未进行计算性能提升的评估。
