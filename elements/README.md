# fuelsim 单元计算库

本目录与 fuelsim 保持在同一个仓库，但能够独立配置、构建和测试。提供轴对称与三维体单元、边界积分和局部接触计算。
库只依赖 C++17 标准库和
ADlite，不查找或链接 PETSc、MPI、Exodus，也不包含全局网格或求解器。

## 单元与文件组织

公共头文件直接放在 `include/`，实现放在 `src/`，局部测试放在 `tests/`。
`include/detail/` 只保存内部共用计算辅助代码。

| 范围 | 头文件 | 主要实现 |
| --- | --- | --- |
| 四节点轴对称 Quad4、CAX4T、CAX4RT | `rz_quad4.hpp`、`cax4t.hpp`、`cax4rt.hpp` | `rz_kernels.cpp`、`cax4t.cpp`、`cax4rt.cpp` |
| 八节点轴对称 CAX8T、CAX8RT | `rz_quad8.hpp` | `rz_quad8.cpp` |
| 三维八节点 C3D8T、C3D8RT | `cartesian3d_hex8.hpp` | `cartesian3d_kernels.cpp` |
| 三维二十节点 C3D20T、C3D20RT | `cartesian3d_hex20.hpp` | `cartesian3d_hex20.cpp` |
| 两节点轴对称边界与局部接触 | `rz_quad4.hpp`、`contact.hpp` | `rz_kernels.cpp` |
| 三节点轴对称边界与局部接触 | `line3_rz_boundary.hpp`、`line3_rz_contact.hpp` | 同名 `.cpp` |
| 三维局部接触与表面积分 | `contact.hpp`、`cartesian3d_hex8.hpp`、`cartesian3d_hex20.hpp` | `cartesian3d_contact.cpp`、`cartesian3d_hex20_contact.cpp`、相应体单元实现 |

各单元沿用具体的局部函数接口，不增加统一对象工厂。输入包括局部坐标、材料、
载荷、节点状态和已接受历史；输出为局部残量、切线矩阵及试探结果。局部自由度
顺序、积分规则、材料更新和自动微分变量数保持各单元现有约定。
边界和接触函数只处理调用方给出的局部几何与候选，整网格候选搜索、唯一归属和
全局自由度映射由 fuelsim 完成。

## CAX4T 调用与状态所有权

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

材料对象、材料函数注册表、坐标类型以及各单元的局部几何和运动学由本库持有。
已有共用类型与局部函数继续位于 `fuelsim` 命名空间；CAX4T 独立入口及新提取的
二次边界入口位于 `fuelsim::elements`。头文件迁移后直接更新调用方，不提供转发头文件。

全局网格、自由度编号、接触候选搜索与归属、问题装配、时间推进及历史状态提交
留在 fuelsim。库没有 fuelsim_core 的反向链接，也不包含父目录源码或测试支持文件。
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

另外四个独立测试为 `fuelsim_rz_element_tests`、`fuelsim_quad8_rz_tests`、
`fuelsim_hex8_tests` 和 `fuelsim_hex20_tests`，从原有测试迁移，继续覆盖局部材料、
体单元及接触残量和切线。二次轴对称边界额外检查压力、两个方向面力、热流和
对流在参考与当前构形下的中心差分导数，并检查直圆柱表面的解析合力或总热率。
全局网格和问题测试仍在仓库根目录的 `tests/` 中。

完整生产计算仍通过 `fuelsim -i <case.fsi>` 验证，保留 B7、B15 和其他现有算例的
参考数据及验收门槛。单元库的局部测试不能替代全局装配、接触和时间推进测试。

## 首次 CAX4T 迁移验证记录

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

## 其余单元迁移验证

本次迁移基线为 `dacf228`，使用相同 GCC 14.3、Release 配置和 ADlite 0.2.3。
原有体单元及局部接触源码除头文件路径外保持不变，只有二十节点曲面接触的至多
四项交点排序改用有界插入排序，解决独立非链接时优化构建中的 GCC 数组越界警告。
三节点轴对称边界从全局装配中提取，积分表达式保持原样。

独立构建的五个局部测试全部通过。八节点轴对称、三维八节点和三维二十节点测试
的标准输出与迁移前逐字一致；该比较只约束现有测试状态和打印精度，不声称对任意
输入逐位一致。拆分后的轴对称局部测试与全局核心测试合并后，56 行数值指标
与迁移前一致，瞬态核心测试的标准输出也一致。新增二次边界解析结果与中心差分
检查通过。

fuelsim Release 构建通过，统一 `ctest --test-dir build -j4 --output-on-failure`
为 **273/273 通过**，包含生产程序外部参考对比、历史事务、重启动和并行检查。
完整回归在沙盒外运行，未修改已有外部参考数据或误差门槛，也未评估性能提升。
