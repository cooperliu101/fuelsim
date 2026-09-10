# fuelsim 单元计算库

本目录与 fuelsim 位于同一仓库，可仅依赖 C++17 和 ADlite 独立构建、测试。
公共头文件直接放在 `include/`；每种体单元在 `src/` 和 `tests/` 中有同名文件。

## 型号与接口

| 型号 | 公共头文件 | 实现 | 测试 |
| --- | --- | --- | --- |
| CAX4T | `cax4t.hpp` | `cax4t.cpp` | `cax4t_tests.cpp` |
| CAX4RT | `cax4rt.hpp` | `cax4rt.cpp` | `cax4rt_tests.cpp` |
| CAX8T | `cax8t.hpp` | `cax8t.cpp` | `cax8t_tests.cpp` |
| CAX8RT | `cax8rt.hpp` | `cax8rt.cpp` | `cax8rt_tests.cpp` |
| C3D8T | `c3d8t.hpp` | `c3d8t.cpp` | `c3d8t_tests.cpp` |
| C3D8RT | `c3d8rt.hpp` | `c3d8rt.cpp` | `c3d8rt_tests.cpp` |
| C3D20T | `c3d20t.hpp` | `c3d20t.cpp` | `c3d20t_tests.cpp` |
| C3D20RT | `c3d20rt.hpp` | `c3d20rt.cpp` | `c3d20rt_tests.cpp` |

入口统一为 `fuelsim::elements::evaluate_<型号>(input, request)`。同一拓扑的
输入与输出具有相同数组尺寸，因此共用具体类型 `Cax4Input/Result`、
`Cax8Input/Result`、`C3d8Input/Result` 和 `C3d20Input/Result`，分别声明在
`cax4_types.hpp`、`cax8_types.hpp`、`c3d8_types.hpp` 和 `c3d20_types.hpp`。
没有模板、对象工厂、旧接口别名或型号之间的互相调用。

型号文件负责本型号的积分规则和专用算法。C3D8T 的选择性体积处理与 C3D8RT
的均匀应变、沙漏计算分别位于各自实现中。CAX8T/CAX8RT 和 C3D20T/C3D20RT
通过各自的 `make_<型号>_geometry` 选择积分阶数，共用形函数和积分点装配，
并在计算入口检查几何所含的材料积分点数量。

默认轴对称 `quad4` 单元已经删除，包括其 Taylor/Rashid 增量算法、体积分派和
专用测试。每个生产输入区域必须显式设置 `element`；不接受 `quad4`，也不将它
作为 CAX4T 的别名。四节点形函数、坐标映射及四节点三维表面仍属于共用几何。

## 调用边界与状态所有权

接口参考 UEL 的局部计算职责边界，不实现 Abaqus 的二进制调用约定。

- 输入借用材料、参考几何、当前节点状态、已接受节点状态和材料历史；库不保存
  这些引用，也不修改输入。历史为空表示稳态热弹性；非空历史表示瞬态材料更新。
- 局部自由度按场排列，分别为 12、20、32、68 个；具体顺序遵循各型号约定。
- `ElementRequest` 控制残量、Jacobian（残量对局部状态的导数矩阵）、试探历史
  及应力输出。Jacobian 按行存储，无额外负号。默认请求残量和历史。
  `stress` 请求控制独立应力输出，轴对称应力从同次局部计算的试探历史提取。
- 轴对称局部实现联合计算残量与历史，因此关闭相应请求只清空不需要的返回字段。
  三维实现可分别执行残量、历史和应力计算。请求选项不改变材料离散方法。
- `include_thermal_time_term` 只控制热容项，不关闭塑性或蠕变积分。材料更新使用
  `time_step`，温度相关材料函数使用该型号规定的节点或积分点温度。
- 返回历史只代表试探结果，接受或丢弃由 fuelsim 决定。非法试探状态通过异常
  交给调用方处理，不隐式夹持温度、材料值或几何。
- 轴对称结果中的 `stored_heat_rate` 和 `generated_heat_rate` 为功率，使用
  对应型号的体积分测度；全局边界、接触及时间积分后的能量汇总由 fuelsim 负责。

## 公共代码与内部代码

`include/` 中的共用类型用于调用接口：坐标、材料、材料函数、几何、边界参数和
接触历史。`src/detail/` 只包含内部实现：

- `axisymmetric_geometry`：四节点轴对称形函数、几何映射和参考积分测度。
- `hex8_geometry`：八节点六面体共用几何和运动学连接。
- `quad8_rz_assembly`、`hex20_assembly`：同一拓扑的共用积分点计算和装配。
- `cartesian_kinematics`：三维矩阵和增量运动学。
- `material_rotation`：三维材料增量更新和客观历史旋转。
- `ad_local_system`、`rz_local_system`、`rz_point`：局部自动微分及节点链辅助计算。

只供一个型号使用的算法留在其 `.cpp` 内。公共头文件不包含 `src/detail/`。

边界与接触按线段或表面拓扑命名，并分别拥有对应的 `.hpp/.cpp`：

- `line2_rz_boundary`、`line3_rz_boundary`；
- `quad4_face_boundary`、`quad8_face_boundary`；
- `line2_rz_contact`、`line3_rz_contact`；
- `quad4_face_contact`、`quad8_face_contact`。

这些函数处理调用方提供的局部几何与候选。全局网格、自由度编号、候选搜索、
唯一归属、全局装配、时间推进、历史提交及失败恢复由 fuelsim 管理。
生产侧的型号选择位于 `src/core/cax4_evaluation.hpp` 和
`src/core/element_evaluation.hpp`；独立测试使用 `tests/support/` 中的测试选择函数。

## 构建与验证

在仓库根目录运行：

```bash
cmake -S elements -B build-elements \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH=/home/cooper/ai_project/fuelsim-dependencies/adlite-0.2.3
cmake --build build-elements --parallel 4
ctest --test-dir build-elements -j4 --output-on-failure
```

八个型号测试和一个共用轴对称边界、接触测试由普通 C++ 自检程序组成。共用断言
在 `tests/support/` 编译一次，各型号分别执行其适用的分支。物理算例继续通过
实际生产程序 `fuelsim -i <case.fsi>`，在 fuelsim 的统一 CTest 中验证。

本次整理以 `0d0454b` 为基线。原有 C3D8、C3D20 局部测试的数值指标以及 CAX8 的节点、积分点探测输出
在接口改造后与基线一致；这只约束已有测试状态和打印精度，不代表任意输入逐位
一致。旧默认 Quad4 的专用检查已删除，不能把它们的旧通过记录用于 CAX4T。
最终构建、回归结果及输入迁移范围见 [整理与验证记录](../docs/element-library-refactor.md)。
