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

## 类型文件的分层

| 文件 | 职责 |
| --- | --- |
| `element_types.hpp` | 应变形式、型号选择枚举和公共计算请求，不包含固定拓扑的数组或几何 |
| `axisymmetric_types.hpp` | 轴对称应力、材料点历史、试探状态和转动类型，不依赖四节点或八节点几何 |
| `cartesian_types.hpp` | 三维应力、材料点历史、历史容器和转动类型，不包含 HEX8 专用型号字段 |
| `cax4_types.hpp`、`cax8_types.hpp` | 对应轴对称拓扑的几何、积分点、局部数组、历史容器及输入输出 |
| `c3d8_types.hpp`、`c3d20_types.hpp` | 对应三维拓扑的几何、积分点、局部数组及输入输出 |

`material.hpp` 使用轴对称和三维公共类型提供材料计算接口。四节点轴对称类型
已从 `axisymmetric_geometry.hpp` 并入 `cax4_types.hpp`，旧头文件删除；
固定为 12 个自由度的数组使用 `Cax4LocalDofs/Values/Residual/Jacobian/AdValues`
名称，不再使用容易误认为通用类型的 `Local*` 名称。

内部体单元计算直接使用对应型号族的 `Input`，不再复制一份材料和载荷数据。
主体持有的区域材料、载荷及型号选择放在 `src/core/element_region_data.hpp`，
不属于单元库公共接口。独立测试的数据准备结构位于
`elements/tests/support/element_test_data.hpp`，不进入生产库。

C3D8T 和 C3D8RT 分别通过型号头文件提供 `diagnose_c3d8t` 和 `diagnose_c3d8rt`。
调用方传入整单元几何、新旧节点状态和应变形式，获得普通双精度的
`C3d8Diagnostics`；不需要材料对象，也不执行材料更新。诊断按需调用，普通残量
和 Jacobian 请求不会计算这些输出。公共的 `c3d8_kinematics.hpp` 已删除，运动学
结构与函数声明全部位于私有的 `src/detail/c3d8_kinematics.hpp`。

诊断包含新旧构形体积和活跃材料点数据。C3D8T 返回八点，C3D8RT 返回一点，
`material_point_count` 之外的预留项为零，不参与统计或历史更新。每点返回几何
应变增量、按行排列的增量转动、热梯度重建系数、热物性求值温度及当前几何测度。
几何应变增量尚未施加型号专用的体积修正，不能当作本构应变或 Abaqus 输出应变。
当前点几何测度也不替代选择性体积处理后的力学权重；主体能量诊断继续采用原有
参考点权重乘整体体积比。

C3D8T 的热梯度使用对应点的当前构形梯度，温度使用对应角点；C3D8RT 的热梯度
和温度按体积平均。小应变使用参考构形，有限应变使用当前构形。完整运动学与
矩阵辅助计算只在库内部执行。CAX4RT 与 C3D8RT 的沙漏能量接口继续接受型号族 `Input`。

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

| 文件名称 | 职责 |
| --- | --- |
| `cax4_geometry.cpp` | CAX4 系列参考形函数、几何映射和积分测度；公共声明在 `cax4_types.hpp`，无需额外私有头文件 |
| `cax8_geometry.*`、`c3d20_geometry.*` | 对应型号族的形函数、混合阶映射和参考积分数据 |
| `c3d8_geometry.*` | C3D8 系列参考几何、节点顺序、积分与沙漏几何数据 |
| `cax8_kinematics.*`、`c3d8_kinematics.*`、`c3d20_kinematics.*` | 对应拓扑的位移梯度、构形转换及运动学；CAX8、C3D20 同时包含当前与中间构形的热几何 |
| `cax8_assembly.*`、`c3d20_assembly.*` | 型号族共用的局部热、力学装配和材料点结果计算 |
| `c3d8_diagnostics.*` | 型号诊断入口共用的整单元输出计算 |
| `cartesian_kinematics.*` | 不依赖节点数量的三维增量运动学，不计算材料切线 |
| `cartesian_material.*` | 三维材料上下文、切线、增量更新和历史张量旋转 |
| `matrix3.*` | 普通双精度与 ADlite 具体类型的三阶矩阵运算，不定义模板 |
| `ad_local_system.hpp` | 不依赖拓扑的自动微分数据转换和结果提取 |
| `line2_rz_local_system.hpp` | 两节点轴对称边界与接触共用的固定 12 自由度数组转换 |
| `line2_rz_geometry.hpp` | 两节点轴对称线段合法性检查，不包含自动微分播种或残量提取 |

`.*` 表示同名 `.hpp/.cpp`。文件按“型号族或共享范围 + 职责”命名，私有目录
保持平铺。只有跨文件使用的函数才建立私有声明；CAX4T 专用的四节点插值留在
`cax4t.cpp` 的匿名命名空间，未使用的自动微分插值重载已删除。

参考构形映射归 `geometry`；新旧构形、应变、转动、当前或中间构形的形函数梯度
及积分测度归 `kinematics`；材料函数、本构切线和历史张量更新归 `material`；
物理项积分、材料与几何导数的组合及闭式节点链装配归 `assembly`。
场温度插值和温度梯度与形函数梯度的收缩仍属于装配，不属于几何映射。
诊断只汇总几何和运动学结果，不更新材料。

CAX8、C3D20 的共用装配均只接收型号族 `Input` 和 `ElementRequest`，不再让调用方
重复传递其中已有的几何、状态、历史和时间步。请求处理集中在共用装配中，型号
入口负责积分点数量检查。私有命名空间与实际共享范围对应：`cax8_detail`、
`c3d8_detail`、`c3d20_detail`、`cartesian_detail` 和 `line2_rz_detail`。
`matrix3` 使用 `cartesian_detail`，不依赖材料或具体单元拓扑。

`tests/detail_boundaries.cmake` 检查私有文件分类、直接包含依赖及公共调用方的
私有头文件引用，并检查几何、运动学和诊断中显式混入材料调用的情况。
这属于源码边界检查，不能替代 C++ 语义审查或数值验证。
逐文件核对结果及保留不同实现的理由见
[私有实现职责核对](../docs/element-detail-boundaries.md)。

C3D20 的普通双精度残量路径和自动微分路径分别保留；文件拆分不改变播种宽度、
材料积分规则、构形检查或运算顺序。公共头文件不包含 `src/detail/`，型号专用
算法继续留在本型号 `.cpp` 内。

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

八个型号测试和一个共用轴对称边界、接触测试由普通 C++ 自检程序组成，另有一个 CMake 源码边界检查。共用断言
在 `tests/support/` 编译一次，各型号分别执行其适用的分支。物理算例继续通过
实际生产程序 `fuelsim -i <case.fsi>`，在 fuelsim 的统一 CTest 中验证。

本次整理以 `0d0454b` 为基线。原有 C3D8、C3D20 局部测试的数值指标以及 CAX8 的节点、积分点探测输出
在接口改造后与基线一致；这只约束已有测试状态和打印精度，不代表任意输入逐位
一致。旧默认 Quad4 的专用检查已删除，不能把它们的旧通过记录用于 CAX4T。
最终构建、回归结果及输入迁移范围见 [整理与验证记录](../docs/element-library-refactor.md)。
