# fuelsim 单元计算库

本目录与 fuelsim 位于同一仓库，使用 C++17 和 ADlite，可独立配置、构建和测试。
`include/` 和 `src/` 均保持平铺，不建立 `detail` 目录。

## 按类型组织

八种体单元分别拥有同名公共头文件、实现文件和测试文件：
`cax4t`、`cax4rt`、`cax8t`、`cax8rt`、`c3d8t`、`c3d8rt`、`c3d20t`、`c3d20rt`。
型号专用计算与辅助函数留在本型号实现中。CAX8T 和 C3D20T 分别承载对应二次
单元的完整算法；CAX8RT 和 C3D20RT 包含对应全积分型号头文件，通过显式积分
规则调用其算法、几何构造及适用的构形检查。其他型号之间不互相调用。
全积分型号原有入口仍固定使用全积分规则；带积分规则的重载严格检查所选规则
与材料积分点数量一致。C3D20 的热积分与力学积分仍按各自原有规则执行。

四种界面类型为 `line2_rz`、`line3_rz`、`quad4_face`、`quad8_face`，同样各有
一组 `.hpp/.cpp` 和测试。每组文件同时维护本类型的边界载荷、热接触和机械接触，
以不同函数区分计算任务，不再按 `boundary/contact` 拆文件。

| 文件 | 归属 |
| --- | --- |
| `src/cax_common.hpp/.cpp` | CAX4T/RT 实际共用的参考几何，以及三个轴对称实现共用的 Hughes-Winget 转动；专用体积处理和沙漏控制仍在型号实现中 |
| `src/c3d_common.hpp/.cpp` | 三维型号共用的矩阵与运动学、Hughes-Winget 转动，以及 HEX8 角点与沙漏常量、减缩热沙漏系数、参考几何和历史几何计算 |
| `include/material.hpp`、`src/material.cpp` | 轴对称和三维材料响应、本构切线、材料历史及其客观旋转的声明与实现 |
| `src/material_functions.cpp` | 材料函数注册、具名参数和函数组合 |
| `src/ad_local_system.hpp` | 多种局部计算共用的自动微分数组初始化与结果提取 |

`src/` 下的头文件都是私有声明，公共调用方不应包含。只在一个实现文件使用的
辅助函数留在该文件中，不为单独的几何、运动学、装配或诊断阶段再建立文件。

型号实现及其共用实现可以依赖材料接口；材料实现、材料类型和材料函数不得反向
包含体单元型号、界面、型号族类型或 `cax_common`、`c3d_common` 头文件。
材料接口只使用基础坐标、材料类型与 ADlite；张量旋转所需的矩阵计算留在材料实现内。

## 公共类型

`element_types.hpp` 保存坐标、应变形式、型号枚举和请求；`material_types.hpp`
保存轴对称和三维材料点应力、历史及转动类型。四个型号族的 `*_types.hpp`
只保存对应拓扑、局部数组和 T/RT 共用的输入输出，不声明计算函数。专用界面数据放在各界面头文件，
共用接触参数和历史放在 `contact_types.hpp`。`fnv_hash.hpp` 同时服务材料函数
签名和主体文件签名，继续作为共用工具保留。

几何构造统一由型号头文件声明 `make_<型号>_geometry`；CAX8 和 C3D20 的
几何构造由对应全积分型号实现，其余型号按实际需要复用私有共用函数。
三维型号头文件还声明 `validate_<型号>_deformation(point, state)`，供主体独立检查
给定参考积分点的当前变形；它不替代单元计算内部对新旧及中间构形的完整检查。
CAX8 的积分点评估位于 `cax8t.cpp`，C3D20 的积分点评估位于 `c3d20t.cpp`；
这些仅供本文件使用的辅助函数保留在匿名命名空间内。

## 局部调用边界

体单元入口为 `fuelsim::elements::evaluate_<型号>(input, request)`，参考 UEL
的局部计算职责，不实现 Abaqus 二进制调用约定。库借用参考几何、材料、新旧
节点状态和已接受材料历史，返回局部残量、Jacobian（残量对局部状态的导数矩阵）、
试探历史及按需应力。全局网格、自由度编号、装配、时间推进、历史提交及失败恢复
均由 fuelsim 主体负责。

局部自由度保持按场排列，分别为 12、20、32、68 个。轴对称自动微分采用宽度 6
的运动学与宽度 5 的材料链，三维分别采用宽度 10 和 7；闭式节点链装配保持原有
规则。不会按整个体单元自由度播种。非法材料和构形状态继续通过异常拒绝。

C3D8 独立诊断接口和诊断结构已经删除。正常材料历史更新返回 `current_volume`、
`committed_volume` 和 `incremental_rotations`，供主体能量统计使用；这些字段
只在请求历史且提供已接受历史时计算。活跃点数量由返回历史长度确定，未使用的
转动项保持零。纯残量或 Jacobian 请求不会额外计算这些历史更新结果。

用于验证的热梯度、型号温度与几何应变重建位于测试支持代码，以普通双精度独立
计算，不依赖单元库私有实现。生产代码不包含测试支持头文件。

## 测试与构建

```bash
cmake -S elements -B build-elements \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH=/home/cooper/ai_project/fuelsim-dependencies/adlite-0.2.3
cmake --build build-elements --parallel 4
ctest --test-dir build-elements -j4 --output-on-failure
```

独立测试包含八种体单元、四种界面及一项源码边界检查。专用检查位于对应类型的
测试文件；重复使用的检查保留在型号族测试支持中，材料准备和测试调用辅助集中
在 `tests/support/test_support.hpp`。原 CAX8 两个型号重复执行的同一界面检查
集中到 `line3_rz_tests.cpp` 执行一次。

`tests/library_boundaries.cmake` 检查平铺结构、类型对应文件、型号之间的调用和
私有头文件访问边界，以及材料文件对型号和界面头文件的直接依赖。它是源码检查，
不代替编译和数值测试。生产算例仍通过
`fuelsim -i <case.fsi>` 运行，并由主体统一 CTest 验证。

## 局部数学计算的模板例外

轴对称三个实现通过 `cax_common` 的普通 `adlite::Scalar` 函数计算面内增量转动
和随转应变。三维五个调用位置通过 `c3d_common.hpp` 中的私有函数模板
`hughes_winget_rotation` 共用相同算法；模板仅接受 `double` 和 `adlite::Scalar`，
普通双精度路径不构造自动微分标量。标量数学函数采用 `using std::函数名` 加
非限定调用，使自动微分参数通过参数相关查找选中 ADlite 重载。两者均接收中间构形位移梯度，不进行播种。
构形合法性、环向应变、体积平均、沙漏控制及闭式节点导数链仍属于各型号。
此例外不允许引入通用矩阵模板、材料模板或单元标量泛型接口。

HEX8 角点符号和原始沙漏模态统一保存于 `c3d_common.hpp`，下标分别为节点—方向
和节点—模态。C3D20 的角点形函数使用同一符号顺序。减缩热沙漏系数由
`c3d_common.cpp` 的两个具体类型重载计算，输入为相应构形的平均梯度和体积；
两个重载调用实现文件中的同一个私有函数模板，公式只保留一份。参考、当前和
中间构形调用双精度接口，自动微分构形调用另一接口；模板同样仅允许 `double`
和 `adlite::Scalar`，数学函数遵循上述参数相关查找规则。
C3D8RT 的闭式几何导数仍在型号文件中计算，参考几何错误与试探构形错误保留
各自的异常类型和构形说明。
