# fuelsim 输入卡 v1

`fuelsim` 只使用一个显式输入文件：

```bash
./build/fuelsim -i case.fsi [PETSc options]
```

输入卡采用 MOOSE 风格的嵌套段和 `key = value`。`#` 开始行内注释；相对
路径以输入卡所在目录为基准。v1 的数值全部使用 SI 单位，不支持 include、
宏、表达式、单位换算、旧键别名或兼容层。未知段、未知键、重复项、非法值
和缺少必填键都会立即报错。

## 问题与单一网格

`[Case]` 的 `problem` 只能是：

- `steady`：稳态热传导和准静态热弹性，使用线性载荷步；
- `transient`：Backward Euler 热传导、准静态力学和积分点非弹性历史。

对应的两个生产 C++ 类型是 `SteadyProblem` 和 `TransientProblem`。M0、M1、
M2 只作为路线和回归名称。

`[Mesh]` 只接受一个 Exodus 文件：

```text
[Mesh]
  type = exodus
  file = model.e
[]
```

文件可以包含多个 Quad4 元素块、节点集和边集。I/O 层保留这些元数据，
求解区域保留选中元素块的原始节点坐标和连接关系；区域不需要是张量积网格。
每个 Quad4 在积分点必须具有正 Jacobian。输入卡不重复定义半径、高度和
离散规模。

Dirichlet 边界可使用任意属于该区域的边集。当前热接触和机械接触仍要求
`primary`、`secondary` 是轴对称圆柱面，并由边集相邻单元自动确定内外侧；
这一接触几何限制独立于体网格是否结构化。

## 时间函数

瞬态卡可定义具体的分段线性时间表：

```text
[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 2.5 5 10
    values = 0 1 0.5 1
  []
[]
```

`times` 和 `values` 是长度相同、至少包含两个值的空格分隔列表。时间必须
有限、非负且严格递增，值必须有限；首末区间外保持端点值，不外推。表内每个
时刻都是执行器必须准确命中的事件。名称必须唯一，引用未知名称会在解析时
报错。稳态卡不接受 `[TimeFunctions]`。

## 自由区域组合

`[Regions]` 下每个子段定义一个物理区域，子段名是区域名。每个区域必须
用且只能用 `block` 或 `block_id` 选择同一 Exodus 文件中的元素块；有块名
时优先使用更易读的 `block`，无块名的 MOOSE 单区域网格可使用非负
`block_id`：

```text
[Regions]
  [pellet]
    block = fuel
    conductivity_inverse_temperature = 3824
    conductivity_constant = 0.61
    young_modulus = 2e11
    poisson_ratio = 0.316
    thermal_expansion = 1e-5
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 2e8
    heat_source_function = power
  []
[]
```

区域数量不固定，因此同一结构可表示单独芯块、单独包壳、芯块—包壳，或
芯块—包壳1—包壳2。每个元素块只能声明一次，区域节点与自由度保持独立。

瞬态问题的每个区域还必须给出 `density`、`specific_heat` 和
`inelastic_model`。可选模型及条件字段为：

| `inelastic_model` | 额外字段 |
| --- | --- |
| `elastic` | 无 |
| `norton_creep` | `creep_coefficient`, `creep_reference_stress`, `creep_exponent` |
| `j2_plasticity` | `yield_stress`, `hardening_modulus` |
| `norton_creep_j2_plasticity` | 上述蠕变与塑性字段全部需要 |

不适用于所选模型的字段会被拒绝。

`heat_source_function` 可选；存在时，当前体积热源为
`volumetric_heat_source * function(time)`。未设置时沿用执行器
`load_ramp_time` 的全局载荷因子。

## 接触

每个 `[Contact/<name>]` 只声明 MOOSE 风格的 `primary` 和 `secondary` 边集。
不得声明 `primary_block` 或 `secondary_block`；fuelsim 通过边集相邻单元读取
块 ID，再匹配到 `[Regions]`：

```text
[Contact]
  [pellet_clad]
    primary = clad_inner
    secondary = pellet_outer

    [thermal]
      gap_conductivity = 0.4
      minimum_gap = 1e-6
    []

    [mechanical]
      formulation = penalty
      penalty = 1e14
    []
  []
[]
```

一个接触对至少包含 `[thermal]` 或 `[mechanical]`，也可以同时包含两者。
当前 RZ 实现要求 primary 是外侧圆柱面、secondary 是内侧圆柱面；两者必须
来自不同区域且各自形成连续的轴向边链。热接触采用 secondary-side STS
积分，机械接触采用 secondary 节点到 primary 线段的唯一 NTS 投影。

`[Contact]` 段本身是必需的，但可以为空，以支持不含接触的单区域问题。

## 边界条件

边界条件同样只引用 Exodus 边集，所属区域由相邻单元自动确定：

```text
[BoundaryConditions]
  [axis]
    type = dirichlet
    boundary = fuel_axis
    field = radial_displacement
    value = 0
  []

  [external_pressure]
    type = pressure
    boundary = clad_outer
    value = 1.5e7
  []
[]
```

`dirichlet` 的 `field` 只能为 `temperature`、`radial_displacement` 或
`axial_displacement`。`pressure` 不接受 `field`，当前只能施加在径向边界。
`traction` 必须声明一个位移 `field`，`value` 是该分量上的有符号表面牵引；
它使用参考 RZ 表面测度积分。`dirichlet`、`pressure` 和 `traction` 都可设置
`scale_with_load = true`，使 `value` 乘以当前执行器载荷因子；默认不缩放。
也可用 `function = <name>` 使 `value` 乘以时间表值；`function` 与
`scale_with_load` 互斥。压力在任一求值时刻都必须非负。

对流热边界写为：

```text
[BoundaryConditions]
  [coolant]
    type = convection
    boundary = clad_outer
    heat_transfer_coefficient = 1000
    ambient_temperature = 500
    coefficient_function = coolant_flow
    ambient_temperature_function = coolant_temperature
  []
[]
```

两个函数均可省略；存在时分别乘以对应基值。换热系数允许时间表计算为零但
不能为负，环境温度必须始终为正。对流项使用参考 RZ 表面测度并作为邻接
Quad4 的 12-DOF ADlite 局部贡献装配，因此残量和温度切线保持一致。
`[BoundaryConditions]` 段本身必需，但可以为空。

## 时间推进、求解与输出

稳态执行器为：

```text
[Executioner]
  type = steady
  load_steps = 20
[]
```

瞬态执行器为：

```text
[Executioner]
  type = transient
  end_time = 20
  initial_time_step = 1
  minimum_time_step = 0.125
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 20
  restart = previous.checkpoint
[]
```

瞬态载荷因子为 `min(time/load_ramp_time, 1)`；`load_ramp_time = 0` 表示从
首步起使用完整载荷。该因子同时控制各区域 `volumetric_heat_source` 和所有
显式设置 `scale_with_load = true` 的边界条件。

时间推进会在不改变名义时间步控制器的前提下截短当前步，以准确命中下一
时间表节点；随后恢复原名义步长增长路径。事件间隔可以小于
`minimum_time_step`，因为事件时刻优先于最小重试步长；若该事件步求解失败，
后续 cutback 仍受最小步长约束。

`restart` 为可选的严格重启动文件。它恢复已提交的节点温度/位移、物理时间、
载荷因子、全部积分点塑性/蠕变历史和已提交应力。文件的版本、字节序、长度、
校验和、网格、材料、边界条件、接触及局部装配拓扑必须与当前问题一致；不
匹配时立即停止。重启动不保存 Newton trial、活动时间步或失败尝试。

`[Solver]` 可设置 `absolute_tolerance`、`relative_tolerance`、
`step_tolerance` 和 `maximum_iterations`；省略时分别为 `1e-8`、`1e-10`、
`1e-12` 和 `40`。其他 PETSc 命令行选项仍可直接覆盖默认行为。

`[Outputs]` 的 `console` 默认为 `true`；可选 `csv` 将同一组命名指标写为
`metric,value` 文件。`exodus` 写出可后处理的场结果；稳态写一个最终步，
瞬态写初始/重启动状态和每个成功提交的时间步。节点变量包括温度、径向与
轴向位移，以及各接触对 secondary 节点上的间隙和压力；单元变量保留四个
积分点的应力、塑性应变、蠕变应变及两种等效应变；全局变量记录载荷因子、
界面总热流和总接触力。未属于所选求解区域或未投影的值写为 `NaN`。

瞬态还可设置 `checkpoint` 和正整数 `checkpoint_interval`：

```text
[Outputs]
  console = true
  csv = summary.csv
  exodus = fields.e
  checkpoint = latest.checkpoint
  checkpoint_interval = 5
[]
```

检查点只在成功提交后按间隔原子替换，并在执行结束或失败退出前再次保存最后
提交态。`checkpoint_interval` 默认是 `1`，不能单独出现；稳态不接受检查点。
结果文件不得覆盖输入网格，也不得与检查点同名。相对路径都以输入卡目录为
基准。重启动后的 Exodus 输出会新建一个结果文件，并以恢复时刻作为第一个
结果步，不尝试修改上一段结果文件。

仓库中的可运行示例为：

- [`steady_single_fuel_moose.fsi`](../verification/fuelsim/steady_single_fuel_moose.fsi)
- [`steady_fuel_cladding.fsi`](../verification/fuelsim/steady_fuel_cladding.fsi)
- [`steady_fuel_cladding_unstructured.fsi`](../verification/fuelsim/steady_fuel_cladding_unstructured.fsi)
- [`transient_heat_moose.fsi`](../verification/fuelsim/transient_heat_moose.fsi)
- [`transient_table_convection_moose.fsi`](../verification/fuelsim/transient_table_convection_moose.fsi)
- [`transient_j2_plastic_moose.fsi`](../verification/fuelsim/transient_j2_plastic_moose.fsi)
- [`transient_norton_creep_moose.fsi`](../verification/fuelsim/transient_norton_creep_moose.fsi)
- [`transient_coupled_displacement_moose.fsi`](../verification/fuelsim/transient_coupled_displacement_moose.fsi)
- [`transient_coupled_traction_moose.fsi`](../verification/fuelsim/transient_coupled_traction_moose.fsi)
- [`transient_fuel_cladding_pcmi.fsi`](../verification/fuelsim/transient_fuel_cladding_pcmi.fsi)

上述十张卡分别驱动 M0、两套 M1、M2.1、M3.1、四套 M2.2 和 M2.3 的
fuelsim-to-MOOSE 对比；测试程序不再直接构造这些案例的材料、载荷路径或
网格选择参数。每个对比读取 MOOSE 最终时刻的全部节点，统一检查温度、
径向位移和轴向位移的三项误差；M1 和 M2.3 还检查全部接触节点的压力三项
误差。
