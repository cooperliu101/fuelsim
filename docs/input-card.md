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

## 自由区域组合

`[Regions]` 下每个子段定义一个物理区域，子段名是区域名，`block` 是同一
Exodus 文件中的元素块名：

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
  heat_source_ramp_time = 20
[]
```

`[Solver]` 可设置 `absolute_tolerance`、`relative_tolerance`、
`step_tolerance` 和 `maximum_iterations`；省略时分别为 `1e-8`、`1e-10`、
`1e-12` 和 `40`。其他 PETSc 命令行选项仍可直接覆盖默认行为。

`[Outputs]` 的 `console` 默认为 `true`；可选 `csv` 将同一组命名指标写为
`metric,value` 文件。

仓库中的可运行示例为：

- [`steady_fuel_cladding.fsi`](../verification/fuelsim/steady_fuel_cladding.fsi)
- [`transient_fuel_cladding_pcmi.fsi`](../verification/fuelsim/transient_fuel_cladding_pcmi.fsi)
