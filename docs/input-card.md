# fuelsim 输入卡 v1

`fuelsim` 通过一个显式输入文件选择物理问题：

```bash
./build/fuelsim -i case.fsi [PETSc options]
```

输入卡使用 MOOSE 风格的嵌套段和 `key = value`。`#` 开始行内注释；值可用
单引号或双引号包围。相对路径以输入卡所在目录为基准。v1 所有数值采用 SI
单位，不执行表达式求值或单位换算。

解析是严格的：未知段、未知键、重复段、重复键、未闭合段、非法数值和缺少
必填键都会报出文件与行号。v1 不提供 include、宏、旧键别名或兼容层。

## 物理问题

`[Case]` 的 `problem` 只能为：

- `steady_fuel_cladding`：稳态热传导与准静态热弹性，使用线性载荷步；
- `transient_fuel_cladding`：Backward Euler 热传导、准静态力学与积分点
  非弹性历史，使用可 cutback 的物理时间步。

生产 C++ 类型对应为 `SteadyFuelCladdingProblem` 和
`TransientFuelCladdingProblem`。M1/M2 仅用于路线和回归命名。

## 网格

`[Mesh]` 当前只接受 `type = exodus`。`file` 指向 `.e` 文件；`[fuel]` 和
`[cladding]` 分别指定元素块及 `radial_inner`、`radial_outer`、`bottom`、
`top` 四个命名边界。半径、高度、节点数和单元数从 Exodus 文件推导，不在
输入卡中重复声明。

I/O 层可以读取一般非结构 Quad4；当前求解装配要求选中的燃料和包壳块各自
能转换为完整的张量积 RZ 网格。

## 材料与本构

燃料和包壳的 `[Materials]` 子段都需要热弹性字段：

```text
conductivity_inverse_temperature
conductivity_constant
young_modulus
poisson_ratio
thermal_expansion
reference_temperature
```

瞬态问题还需要 `density`、`specific_heat` 和 `inelastic_model`。模型及其
条件字段为：

| `inelastic_model` | 额外字段 |
| --- | --- |
| `elastic` | 无 |
| `norton_creep` | `creep_coefficient`, `creep_reference_stress`, `creep_exponent` |
| `j2_plasticity` | `yield_stress`, `hardening_modulus` |
| `norton_creep_j2_plasticity` | 上述蠕变与塑性字段全部需要 |

为避免输入歧义，不适用于所选模型的字段也会被拒绝。

## 物理、时间推进和求解器

`[Physics]` 定义初始/外边界温度、最终体积热源、气隙导热、最小有效热隙和
接触罚参数。瞬态热源斜坡由 `heat_source_ramp_time` 定义。

稳态 `[Executioner]` 使用：

```text
type = steady
load_steps = <positive integer>
```

瞬态 `[Executioner]` 使用：

```text
type = transient
end_time
initial_time_step
minimum_time_step
maximum_time_step
growth_factor
cutback_factor
maximum_cutbacks
```

`[Solver]` 可配置 PETSc SNES 的 `absolute_tolerance`、
`relative_tolerance`、`step_tolerance` 和 `maximum_iterations`。`-i` 与输入
文件名会在 PETSc 初始化前从参数中移除，其他 PETSc 命令行选项继续生效。

`[Outputs]` 的 `console` 控制键值摘要；可选 `csv` 将相同摘要写为
`metric,value` 文件。

仓库中的可运行示例为：

- [`steady_fuel_cladding.fsi`](../verification/fuelsim/steady_fuel_cladding.fsi)
- [`transient_fuel_cladding_pcmi.fsi`](../verification/fuelsim/transient_fuel_cladding_pcmi.fsi)
