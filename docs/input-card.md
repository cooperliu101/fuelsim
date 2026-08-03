# fuelsim 输入卡 v1

`fuelsim` 只使用一个显式输入文件：

```bash
./build/fuelsim -i case.fsi [PETSc options]
```

输入卡采用 MOOSE 风格的嵌套段和 `key = value`。`#` 开始行内注释；相对
路径以输入卡所在目录为基准。v1 的数值全部使用 SI 单位，不支持 include、
宏、表达式、单位换算、旧键别名或兼容层。未知段、未知键、重复项、非法值
和缺少必填键都会立即报错。

每张卡都必须显式声明格式版本和物理问题：

```text
[Case]
  version = 1
  problem = transient
[]
```

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

Dirichlet 和接触边界可使用任意属于所选区域的边集。每个接触面必须是一条
不分叉的开放 Line2 边链；圆柱侧面、水平芯块端面和斜面均可使用。接触两侧
的分段无须匹配，但 secondary 的投影必须被 primary 边链完整且唯一地覆盖。

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
    strain = small
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

每个区域必须显式设置 `strain = small` 或 `strain = finite`。`small` 使用
参考构形小应变弱式；`finite` 使用 MOOSE 默认的增量 Taylor 应变与 Rashid
转动，并用 Cauchy 应力、当前构形梯度和当前 RZ 测度装配力学内力。热传导
与热容仍使用参考构形。有限应变区域的 pressure 是当前构形 follower load；
traction 仍是参考构形 dead load。区域发生非正 Jacobian、非正环向伸长或
非正当前半径时会拒绝 Newton 试探态，不做隐式夹持。

瞬态问题的每个区域还必须给出 `density`、`specific_heat` 和
`inelastic_model`。可选模型及条件字段为：

| `inelastic_model` | 额外字段 |
| --- | --- |
| `elastic` | 无 |
| `norton_creep` | `creep_coefficient`, `creep_reference_stress`, `creep_exponent` |
| `j2_plasticity` | `yield_stress`, `hardening_modulus` |
| `norton_creep_j2_plasticity` | 上述蠕变与塑性字段全部需要 |

不适用于所选模型的字段会被拒绝。

热弹性参数还可设置
`young_modulus_temperature_coefficient`、
`poisson_ratio_temperature_coefficient` 和
`thermal_expansion_temperature_coefficient`；所选非弹性模型可对应设置
`creep_coefficient_temperature_coefficient`、
`creep_reference_stress_temperature_coefficient`、
`creep_exponent_temperature_coefficient`、
`yield_stress_temperature_coefficient` 和
`hardening_temperature_coefficient`。所有斜率默认为零，定义为
`property(T)=property_ref+slope*(T-reference_temperature)`，单位为对应属性
每 K。参数以 `adlite::Scalar` 活跃求值，运行中越过物理定义域会拒绝当前
Newton 试探值而不会夹持。它们是线性算法接口，不是已标定的真实材料模型。

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
两者必须来自不同区域且各自形成一条不分叉的开放边链。热接触按投影重叠
区间切分 secondary-side STS 积分，机械接触采用 secondary 节点到 primary
线段的唯一 NTS 投影；二维法向同时装配径向和轴向反力。当前是小滑移候选面
策略：每个从节点预建参考最近主段及相邻段，位移不得跨越更多主段。
Newton 试探状态一旦离开候选窗口会作为物理域错误交给回溯线搜索；若仍无法
恢复则拒绝当前载荷步或时间步，不再静默返回零接触力。主面链首尾保留端点
支承；内部顶点使用半开区间，任何时刻只允许一条主段拥有该从节点。热接触
采用相同的候选窗口越界策略。

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
`axial_displacement`。`pressure` 不接受 `field`，可施加在任意不退化的
Line2 外边界；方向取边界相邻 Quad4 的外法向。小应变区域使用参考 RZ 表面，
有限应变区域使用当前半径、当前法向和当前表面测度。`traction` 必须声明
一个位移 `field`，`value` 是该分量上的有符号参考构形表面牵引。
`dirichlet`、`pressure` 和 `traction` 都可设置
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
  cutback_factor = 0.5
  maximum_cutbacks = 12
  minimum_load_increment = 1e-6
[]
```

后三项可省略并使用上示默认值。名义载荷步失败时，执行器缩小从最近成功载荷
到目标载荷的增量；成功的中间状态成为下一次尝试的初值。最小载荷增量会实际
尝试一次后才报告失败。

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
  target_nonlinear_iterations = 6
  iteration_window = 2
  time_error_relative_tolerance = 2e-4
  temperature_time_absolute_tolerance = 1e-3
  displacement_time_absolute_tolerance = 1e-10
  time_error_safety_factor = 0.9
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

`target_nonlinear_iterations` 和 `iteration_window` 可选。目标为正时，成功步的
非线性迭代数低于 `target-window` 会按 `growth_factor` 增长下一名义步长，
高于 `target+window` 会按 `cutback_factor` 缩短下一步，窗口内保持不变；
结果始终限制在最小/最大步长内。目标省略或为零时，保持每个成功步均增长的
原行为，此时窗口必须为零。窗口必须小于目标。

`time_error_relative_tolerance` 省略或为零时不做时间离散误差控制。设为正值
后，每个候选步从同一 committed 状态计算一个 Backward Euler 全步和两个
半步；逐场归一化 L2 差的最大值大于 1 时完整回滚并缩步，成功时采用两个
半步的结果。温度和两个位移场分别使用上述绝对容差，安全系数必须位于
`(0,1)`。该估计器增加到约三倍的非线性求解工作量，但直接控制时间截断误差，
并在进度和拒步诊断中输出 `time_error_estimate`。

每次未收敛尝试都会记录尝试终点、步长、cutback 序号、非线性迭代数、PETSc
收敛原因、残量范数、失败类别和物理域消息；最终停止原因区分 `completed`、
`maximum_cutbacks` 和 `minimum_time_step`。拒绝步仍完整回滚 committed
状态；最小时间步会实际求解一次，只有该次也失败才终止。

`restart` 为可选的严格重启动文件。它恢复已提交的节点温度/位移、物理时间、
载荷因子、全部积分点塑性/蠕变历史和已提交应力。文件的版本、字节序、长度、
校验和、网格、材料、边界条件、接触及局部装配拓扑必须与当前问题一致；不
匹配时立即停止。重启动不保存 Newton trial、活动时间步或失败尝试。
格式 v2 还保存成功提交后控制器给出的下一名义时间步，因此自适应计算从检查
点继续时不会重新使用输入卡的初始步长。

`[Solver]` 可设置非线性 `absolute_tolerance`、`relative_tolerance`、
`step_tolerance` 和 `maximum_iterations`；省略时分别为 `1e-8`、`1e-10`、
`1e-12` 和 `40`。线性选项为：

- `linear_solver = automatic|direct|gmres`；
- `preconditioner = automatic|lu|block_jacobi|field_split|hypre`；
- `linear_relative_tolerance`，默认 `1e-8`；
- `maximum_linear_iterations`，默认 `500`。
- `backtracking_fallback`，默认 `true`；BASIC 失败后从原始初值用 BT 重试；
- `residual_reduction_tolerance`，默认 `1e-6`，用于总残量和分场残量复核；
- `temperature_residual_absolute_tolerance`，默认 `1e-8 W`；
- `mechanical_residual_absolute_tolerance`，默认 `1e-4 N`，同时用于径向和轴向；
- `field_residual_scaling`，默认 `false`，可选启用热/力分组的自动行缩放。

`automatic` 使用直接 LU：单 rank 采用 PETSc LU，多 rank 采用 PETSc 的 MUMPS
分解。选择 `block_jacobi`、`field_split` 或 `hypre` 会自动选 GMRES；
`field_split` 按固定 `[T(:)]` 和 `[ur(:), uz(:)]` 建立乘法场分裂。具体 PETSc
命令行选项仍在上述设置之后生效，可用于选择 HYPRE 子类型和场分裂子 KSP。
非线性默认使用 PETSc BASIC 全步。BASIC 失败时，默认从本次求解的原始初值
自动用 BT 回溯重试；BT 仍失败才由稳态载荷二分或瞬态 cutback 恢复。
PETSc `-snes_linesearch_type` 仍可覆盖具体类型；输入卡可关闭自动回退。
即使 PETSc 返回正收敛原因，fuelsim 仍按配置的绝对/相对门槛复核最终残量，
步长停滞不能单独算作成功。因热—力残量单位不同，总残量和温度/径向/轴向
三个分场均须达到绝对门槛或至少按 `residual_reduction_tolerance` 相对降低，
并设数值噪声底线。自动行缩放是实验性可选项：残量和 Jacobian 同行缩放，
但其初始残量尺度随载荷步改变，默认验证路径只启用分场诊断和复核。
两个分场绝对门槛具有明确物理单位；极小载荷或不同量级模型应在输入卡中按
所需平衡精度显式收紧或放宽，不能用混合单位的总残量容差替代。

`[Outputs]` 的 `console` 默认为 `true`；可选 `csv` 将同一组命名指标写为
`metric,value` 文件。`exodus` 写出可后处理的场结果；稳态写一个最终步，
瞬态按 `exodus_interval` 写初始/重启动状态、成功提交步及最终状态。节点变量包括温度、径向与
轴向位移，以及各接触对 secondary 节点上的间隙和压力；单元变量保留四个
积分点的应力、塑性应变、蠕变应变及两种等效应变；全局变量记录载荷因子、
界面总热流和总接触力。未属于所选求解区域或未投影的值写为 `NaN`。

瞬态还可设置输出频率、工程标量时程、进度和检查点：

```text
[Outputs]
  console = true
  csv = summary.csv
  exodus = fields.e
  exodus_interval = 10
  history = engineering_history.csv
  history_interval = 1
  progress_interval = 10
  checkpoint = latest.checkpoint
  checkpoint_interval = 5
[]
```

`exodus_interval`、`history_interval` 和 `progress_interval` 均按成功步计数，
默认为 `1`；无论频率如何，结束或失败时仍写出最后 committed 状态。工程时程
按区域名称记录最高温度及最大塑性/蠕变等效应变，按接触名称记录最小间隙、
最大压力、总热流和总反力。检查点只在成功提交后按间隔原子替换，并在执行
结束或失败退出前再次保存最后提交态。

CSV、Exodus、工程时程和检查点彼此不得同名，也不得覆盖输入卡、输入网格或
重启动检查点。相对路径都以输入卡目录为基准。重启动后的 Exodus 和工程
时程自动使用首个空闲的 `.partN` 文件，并以恢复时刻作为第一条记录，不修改
上一段结果。

全局场级 Jacobian 诊断入口为：

```bash
fuelsim -i case.fsi --check-jacobian
```

该模式不执行完整载荷路径；稳态在完整载荷初值、瞬态在下一物理时间步的
committed 初值上装配解析方向导数，并与中心差分比较。输出按温度、径向位移
和轴向位移分别给出残量 L2、解析/差分方向导数 L2、差值 L2、相对 L2 和最大
绝对差。Dirichlet 行采用与 PETSc 回调完全相同的 `F_i=x_i-g_i`。

仓库中的可运行示例为：

- [`steady_single_fuel_moose.fsi`](../verification/fuelsim/steady_single_fuel_moose.fsi)
- [`steady_fuel_cladding.fsi`](../verification/fuelsim/steady_fuel_cladding.fsi)
- [`steady_fuel_cladding_unstructured.fsi`](../verification/fuelsim/steady_fuel_cladding_unstructured.fsi)
- [`steady_two_pellet_contact_moose.fsi`](../verification/fuelsim/steady_two_pellet_contact_moose.fsi)
- [`transient_heat_moose.fsi`](../verification/fuelsim/transient_heat_moose.fsi)
- [`transient_table_convection_moose.fsi`](../verification/fuelsim/transient_table_convection_moose.fsi)
- [`transient_j2_plastic_moose.fsi`](../verification/fuelsim/transient_j2_plastic_moose.fsi)
- [`transient_norton_creep_moose.fsi`](../verification/fuelsim/transient_norton_creep_moose.fsi)
- [`transient_coupled_displacement_moose.fsi`](../verification/fuelsim/transient_coupled_displacement_moose.fsi)
- [`transient_coupled_traction_moose.fsi`](../verification/fuelsim/transient_coupled_traction_moose.fsi)
- [`transient_fuel_cladding_pcmi.fsi`](../verification/fuelsim/transient_fuel_cladding_pcmi.fsi)
- [`transient_finite_strain_pcmi.fsi`](../verification/fuelsim/transient_finite_strain_pcmi.fsi)

上述十二张卡分别驱动 M0、两套 M1、M2.1、M3.1、四套 M2.2、M2.3、
M3.3 和 M4.1 的
fuelsim-to-MOOSE 对比；测试程序不再直接构造这些案例的材料、载荷路径或
网格选择参数。每个对比读取 MOOSE 最终时刻的全部节点，统一检查温度、
径向位移和轴向位移的三项误差；M1、M2.3、M3.3 和 M4.1 还检查全部接触
节点的压力三项误差。
