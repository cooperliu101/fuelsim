# fuelsim 输入卡 v3

`fuelsim` 只使用一个显式输入文件：

```bash
./build/fuelsim -i case.fsi [PETSc options]
```

输入卡采用 MOOSE 风格的嵌套段和 `key = value`。`#` 开始行内注释；相对
路径以输入卡所在目录为基准。v3 的数值全部使用 SI 单位，不支持 include、
宏、表达式、单位换算、旧键别名或兼容层。未知段、未知键、重复项、非法值
和缺少必填键都会立即报错。

每张卡都必须显式声明格式版本和物理问题：

```text
[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
```

## 问题与单一网格

`[Case]` 的 `problem` 只能是：

- `steady`：稳态热传导和准静态热弹性，使用线性载荷步；
- `transient`：Backward Euler 热传导、准静态力学和积分点非弹性历史。

对应的两个生产 C++ 类型是 `SteadyProblem` 和 `TransientProblem`。M0、M1、
M2 只作为路线和回归名称。

`geometry` 必须显式选择 `axisymmetric_rz` 或 `cartesian_3d`。前者要求二维
Quad4 网格，使用 `[T(:), ur(:), uz(:)]`；后者要求三维 HEX8 或 HEX20 网格，
使用 `[T(:), ux(:), uy(:), uz(:)]`。HEX20 采用二阶 20 节点位移和一阶八角点
温度，同一文件不能混合 HEX8 与 HEX20。版本 2 和省略几何的输入都会被拒绝。

`[Mesh]` 只接受一个 Exodus 文件：

```text
[Mesh]
  type = exodus
  file = model.e
[]
```

文件可以包含多个 Quad4、HEX8 或 HEX20 元素块、节点集和边集；一个输入卡不能
混合这些拓扑。输入输出层保留这些元数据，
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

## 组合式材料函数

`[Materials]` 下的每个子段定义一个可被多个区域引用的材料。材料由热物性、
弹性、本征应变、蠕变和塑性函数组合。热物性和弹性函数必须存在；本征应变
可以有任意多个具名实例，蠕变和塑性函数可以省略：

```text
[Materials]
  [fuel]
    [thermal]
      function = inverse_temperature_thermophysical
      conductivity_inverse_temperature = 3824
      conductivity_constant = 0.61
      density = 10970
      specific_heat = 300
    []

    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.316
    []

    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []

    [creep]
      function = norton
      coefficient = 1e-20
      reference_stress = 1e8
      stress_exponent = 4.5
    []

    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 2e8
      hardening_modulus = 1e9
    []
  []
[]
```

每个注册函数同时声明严格的参数名称和 SI 单位。解析器先读取 `function`，再按
该函数的参数表验证本段其余键。所有参数必须显式给出；未知参数、缺少参数、
重复参数和非有限值都会报告输入路径及行号。单位字符串用于接口说明和错误
信息，不执行单位换算。

内置函数包括：

- 热物性：`constant_thermophysical`、`inverse_temperature_thermophysical`；
- 弹性：`constant_isotropic`、`linear_temperature_isotropic`；
- 本征应变：`isotropic_thermal_expansion`、
  `linear_temperature_isotropic_thermal_expansion`；
- 蠕变：`norton`、`linear_temperature_norton`；
- 塑性：`linear_isotropic_hardening`、
  `linear_temperature_isotropic_hardening`。

本征应变表示材料在无外力时产生的应变，例如热膨胀、肿胀和致密化。同一材料
的多个本征应变实例逐分量相加。蠕变函数返回等效蠕变速率，塑性函数返回给定
等效塑性应变下的流动应力。fuelsim 统一执行 J2 关联流动、Backward Euler
局部积分、塑性—蠕变全隐式耦合和有限应变客观旋转。

用户可创建 `MaterialFunctionRegistry`，注册普通 C++17 函数，再调用
`read_case_input(path, registry)` 读取引用这些函数的输入卡。函数的活跃
输入和输出使用具体类型 `adlite::Scalar`，因此温度和力学链式导数进入局部
Jacobian。注册函数采用两阶段接口：绑定函数在读取材料时按具名参数读取一次，
返回捕获已绑定参数的求值函数；积分点只调用求值函数，不再执行参数字符串查找。
函数必须是无副作用的纯函数，不能保存积分点 trial 状态或访问 PETSc 和全局解向量。

当前组合接口限定为各向同性弹性和各向同性本征应变。轴对称后端还支持 J2
关联塑性及沿最终 J2 方向的等效蠕变；三维阶段 B 只接受弹性材料。它不声明
各向异性、运动硬化、损伤、非 J2 屈服面或非关联流动已经受支持。

## 自由区域组合

`[Regions]` 下每个子段定义一个物理区域，子段名是区域名。每个区域必须
用且只能用 `block` 或 `block_id` 选择同一 Exodus 文件中的元素块；有块名
时优先使用更易读的 `block`，无块名的 MOOSE 单区域网格可使用非负
`block_id`：

```text
[Regions]
  [pellet]
    block = fuel
    material = fuel
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2e8
    heat_source_function = power
  []
[]
```

区域数量不固定，因此同一结构可表示单独芯块、单独包壳、芯块—包壳，或
芯块—包壳1—包壳2。每个元素块只能声明一次。三维模式中，不同元素块若在
Exodus 网格中引用同一源节点，会原生映射到同一套全局温度和位移自由度；不
共享的区域节点仍保持独立。二维轴对称区域仍保持各区域独立节点的规则。

每个区域必须显式设置 `strain = small` 或 `strain = finite`。`small` 使用
参考构形小应变弱式。瞬态 `finite` 使用 MOOSE 默认的增量 Taylor 应变与
Rashid 转动；稳态没有 committed 材料历史，从参考构形 `F_old=I` 对当前总
变形做一次 Taylor 更新，不能解释为随稳态载荷步累计的增量材料路径。两者
都用 Cauchy 应力、当前构形梯度和当前体积测度装配力学内力。热传导与热容
仍使用参考构形。pressure 的构型由边界上的 `configuration` 选择；省略该字段时，程序根据所属区域的应变形式
自动采用推荐值：小应变为参考构形、有限应变为当前构形。显式选择非推荐组合时会提示，但仍按输入执行。
traction 的构形选择规则见下文。区域发生非正 Jacobian、非正环向
伸长或非正当前半径时会拒绝 Newton 试探态，不做隐式夹持。

每个区域必须用 `material` 引用 `[Materials]` 中已经定义的材料。旧版把导热率、
弹性和非弹性参数直接写在区域内的格式不再接受。

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
      law = gas_gap
      gap_conductivity = 0.4
      minimum_gap = 1e-6
    []

    [mechanical]
      formulation = penalty
      discretization = surface_to_surface
      sliding = small
      penalty_factor = 0.25
      mu = 0.3
    []
  []
[]
```

一个接触对至少包含 `[thermal]` 或 `[mechanical]`，也可以同时包含两者。
热接触省略 `law` 时使用 `gas_gap`，其导热系数为：

```text
h = gap_conductivity / max(g, minimum_gap)
```

需要对齐 Abaqus 的间隙导热表局部斜率时，可以选择线性仿射定律：

```text
[thermal]
  law = affine
  conductance = 50
  clearance_derivative = -1000
  pressure_derivative = 2e-6
  temperature_derivative = 0.1
  reference_temperature = 350
[]
```

对应的积分点定律为：

```text
T_average = (T_secondary + T_primary) / 2
p = penalty * max(-g, 0)
h = conductance
  + clearance_derivative * g
  + pressure_derivative * p
  + temperature_derivative * (T_average - reference_temperature)
q = h * (T_secondary - T_primary)
```

`conductance` 的单位为 `W/(m2 K)`，三个导数分别相对于间隙、压力和平均温度。
压力导数非零时必须同时定义罚函数机械接触，使热学和力学使用同一个接触压力。
参考状态的 `conductance` 必须非负；试探态计算得到负值或非有限值时，程序把它
作为物理域错误交给线搜索或时间步缩小重试，不会夹持为零。该仿射定律准确表示
Abaqus 间隙导热表中的一个线性单元；当前没有实现任意表格的分段插值、外推和
截断语义。

机械接触的 `formulation` 必须显式选择 `penalty` 或
`augmented_lagrangian`。`penalty` 的单位为 `Pa/m`；如果省略，程序用两侧
边界单元的材料刚度和法向网格尺度自动计算。每个边界相邻 Quad4 的法向尺度
为参考平面面积除以边长，每侧取最小值，并按下式组合：

```text
k_interface = 1 / (h_primary / E_primary + h_secondary / E_secondary)
penalty = penalty_factor * k_interface
```

`discretization` 可显式选择 `node_to_surface` 或 `surface_to_surface`。省略时，
HEX20 机械接触采用 `surface_to_surface`，HEX8 和轴对称 RZ 采用
`node_to_surface`。三维 HEX8 和 HEX20 的 `surface_to_surface` 当前只支持罚函数形式；
HEX8 必须显式选择该离散。HEX8 的 `sliding = small` 只允许两侧均为小应变，
`sliding = finite` 可用于小应变或有限应变区域。轴对称 RZ 机械接触仍只采用节点到线段离散。
`sliding` 可选择 `small` 或 `finite`，省略时为 `small`；该键只允许与显式的
`discretization = surface_to_surface` 同时使用。小滑移始终保留参考构形确定的
primary 面归属。对于两侧均为小应变的线性罚接触，小滑移采用 Abaqus 对标识别出的
节点中心约束。每个 HEX8 secondary QUAD4 面形成四个等正面积约束，其中心为父坐标
`(±0.5,±0.5)`，固定平均矩阵为：

```text
A4 = (1/16) * [9 3 1 3; 3 9 3 1; 1 3 9 3; 3 1 3 9]
```

每个 HEX20 secondary QUAD8 面形成八个正面积约束，使用已经识别的二次平均矩阵；
四个角点约束的中心为 `(±0.75,±0.75)`，四个中边节点约束的中心为 `(0,±0.5)`
和 `(±0.5,0)`。每个约束在相应父坐标位置计算独立法向，因此多个平面面片能够逐面
改变法向，真正的二次曲面也能在同一个面内改变法向。共享 secondary 节点会汇总
相邻面的约束面积和法向。

问题构造时仍把 secondary 面按参考态 primary 面所有权递归分片，用一致的测试函数
把每个节点中心约束守恒地投影到可能非匹配的 primary 面。每个约束的罚残量采用同一
间隙梯度向两侧装配作用力和反作用力，切线为该梯度的外积。该规则不是经典 mortar
离散，也不引入 mortar 乘子。两侧残量满足：

```text
R_secondary_i = integral(N_secondary_i * traction dA)
R_primary_j   = -integral(N_primary_j * traction dA)
```

当 `mu > 0` 时，同一个节点中心约束还在其固定法向的切平面内累计
secondary 相对 primary 的增量位移，并执行三维 Coulomb 粘着或滑移返回。输入使用与
Abaqus/Standard `*Friction, slip tolerance=` 相同语义的无量纲 `slip_tolerance`。程序以
secondary 各接触面参考面积的算术平均值 `A_mean` 计算特征长度 `sqrt(A_mean)`，并形成
最大弹性滑移距离 `delta_e = slip_tolerance*sqrt(A_mean)`；粘着刚度为
`mu*p/delta_e`。省略或设为零时，
`surface_to_surface` 接触采用 Abaqus 默认值 `0.005`。该状态按每个节点中心约束进入
提交、回滚和检查点事务。
`slip_tolerance` 允许用于 HEX8 或 HEX20 `surface_to_surface` 接触。两侧小应变且选择
小滑移时使用上述固定切平面的节点中心约束；任一侧采用有限应变的小滑移继续使用参考分片后的
3×3 积分路径只适用于 HEX20。HEX8 小滑移面对面接触若任一侧采用有限应变，会在问题构造时
明确报错。

有限滑移对小应变和有限应变使用同一个当前构形表面积分算法。HEX8 使用四个
父坐标 `(±0.5,±0.5)` 节点中心点和 QUAD4 当前几何；HEX20 使用 3×3 积分点和
QUAD8 当前几何。每次状态验证都为每个 secondary 积分点在完整 primary 面集合中
重新搜索唯一最近投影；超过
64 个 primary 面时使用可重整的空间搜索树。跨内部边时只保留一个所有者，积分点滑出
完整 primary 表面时明确拒绝当前 Newton 状态。构造期仍为每个潜在面配对预留稀疏
耦合，因此跨面不会改变矩阵非零结构。

有限滑移摩擦以 primary 面当前构形的随动正交切向基保存历史。弹性滑移先从已提交的
接触切向基客观运输到当前切平面，再用当前投影点和已提交接触点之间的曲面滑移增量
更新；共同刚体转动不会被误算为新增滑移。法向、两个切向分量和表面测度均进入
ADlite 一致切线。结果文件中的节点压力仍是恢复量；有符号等效节点量只用于输出，
不作为独立节点罚刚度。有限滑移目前只用于三维 HEX8 或 HEX20 的
`surface_to_surface` 罚接触，不用于轴对称 RZ、`node_to_surface` 或增广拉格朗日接触。

H20.30 使用同一组 Exodus 坐标分别比较平面面片拼成的四分之一圆柱和中边节点位于
真实圆弧上的二次圆柱。真正二次曲面算例的径向位移三项误差分别为 `0.113348%`、
`0.063526%` 和 `0.289203%`，节点径向接触力三项误差分别为 `0.047713%`、
`0.056410%` 和 `0.106361%`，均小于 `1%`；径向合力误差为 `0.0000723%`。

`quad8_nodal_area_rule` 只适用于显式选择的 HEX20 `node_to_surface` 对比路径。
`positive_lumped` 先计算每个节点的平方形函数积分，再按当前面面积归一化：

```text
raw_area_i = integral(N_i * N_i dA)
nodal_area_i = face_area * raw_area_i / sum(raw_area)
```

因此每个节点面积严格为正，并且八个节点面积之和等于当前接触面面积。对于平直、
规则的单位 QUAD8 面，四个角点面积各为 `3/76`，四个边中点面积各为 `4/19`。
可显式选择 `consistent_shape`，保留原来的 `integral(N_i dA)` 有符号面积，仅用于
与 MOOSE 的传统 node-face `MechanicalContactConstraint` 且
`normalize_penalty = true` 的结果对比。该旧规则在规则 QUAD8 面的角点面积为
`-1/12`、边中点面积为 `1/3`，所以角点不具有正的 Coulomb 摩擦容量。
`consistent_shape` 只允许用于 HEX20 `node_to_surface`；二维 RZ、HEX8 和
`surface_to_surface` 接触会在问题构造时拒绝它。

H20.19 另用 MOOSE 的双基函数 mortar 面积分作为独立排序参考。在同一二单元纯法向
压缩算例中，`consistent_shape` 的法向位移和合力比 `positive_lumped` 更接近 mortar，
但这不等于允许把负的一致节点面积直接乘入 node-face 罚刚度。mortar 在面分段上积分
分布式约束，并不依赖这个有符号节点罚刚度。新的 HEX20 生产默认
`surface_to_surface` 的法向位移相对 L2、相对绝对峰值、最大逐点相对误差和合力误差
分别为 `0.320104%`、`0%`、`0.657040%` 和 `0.507750%`，四项均在该离散差异算例的
`1%` 门槛内。

H20.21 又以 Abaqus/Standard 全积分 C3D20 surface-to-surface 罚接触复核同一终态
闭合量。`consistent_shape` 的法向位移相对 L2、相对绝对峰值、最大逐点相对误差
和法向合力误差分别为 `0.073786%`、`0.000000834%`、`0.184464%` 和
`0.001608%`；`positive_lumped` 对应为 `11.0014%`、`0.000000834%`、
`28.8039%` 和 `7.74441%`。新的 `surface_to_surface` 路径对应为
`0.0422346%`、`0.000000834%`、`0.105816%` 和 `0.000931639%`。Abaqus 与
MOOSE mortar 都支持生产面约束路径，同时仍证明不了负节点面积的 node-face 罚刚度安全。

H20.23 在同一 Abaqus C3D20 模型上增加 `mu = 0.001` 和切向位移，并把两边的
`slip_tolerance` 都设为 `1e-6`。该接触面的特征长度为 `0.01 m`，程序内部得到
`delta_e = 1e-8 m`，八个节点中心约束都进入滑动。法向、切向和离开对称面的第三
位移分量三项误差中的最大值分别为 `0.00000464%`、`0.00000413%` 和
`0.00000472%`；理论零值对称面上的八个 Abaqus 舍入噪声点不进入相对误差，改用
绝对差单独核算，其最大绝对差为 `1.1352e-18 m`。法向和切向合力误差分别为
`0.000000589%` 和 `0.0255274%`。局部测试还
覆盖粘着与滑动切线、历史
提交与回滚以及两侧残量严格守恒。

H20.33 在 H20.24 的非匹配平面网格上增加 `mu = 0.3` 和七步反向路径，逐步比较
三个位移分量、带符号法向与切向节点力、两个切向滑移分量以及法向与切向合力。
路径覆盖全粘着、同一时刻粘着与滑动并存、正向滑动、卸载、反向滑动和重新粘着；
第四步后的检查点重新启动与不间断路径在后续每步的节点状态和全部接触历史上完全
一致。非匹配转变期间法向位移的小参考值使用显式 `12%` 逐点门槛，节点力使用显式
`1.5%` 逐点门槛；聚合误差仍小于 `1%`，并且没有增加分母下限。

H20.35 在真正二次圆柱面上同时施加周向和轴向切向运动，37 个约束全部保持粘着。
三个位移分量、带符号径向法向节点力、周向与轴向切向节点力、两个 Abaqus 局部切向
滑移分量和三个合力均逐项对比。除一个 `-6.26679e-8 m` 的 X 位移参考点采用显式
`1.5%` 逐点门槛外，其余三项场误差均小于 `1%`；该点的绝对差为
`8.54026e-10 m`，零参考值仍单独统计。Abaqus 数据文件报告的特征接触长度为
`0.44311 m`，因此两边直接输入 `slip_tolerance = 1e-5`，程序内部换算出的
`delta_e = 4.4311e-6 m`；它不是拟合得到的材料系数。

H20.36 在同一真正二次圆柱面上施加 `2, 4, 12, 24, 4, -20, -18 um` 的七步轴向
位移路径，覆盖全粘着、粘着与滑动并存、正向滑动、卸载、反向滑动和全部重新粘着。
37 个约束始终闭合；检查点重新启动与不间断路径逐项一致，粘着与滑动接触 Jacobian
方向导数误差分别为 `2.15e-11` 和 `2.59e-10`。径向位移、轴向节点力和近零轴向
滑移分别使用只限本算例的 `1.4%`、`2%` 和 `6%` 逐点鉴定门槛，实际最大值为
`1.32187%`、`1.84225%` 和 `5.45628%`；对应聚合误差均小于 `1%`，没有增加分母
下限。未施加驱动的周向场仍输出三项相对误差，但由于局部参考值穿过零，使用相对
轴向驱动峰值的最大绝对差验收，门槛为 `1%`。该算例只鉴定曲面上的单一切向驱动。

H20.38 使用同一追踪 Exodus 坐标建立两个 C3D20 单元的受控表面到表面接触。第一步在
两个切向分量同时非零的条件下使全部九个 Fuelsim 表面积分点进入滑动；第二步把已经
压入和滑动的完整接触构形共同刚体转动 `0.35 rad`，不增加物理相对滑移。Fuelsim 的
法向合力对 Abaqus 的相对 L2、相对绝对峰值和最大逐点相对误差分别为
`0.0035795%`、`0.00506231%` 和 `0.00506231%`；双切向合力对应为 `0.0262916%`、
`0.0275768%` 和 `0.0490877%`。Fuelsim 合力转动最大绝对误差为 `1.85e-10 N`，双分量
弹性滑移历史转动最大绝对误差小于 `9.9e-17 m`，两侧残量不平衡小于 `8.6e-14 N`。
Abaqus 自身的法向和切向合力转动三项误差也都小于 `0.1%`，两个非零总滑移分量在转动
前后最大变化为 `2.32e-8 m`。另外，局部内核在真正二次圆柱面上验证相同的双分量历史
转移及其自动微分切线。这个组合证据鉴定小滑移锚点随面转动的客观性，但不外推为任意
大滑移、任意曲率或 Abaqus 专有节点平均算子的有限应变逐项等价。

`penalty_factor` 是无量纲可选值，默认 `1`，必须有限且大于零。显式 `penalty`
和 `penalty_factor` 互斥，不能同时出现。自动选择只是网格与材料一致的起点；
生产工况仍须用穿透、接触力和网格收敛证明其适用性。

`augmented_lagrangian` 在固定乘子下完成一次 Newton 求解，再更新节点法向乘子，
直到所有被捕获节点的绝对间隙小于 `penetration_tolerance`。该容差默认
`1e-8 m`，必须有限且大于零；`maximum_augmented_iterations` 默认 `20`，必须
大于零。达到上限仍未满足约束时，本次载荷步或时间步失败并完整回滚。
这两个字段只允许用于 `augmented_lagrangian`，在 `penalty` 形式中会被拒绝。
乘子仅在内层 Newton 收敛后更新，瞬态提交、失败重试、step-doubling 和检查点
均把它作为事务状态处理。

`mu` 是可选的 Coulomb 摩擦系数，必须为非负有限值，默认值为 `0`。默认值
保持原无摩擦残量和解逐项不变。`mu > 0` 时，切向罚刚度与最终确定的
`penalty` 使用同一个 `Pa/m` 数值；节点到面路径在 secondary 节点、表面到面路径在
3×3 积分点用本步相对切向位移形成弹性预测牵引，再将其限制在
`mu * pressure`。上限以内为粘着，达到上限并继续
同向运动时为滑移，反向运动可重新进入粘着。弹性切向滑移和粘滑标志只在
收敛载荷步或时间步提交，失败重试从同一 committed 状态重算；瞬态检查点
保存这两个量，并且不读取旧检查点格式。

轴对称 RZ 接触的两侧必须来自不同区域且各自形成一条不分叉的开放边链。热接触
按参考投影重叠区间切分 secondary-side STS 积分，机械接触采用 secondary 节点
到 primary 线段的唯一 NTS 投影；二维法向同时装配径向和轴向反力。构造时分别
为每个热积分点和每个机械 secondary 节点预留整条 primary 链的潜在稀疏耦合，
残量和 Jacobian 评估前按当前构形选择距离最近的唯一有效线段，因此两者都可以
跨越任意数量的链内线段而不重建 PETSc 工作区。内部顶点使用半开区间，任何时刻
只允许一条 primary 段拥有同一积分点或节点；整条链的首端和末端可以归属其端点。
热积分点滑出完整 primary 链时，Newton 试探状态会作为物理域错误交给回溯线搜索；
若仍无法恢复则拒绝当前载荷步或时间步，不会夹持到链端、静默返回零热流或继续
使用陈旧候选。机械接触仍保留参考链首尾节点的物理端点支承，其他机械节点失去
全部有效投影时同样拒绝状态。动态候选段由当前几何确定，不写入检查点。唯一活动
热候选向两侧装配严格相反的残量，保证离散热守恒。

三维笛卡尔接触同时支持 HEX8 四节点面和 HEX20 八节点二次面。热接触都在
secondary 面的 2×2 四个积分点上计算；HEX20 温度仍只使用四个角点的一阶形函数，
但当前面坐标和投影使用八节点二次几何。HEX8 默认机械路径采用四个 secondary 面
节点到 primary 面的正交投影；显式选择小滑移 `surface_to_surface` 时采用上述四节点
中心平均约束。选择有限滑移时，无摩擦平行 primary 面采用当前构形节点中心平均约束，
摩擦或非平行 primary 面采用四个当前构形节点中心积分点。HEX20 默认采用小滑移
节点中心平均约束，并保留八节点
`node_to_surface` 显式对比路径；有限滑移和有限应变小滑移使用 secondary 面 3×3
积分点到 primary 二次面的正交投影。需要动态搜索的每个点都预留该接触
对全部 primary 面候选，并在当前构形中选择唯一有效面。内部公共边只允许一个面
拥有投影，投影跨边时所有权唯一转移。无摩擦 HEX8 有限滑移平均表面约束滑出完整
primary 表面时自然释放并贡献严格零残量；热接触、节点到面、小滑移、摩擦有限滑移及
HEX20 积分路径仍拒绝失投影的当前 Newton 状态。热流、法向力和三维切向力均向两侧
装配严格相反的贡献。三维 Coulomb 摩擦
保存全局三分量切向弹性滑移向量，因此可以表示接触面的两个独立切向方向。当前
三维机械接触只接受 `formulation = penalty`；选择 `augmented_lagrangian` 会在
问题构造时明确报错。小应变、小滑移表面到面路径按每个唯一 secondary 节点保存一份
法向和摩擦历史；无摩擦平行 primary 面的 HEX8 有限滑移平均路径同样按唯一 secondary
节点保存状态，摩擦或非平行 primary 面路径按每个 secondary 面四个点保存历史，HEX20
的积分路径按九个积分点保存历史。输出层再恢复等效节点力和
节点压力。正集总与有符号一致面积只保留给显式 `node_to_surface` 对比路径。

接触两侧允许零初始间隙：参考构形中 secondary 节点可以恰好骑在 primary
线段上（例如初始贴合的芯块—包壳），构造不再要求处处为正的参考间隙。
此时法向朝向无法由间隙符号决定，程序改用单元材料侧拓扑：计算 secondary
边所属 Quad4 父单元质心相对 primary 线段的有符号位置（沿用间隙的基准
法向 `(tangent_z, -tangent_r)/length`），并取法向背离 secondary 材料、
指向 primary 一侧，使节点向 secondary 材料内部移动时间隙增大。该约定与
把同一几何的间隙打开任意小量后按间隙符号得到的朝向完全一致。若零间隙下
两侧父单元质心位于 primary 线段同侧（材料重叠的畸形网格），或质心恰好
落在线上（退化单元），构造仍会明确报错，不做隐式夹持。

Exodus 节点结果除间隙和法向压力外，还输出接触切向牵引、弹性切向滑移和
以 `0/1` 表示的滑移标志。轴对称结果中的切向量按有序 Line2 边链方向保留
符号；三维结果输出切向牵引与弹性切向滑移的向量模。全局结果与工程历史另输出
secondary 侧切向合力；轴对称为有符号标量，三维为合力向量的模。

`[Contact]` 是可选段。只有定义至少一个接触对时才需要写出该段；不含接触的算例应省略它。

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

轴对称 `dirichlet` 的 `field` 可为 `temperature`、`radial_displacement` 或
`axial_displacement`；三维可为 `temperature`、`displacement_x`、
`displacement_y` 或 `displacement_z`。`pressure` 不接受 `field`，方向取父单元
外法向。压力边界可以在小应变和有限应变区域中分别设置
`configuration = reference` 或 `configuration = current`；前者使用参考表面，后者使用当前法向和当前表面测度，
轴对称区域还会使用对应的参考半径或当前半径。小应变区域推荐参考构形，有限应变区域推荐当前构形；选择
其他组合时程序会给出提示，但仍按用户选择装配。
`traction` 必须声明一个位移 `field`，`value` 是该全局位移分量上的有符号表面牵引；也可以在小应变和有限应变
区域中设置 `configuration = reference` 或 `configuration = current`；省略该字段时同样按区域应变形式自动采用小应变参考构形或有限应变当前构形。当前构形仍固定为所选全局分量方向，
但周长和边长使用当前构形并进入 AD Jacobian；小应变推荐参考构形，有限应变推荐当前构形，非推荐组合会提示
但仍按输入执行。
`dirichlet`、`pressure` 和 `traction` 都可设置
`scale_with_load = true`，使 `value` 乘以当前执行器载荷因子；默认不缩放。
也可用 `function = <name>` 使 `value` 乘以时间表值；`function` 与
`scale_with_load` 互斥。压力的物理范围由用户负责，程序不强制检查。

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

两个函数均可省略；存在时分别乘以对应基值。换热系数和环境温度的物理范围由
用户负责，程序不强制检查。
对流项使用参考表面测度。三维 HEX20 的一阶温度对流边界采用 2×2（4 点）面内
积分；同一二次面上的压力和分量牵引仍采用 3×3（9 点）积分。轴对称边界作为
12 自由度贡献装配，三维四节点面作为 16 自由度贡献装配，HEX20 二次面作为 28
自由度贡献装配；这些贡献的残量和温度
切线都由 ADlite 保持一致。
`[BoundaryConditions]` 是可选段。没有边界条件时应省略它。

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
  strain_history_time_absolute_tolerance = 1e-10
  stress_history_time_absolute_tolerance = 1
  time_error_safety_factor = 0.9
  include_thermal_time_term = true
  restart = previous.checkpoint
[]
```

`include_thermal_time_term` 默认为 `true`。设为 `false` 时，瞬态温度方程不加入
`rho*cp*(T_new-T_old)/dt` 热容时间项，但仍保留热传导、热源和力学瞬态材料更新。

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
原行为，此时窗口必须为零。窗口必须小于目标。若当前接受步经历过 cutback，
下一步最多保持接受后的实际步长，不在同一步成功后立即放大回刚失败的尺度；
后续无拒步的成功步才重新允许增长。

`time_error_relative_tolerance` 省略或为零时不做时间离散误差控制。设为正值
后，每个候选步从同一 committed 状态计算一个 Backward Euler 全步和两个
半步；温度、两个位移场、弹性/塑性/蠕变张量、两个等效应变及应力的归一化
L2 差最大值大于 1 时完整回滚并缩步，成功时采用两个半步的结果。温度、位移、
应变历史和应力分别使用上述绝对容差，安全系数必须位于 `(0,1)`。该估计器
增加到约三倍的非线性求解工作量，但可识别节点场不敏感而材料历史不准确的
时间步，并在进度和拒步诊断中输出总估计及每个分量。接受步的热率诊断取两个
半步的时间平均，功、能量变化和耗散取两个半步之和，因此对应完整控制步，
不只对应第二个半步。

生产默认保持 `time_error_relative_tolerance = 0`。原因不是误差估计器不可靠，
而是它约需三倍非线性求解工作量，且位移、应变和应力绝对容差必须根据工况的
物理尺度确定，不能由程序给出通用值。固定步长 MOOSE 验证输入也因此保持原样。
对长时蠕变、塑性累积或接触状态变化工况，建议先用固定步长做至少三档全局
收敛研究，再以 `1e-3` 作为相对容差的起始筛选值，并为温度、位移、应变和
应力分别填写有物理意义的绝对容差；若接受解不能满足工况自己的全局误差目标，
继续收紧相对和绝对容差。`1e-3` 只是起始建议，不是精度保证或默认门槛。

每次未收敛尝试都会记录尝试终点、步长、cutback 序号、非线性/KSP 迭代数、PETSc
收敛原因、残量范数、失败类别和物理域消息；最终停止原因区分 `completed`、
`maximum_cutbacks` 和 `minimum_time_step`。拒绝步仍完整回滚 committed
状态；最小时间步会实际求解一次，只有该次也失败才终止。

`restart` 为可选的严格重启动文件。它恢复已提交的节点温度/位移、物理时间、
载荷因子、全部积分点塑性/蠕变历史和已提交应力。文件的版本、字节序、长度、
校验和、网格、材料、边界条件、接触及局部装配拓扑必须与当前问题一致；不
匹配时立即停止。重启动不保存 Newton trial、活动时间步或失败尝试。
当前格式 v7 还保存几何类型、接触摩擦历史、法向增广乘子、成功提交后控制器给出的下一名义时间步和
最后一个完整接受步的守恒/耗散摘要，因此自适应计算从检查点继续时不会重新使用输入卡的初始步长，
零步重启动结束也不会把上一接受步诊断伪装成全零。版本 6 及更早格式会被
明确拒绝，不提供跨版本兼容层。

`[Solver]` 是可选段；省略时使用下列全部默认值。存在时可设置非线性 `absolute_tolerance`、`relative_tolerance`、
`step_tolerance` 和 `maximum_iterations`；省略时分别为 `1e-8`、`1e-10`、
`1e-12` 和 `40`。线性选项为：

- `linear_solver = automatic|direct|gmres`；
- `preconditioner = automatic|lu|block_jacobi|field_split|hypre`；
- `direct_factorization = automatic|mumps`；默认 `automatic` 在单 rank 使用 PETSc
  LU、多 rank 使用 MUMPS，`mumps` 则在所有进程数明确选择 MUMPS；
- `linear_relative_tolerance`，默认 `1e-8`；
- `maximum_linear_iterations`，默认 `500`；
- `jacobian_lag`，默认 `1`；设为大于一的整数时，在该数量的非线性迭代内复用
  已装配的 Jacobian 矩阵；
- `backtracking_fallback`，默认 `true`；BASIC 失败后从原始初值用 BT 重试；
- `residual_reduction_tolerance`，默认 `1e-6`，用于总残量和分场残量复核；
- `temperature_residual_absolute_tolerance`，默认 `1e-8 W`；
- `mechanical_residual_absolute_tolerance`，默认 `1e-4 N`，同时用于径向和轴向；
- `field_residual_scaling`，默认 `false`，可选启用热/力分组的自动行缩放。
- `temperature_residual_scale` 与 `mechanical_residual_scale`，默认均为 `0`；
  成对设为正数时作为跨求解固定的物理残量特征尺度，与自动行缩放互斥。

`automatic` 使用直接 LU；具体分解器由 `direct_factorization` 选择。选择
`block_jacobi`、`field_split` 或 `hypre` 会自动选 GMRES；
`field_split` 按问题提供的热学字段和全部力学字段建立乘法场分裂。具体 PETSc
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
PETSc 若以 `MAX_IT`、`LINE_SEARCH` 或 `LOCAL_MIN` 等负原因停止，fuelsim 只在
重新计算的总残量和三个分场残量都通过同一套复核时才接受该状态；原始负原因
仍保留在诊断中。输出同时报告本次和全程 KSP 迭代数。

工程规模选择证据目前只覆盖同一 M1 物理模型的 23,010 和 45,630 自由度
结构网格、两个 MPI 进程以及当前 PETSc 构建。两种规模的 20 步工况中，直接
MUMPS 和未缩放 HYPRE 都完成求解；预热内部载荷路径时间分别为
`10.0455/21.3292 s` 和 `87.6334/84.4195 s`。HYPRE 相对直接分解慢
`8.72/3.96` 倍，累计线性迭代分别为 `7,010/2,960`，而直接分解均为
`62`。因此，在这两个已测规模上，只要直接分解的内存可接受，`automatic`
或显式 `direct + lu` 是有实测支持的首选；不能仅因网格变大就假定 HYPRE
更快。

同一 23,010 自由度工况的两步筛查中，块 Jacobi 和乘法场分裂无论是否开启
自动分场残量缩放，都未完成第一个载荷步；累计线性迭代为
`1,270–2,600`，末次残量为 `0.685–0.906`。它们保留为需要针对具体问题
调参和完整验证的实验选择，不是工程默认。自动分场缩放曾把 HYPRE 两步筛查
从 557 次线性迭代降到 435 次，但完整 20 步在第 14 步失败，所以这个短程
改善也不能作为推荐设置。选择任何迭代组合时，必须检查
`load_steps_completed`、最终总残量和分场残量、非线性与线性迭代数，并在目标
网格和完整载荷路径上与直接解核对；`benchmarks/README.md` 给出完整命令和
首轮/预热计时。

多进程残量和 Jacobian 回调只收集本进程贡献以及完整主面链接触搜索依赖所需
的影子自由度，不再每次复制完整试探态。诊断中的 `global_state_dofs`、
`maximum_shadow_state_dofs`、`total_shadow_state_dofs` 和
`total_remote_shadow_state_dofs` 分别给出全局规模、单进程最大影子规模、所有
进程影子槽总数和每次回调所需的远程值总数。求解结束仍有一次完整状态收集，
以保持 replicated committed 状态、材料历史和输出事务；因此这些指标只描述
回调通信和明确的影子缓冲区，不代表总进程内存。

求解诊断还输出 `memory.*` 和 `aggregate_memory.*` 的驻留内存字段，单位为字节，
由 PETSc 的当前/最大内存接口采样。`initial_resident_bytes`、
`setup_resident_bytes`、`solve_resident_bytes` 和 `final_resident_bytes` 分别是
进入本次求解、PETSc 工作区设置完成、`SNESSolve` 返回和最终状态收集完成时的
最大 rank 当前驻留集；`minimum_peak_resident_bytes`、
`maximum_peak_resident_bytes` 和 `total_peak_resident_bytes` 是本次求解期间各
rank 峰值驻留集的最小值、最大值和总和。瞬态的 `aggregate_memory.*` 取所有
时间步和重试中观察到的最大阶段值，适合判断直接分解的内存上界。

`[Outputs]` 是可选段；省略时 `console` 默认为 `true`，且不写文件。可选 `csv` 将最终命名指标写为
`metric,value` 汇总文件。周期性的 `progress.*` 只写控制台，不混入最终 CSV。
`exodus` 写出可后处理的场结果；稳态写一个最终步，
瞬态按 `exodus_interval` 写初始或重启动状态、成功提交步及最终状态。轴对称
节点变量包括温度、径向与轴向位移，以及各接触对 secondary 节点上的间隙和
压力；单元变量保留四个积分点的应力与非弹性历史。三维节点变量为温度、三个
笛卡尔位移，以及每个接触对 secondary 节点上的间隙、压力、切向牵引模、弹性
切向滑移模和滑移标志；单元变量为八个积分点的六分量应力。全局变量记录载荷
因子、界面总热流和总接触力。未属于所选求解区域或未投影的值写为 `NaN`。

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
最大压力、总热流和总反力。每步还记录生成/储存/对流/界面/Dirichlet 热率、
全局热平衡、内力/压力牵引/约束反力/接触功平衡、弹性能变化及塑性/蠕变耗散。
检查点只在成功提交后按间隔原子替换；若该间隔写入已经是最终提交态，结束时
不重复写同一文件，否则在结束或失败退出前补写最后提交态。

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
- [`steady_hex8_thermoelastic.fsi`](../verification/fuelsim/steady_hex8_thermoelastic.fsi)
- [`steady_fuel_cladding.fsi`](../verification/fuelsim/steady_fuel_cladding.fsi)
- [`steady_fuel_cladding_unstructured.fsi`](../verification/fuelsim/steady_fuel_cladding_unstructured.fsi)
- [`steady_augmented_contact_moose.fsi`](../verification/fuelsim/steady_augmented_contact_moose.fsi)
- [`steady_two_pellet_contact_moose.fsi`](../verification/fuelsim/steady_two_pellet_contact_moose.fsi)
- [`transient_heat_moose.fsi`](../verification/fuelsim/transient_heat_moose.fsi)
- [`transient_table_convection_moose.fsi`](../verification/fuelsim/transient_table_convection_moose.fsi)
- [`transient_j2_plastic_moose.fsi`](../verification/fuelsim/transient_j2_plastic_moose.fsi)
- [`transient_j2_unload_reload_moose.fsi`](../verification/fuelsim/transient_j2_unload_reload_moose.fsi)
- [`transient_norton_creep_moose.fsi`](../verification/fuelsim/transient_norton_creep_moose.fsi)
- [`transient_coupled_displacement_moose.fsi`](../verification/fuelsim/transient_coupled_displacement_moose.fsi)
- [`transient_coupled_traction_moose.fsi`](../verification/fuelsim/transient_coupled_traction_moose.fsi)
- [`transient_fuel_cladding_pcmi.fsi`](../verification/fuelsim/transient_fuel_cladding_pcmi.fsi)
- [`transient_finite_strain_pcmi.fsi`](../verification/fuelsim/transient_finite_strain_pcmi.fsi)
- [`steady_finite_follower_pressure.fsi`](../verification/fuelsim/steady_finite_follower_pressure.fsi)
- [`transient_noncoaxial_finite_strain.fsi`](../verification/fuelsim/transient_noncoaxial_finite_strain.fsi)

上述十六张卡分别驱动 M0、两套 M1、M5.4、M2.1、M3.1、五套 M2.2、M2.3、
M3.3、M4.1、M4.2 和 M4.3 的
fuelsim-to-MOOSE 对比；测试程序不再直接构造这些案例的材料、载荷路径或
网格选择参数。每个对比读取 MOOSE 最终时刻的全部节点，统一检查温度、
径向位移和轴向位移的三项误差；M1、M2.3、M3.3 和 M4.1 还检查全部接触
节点的压力三项误差。
