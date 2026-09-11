# fuelsim Agent Guide

本文件适用于整个仓库。

## 沟通与描述方式

面向用户的所有说明、审查结论和进度汇报必须使用完整、平实的句子，不得使用
两字缩略语或只有上下文里才懂的内部简称。具体要求：

- 每个术语第一次出现时解释含义（例如不说"拒步"，而说"拒绝这个时间步、
  缩小步长重试"）。
- 结论与建议分开陈述；未能亲自验证的内容必须说明验证方式和边界。

## 目标与当前范围

`fuelsim` 使用 C++17 开发核燃料性能有限元程序。当前生产入口读取一个
Exodus 文件，显式选择二维轴对称 RZ 或三维 Cartesian 几何，使用相应网格
中的任意数量命名区域，并选择稳态或瞬态计算。

三维 HEX8（c3d8t或c3d8rt）接触支持 node-to-surface（NTS，节点到面）和
surface-to-surface（STS，面到面），C3D20T 机械接触只启用 STS。C3D20T
的旧 NTS 实现暂时保留在源码中，但问题构造会明确拒绝该路径。全局自由度均采用
field-major 排列：

```text
RZ:       [T(:), ur(:), uz(:)]
Cartesian:[T(:), ux(:), uy(:), uz(:)]
```

## 依赖与 C++ 约束

- 仅使用 C++17。
- 不定义项目自己的 C++ 类模板、表达式模板或标量泛型层。函数模板仅允许以下私有数学实现：
  `elements/src/c3d_common.hpp` 中的 `hughes_winget_rotation`；
  `elements/src/c3d_common.cpp` 中的 `determinant_impl`、`inverse_impl`、`multiply_impl`、
  `reduced_hex8_thermal_hourglass_coefficients_impl`；
  `elements/src/material.cpp` 中的 `rotate_cartesian_tensor_impl`。
  标量类型仅允许 `double` 和 `adlite::Scalar`；矩阵运算固定为 3×3，
  张量旋转固定为六分量三维对称张量，不得扩展为通用矩阵库、材料模型模板、
  单元泛型接口或自动微分播种层。材料张量旋转实现不得依赖型号共用头文件。
  模板内的标量数学函数采用 `using std::函数名` 加非限定调用，通过参数相关查找
  让 `adlite::Scalar` 使用 ADlite 重载、`double` 使用标准库，不得剥离导数来调用数学函数。
  轴对称 Hughes-Winget 计算使用 `cax_common` 中的普通具体类型函数。
- 允许使用 `std::vector`、`std::array` 等标准库模板。
- 自动微分只能使用用户的 ADlite 软件包和具体类型
  `adlite::Scalar`。
- 除 ADlite 外，唯一允许的外部数值依赖为 PETSc；Exodus 仅作为直接网格
  I/O 依赖。
- MPI、BLAS、LAPACK、MUMPS、Hypre 等只能作为 PETSc 的传递依赖，不得由
  fuelsim 单独发现或链接。
- 不引入 Eigen、Boost、fmt、JSON/YAML、CLI、日志或第三方测试框架。
- 测试使用普通 C++ 自检程序和 CTest。
- 公共头文件不得暴露 PETSc 类型；PETSc 保持在 `fuelsim_solver` 实现层。
- 类成员变量统一采用 MOOSE 风格的前缀下划线，如 `_parameters`；不得使用
  `parameters_` 后缀。

## 数值契约

- 局部自由度顺序按单元拓扑固定：RZ Quad4 为
  `[T0..T3, ur0..ur3, uz0..uz3]`，三维 HEX8 为
  `[T0..T7, ux0..ux7, uy0..uy7, uz0..uz7]`，混合阶 HEX20 为
  `[T0..T7, ux0..ux19, uy0..uy19, uz0..uz19]`。
- ADlite 只按体单元或界面局部量播种，禁止按全局自由度播种。RZ 运动学链按
  4 个面内位移梯度分量、积分点径向
  位移和积分点温度（宽度 6）播种，本构关系按 4 个应变分量加温度（宽度 5）
  播种并经 `adlite::compose` 挂回运动学链，单元 12×12 Jacobian 再由参考形
  函数梯度与形函数的线性闭式链组装。HEX8 运动学链按 9 个位移梯度分量加
  积分点温度（宽度 10）播种，本构关系按 6 个应变分量加温度（宽度 7）播种，
  同样经 `adlite::compose` 挂回并由闭式链组装 32×32 Jacobian；HEX20 使用
  相同的宽度 10 运动学链和宽度 7 本构链，并闭式组装 68×68 Jacobian。
- RZ 积分测度为完整的 `2*pi*r*detJ*w`。
- 每个区域必须显式设置 `element`；轴对称只接受 `cax4t|cax4rt|cax8t|cax8rt`。
  默认 `quad4` 型号及其专用实现已经删除，不提供别名或隐式回退。
  `cax4t` 使用四个材料积分点和 Abaqus 轴对称选择性体积处理：平均完整体积应变，
  环向分量单独平均，剩余修正只平均分配到两个面内正应变，不能套用三维的三等分规则。
  有限应变使用增量中间构形的 Hughes-Winget 应变与转动；环向变形梯度按参考体积
  平均，完整应变迹按中间构形体积平均。对应的面内平均应力、环向应力、当前体积
  及跨积分点导数必须同时进入内力和 Jacobian。仍只使用宽度 6 的运动学和宽度 5
  的本构播种，经闭式节点链装配，不得扩大为完整单元自由度播种。
  热膨胀使用四个角点温度的算术平均值；其他力学材料物性仍使用积分点温度。
  平均温度对四个温度自由度的闭式导数以及有限应变下新旧热应变之差必须保留，
  不得把这一规则推广为所有材料函数都使用平均温度。
  `cax4t` 热传导使用 2×2 Gauss 梯度，但导热系数在对应角点温度处求值，
  必须保留该角点温度的独立导数。热容采用 `integral(N_i*dV)` 节点对角权重，
  密度、比热及温度变化率在节点求值。小应变热算子使用参考构形；有限应变
  热传导使用增量中间构形梯度，积分权重为参考积分权重乘单元整体当前/参考
  体积比，热容与体热源使用逐积分点当前构形测度。几何导数和整体体积的跨点
  导数必须进入 Jacobian。上述规则由 B15.0/B15.1 独立热算子探测鉴定。
- `cax4rt` 使用体积平均梯度和一个活跃材料积分点，小应变和有限应变均支持。
  热膨胀使用轴对称体积平均温度 `sum(T_i*integral(N_i*dV))/integral(dV)`，
  不能套用 CAX4T 的四节点算术平均规则。逐节点升温的独立 Abaqus 证据见
  `verification/abaqus/thermal_expansion_probe/README.md`。该探测使用常数膨胀
  系数和零位移约束，不能单独鉴定大变形后的平均构形，也不推广到其他一阶单元。
  有限应变独立平均完整体积增量与环向增量，采用 Hughes-Winget 客观转动。
  力学沙漏能量为 `0.5*C*|Fbar^T*a|^2`，小应变以单位矩阵代替 `Fbar`；
  `a=sum(gamma_i*u_i)`，`gamma_i=Gamma_i-gradbar_i·sum(Gamma_j*X_j)`，
  `Gamma=[1,-1,1,-1]`，`C=0.005*G_initial*V_ref*sum(|gradbar_i|^2)/(sum(gamma_i*Gamma_i))^2`。
  初始剪切模量来自初始温度，不能按验证误差调整。能量、残量及其完整导数必须一致。
  热传导采用均匀梯度与一个沙漏模态；B9.0/B9.1/B9.2/B9.15 原生算子识别得到
  热模态系数 `V*(|gradbar_0|^2+|gradbar_1|^2)/(12*pi)`，保留局部节点顺序。
  热容使用 `integral(N_i*dV)` 作为节点对角权重；体热源同样按
  `integral(N_i*dV)` 分配到节点，不使用中心体积均分。B9.16/B9.17 的独立
  源项开关探测分别验证参考与当前构形分配，源项几何导数和守恒诊断必须一致。
  小应变使用参考构形，有限应变使用当前构形，并闭式保留热残量对位移的几何导数。
  仍只使用宽度 6 的运动学和宽度 5 的本构播种。材料历史沿用容量为四的存储，
  只有第零项活跃，其余项始终为零；输出明确记录 `material_point_count=1`，
  不活跃积分点字段为 NaN，不能作为额外积分点参与统计或更新。
- 应变和应力分量顺序为 `[rr, zz, hoop, rz]`，`rz` 是张量剪应变。
- `cax8t` 和 `cax8rt` 必须使用 QUAD8 网格，局部自由度为
  `[T0..T3, ur0..ur7, uz0..uz7]`，温度只在角点求解。两者分别使用 3×3 和
  2×2 体积分，材料温度在各积分点由角点插值；不使用 CAX4T 的体积处理或
  CAX4RT 的沙漏控制。仍只用宽度 6 的运动学和宽度 5 的本构播种，闭式装配
  20×20 Jacobian。小应变体热算子使用参考构形；有限应变热传导使用中间构形
  梯度和当前二次几何测度，热容使用当前二次几何的一致矩阵，体热源使用温度
  角点形成的线性当前几何测度，完整保留热残量对位移的导数。
  CAX8RT 仅有四个活跃材料积分点，历史更新、时间误差估计、能量诊断和检查点
  只使用这四点；输出 `material_point_count=4`，预留的 `q4`—`q8` 字段必须为
  NaN。以 B12.0—B12.14 的独立 Abaqus CAX8RT 参考按小于 0.1% 验收。
- fuelsim的c3d8t 力学采用与 Abaqus C3D8T 一致的选择性减缩体积积分。小应变时偏应变
  保留八点积分，每个积分点的应变迹替换为按参考体积加权的单元平均迹。平均迹
  对节点位移的闭式链加入 32×32 Jacobian，不得扩大宽度 10 的运动学播种或宽度
  7 的本构播种。有限应变时采用由 B5.19 鉴定的 C3D8T 有限应变选择性体积处理，
  体积平均量、当前积分点体积及其闭式位移链必须同时进入残量和 Jacobian。
  RZ Quad4 和三维 HEX20 不采用这条 HEX8 专用规则。
- 三维 HEX8 选择 `element = c3d8rt` 时采用 Abaqus C3D8RT 一点减缩积分规则，
  不采用上一条 C3D8T 选择性体积积分。体积平均形函数梯度形成一个均匀应变和
  一个体积加权材料温度，每个单元只保存一个材料积分点。热传导由均匀梯度项和
  四个正交沙漏模态组成；均匀体热源使用中心 Jacobian 的 `8*detJ_center` 并平均
  分到八个节点；Backward Euler 热容使用一致热容矩阵行和作为节点对角权重。
  小应变使用参考构形权重，有限应变使用当前构形的体积、平均梯度和节点热容
  权重，并保留热残量对位移的几何 Jacobian。
- C3D8RT 力学沙漏控制固定使用 Abaqus 默认总刚度系数 `0.005` 乘初始温度剪切
  模量，不提供输入覆盖，也不得根据验证误差拟合。有限应变沙漏能量必须把参考
  模态位移通过单元平均变形梯度推前，完整保留 `F*D*transpose(F)` 及变形梯度
  自身导数。C3D8RT 仍只使用宽度 10 的运动学链和宽度 7 的本构链；有限应变
  当前构形和中间构形的体积、平均梯度、热沙漏量及节点热容权重对 24 个位移
  自由度的导数必须闭式计算，本构切线再通过这些几何链闭式组装完整 32×32
  Jacobian。禁止按 24、8 或完整 32 自由度做宽恒等播种，也禁止按 Jacobian 列
  重复完整残量。残量必须使用独立的普通双精度路径，并与 Jacobian 调用返回的
  残量逐项相同。有限应变还必须分别保持 committed、midpoint 和 current 构形
  Jacobian 为正。
- 小应变区域在参考构形装配力学。轴对称有限应变遵循前述各 CAX 型号规则。
- 三维有限应变区域使用 Abaqus 风格的增量应变、Hughes-Winget 客观转动、
  Cauchy 应力、当前构形形函数梯度和当前体积测度装配力学。轴对称试探态必须
  保持面内变形 Jacobian、`F_hoop` 和当前半径为正；三维试探态必须保持
  committed、incremental 和 current 构形 Jacobian 为正，HEX8 有限应变选择性
  体积积分还必须保持 midpoint 构形 Jacobian 为正。非法态必须作为 domain
  error 进入线搜索或拒绝这个时间步、缩小步长重试，不得夹持。
- `cax4t` 使用前述经原生识别的热算子规则。
  `cax4rt` 使用前述减缩积分专用热算子。三维 HEX20
  小应变的上述三个体热算子也使用参考构形。有限应变 C3D20T 按 Abaqus 识别结果
  使用 `3*3*3` 积分：热传导的试函数梯度和温度梯度使用增量中间构形，积分测度
  使用完整二次位移几何的当前构形；一致热容矩阵使用完整二次位移几何的当前构形；
  体热源使用八个温度角点形成的线性当前构形测度，不受位移中间节点影响。上述
  有限应变体热残量必须保留对位移的几何 Jacobian。三维 HEX8 小应变的热传导、
  体热源、角点热容、表面热流和对流使用参考
  构形；有限应变的上述五类热算子均使用当前构形。三维 HEX8 有限应变热残量
  必须保留对位移的几何 Jacobian。
- 有限应变材料必须在中间构形更新，并在步末以增量转动客观旋转应力以及
  弹性、塑性和蠕变张量历史；等效塑性和等效蠕变标量不得旋转。
- MOOSE 默认 `ADComputeMultipleInelasticStress` 只客观旋转应力、弹性应变和
  `combined_inelastic_strain`；其具体模型暴露的 `plastic_strain` 与
  `creep_strain` 不旋转，不能作为 fuelsim 分机制客观张量的逐分量参考。
  非共轴对标必须比较 MOOSE 的总非弹性张量、应力、弹性张量和两个等效标量；
  fuelsim 分机制张量另由局部客观性测试约束。
- 历史变量使用 `double` 保存；只有 trial state 使用 ADlite。
- 小应变三维 HEX20 热容使用参考构形一致质量矩阵；有限应变三维
  HEX20 使用上一条规定的当前构形一致热容矩阵。三维 C3D8T 按 Abaqus
  一阶热单元规则在八个自然坐标角点做节点积分，第 `i` 个节点的热容残量为
  `detJ_i*rho(T_i)*cp(T_i)*(T_i_new-T_i_old)/dt`。小应变的 `detJ_i` 为参考构形
  角点 Jacobian 行列式，有限应变则为当前构形角点 Jacobian 行列式。热容的
  温度—温度块按节点对角集总；有限应变的热残量—位移几何块一般不为零，因此
  不能把整个 32×32 Jacobian 称为对角矩阵。三种体单元都不包含位移惯性。
- Backward Euler step-doubling 必须从同一完整 committed 状态比较节点场、
  应力和全部弹性/塑性/蠕变历史；接受解使用两个半步，拒绝时完整恢复。
- step-doubling 接受两个半步时，热率诊断按两个半步时间平均，功、能量变化
  和耗散按两个半步求和；不得只保留第二个半步的诊断。
- 失败时间步必须同时回滚 nodal state、积分点历史、物理时间和热源；失败
  `SolveResult.state` 不得作为重试初值。
- 轴对称 J2 内积的 `rz` 项必须乘 2；塑性与蠕变应变增量必须无迹。
- Norton 首版为 `rate=A*(q/q_ref)^n` 的等温后向 Euler 更新；J2 首版为
  线性各向同性硬化径向返回。
- 耦合分支必须先求 creep-only 松弛应力；只有该应力仍超过当前屈服应力
  时才激活塑性。活跃时塑性和蠕变共用最终 J2 应力方向并同时满足应力平衡、
  Norton 后向 Euler 方程和塑性一致性条件，禁止一次性算子分裂。
- Norton 与耦合标量根使用对数域保护，等效应力使用 ADlite `hypot`；
  不得通过直接平方和或显式构造超范围幂律系数破坏极端尺度。
- PETSc Dirichlet 约束使用 `F_i=x_i-g_i` 和只清行的
  `MatZeroRows(..., diagonal=1)`。
- 固定残量尺度必须成对提供温度与力学物理尺度，对残量和 Jacobian 同行缩放，
  跨载荷步保持不变，并与基于初始残量的自动缩放互斥。
- PETSc Vec、Mat、SNES 和局部贡献装配支持多个 MPI rank；每个贡献只由一个
  rank 计算，不得退化为每个 rank 重复装配全模型。当前 Exodus 读取、问题
  几何和 committed 状态仍在各 rank 复制；PETSc 回调只收集本 rank 贡献及
  接触搜索依赖所需的影子自由度，输出文件只由 rank 0 写入。
- 不隐式夹持异常材料值或几何值；非法结构输入应明确报错。

- 三维 Cartesian HEX8 接触在当前构形上支持 NTS 和 STS；C3D20T 机械接触只启用
  STS，旧 NTS 实现保留但问题构造必须明确拒绝。每个 secondary 节点或积分约束
  在一次状态验证中至多选择一个有效primary 面候选；内部面边界不得重复归属，首次装配必须为所有潜在候选预留稀疏零块。热流和机械反力在两侧必须严格离散守恒。
- 三维热接触、HEX8 NTS 机械接触、小滑移机械接触和带摩擦的有限滑移机械接触失去
  全部有效投影时，必须由 `validate_state` 拒绝当前 Newton 状态。只有无摩擦
  HEX8 有限滑移 STS 的平均表面约束可以在滑出对面后自然释放，并贡献严格零
  残量；不得把这个例外扩大到其他接触离散。
- pressure 和 traction 省略 `configuration` 时随应变形式采用推荐构形，即
  小应变使用参考构形、有限应变使用当前构形；显式指定 `reference` 或 `current`
  时必须遵从指定构形。pressure 使用所选构形的法向和表面测度；分量 traction
  的方向始终固定为全局分量，但使用所选构形的表面测度。当前构形的周长、边长
  或三维面测度必须进入 ADlite 几何切线。
- 三维 HEX8 C3D8T 和 C3D8RT 的 surface heat flux 与 convection 同样接受
  `configuration = reference|current`；省略时小应变使用参考构形、有限应变
  使用当前构形。当前构形的表面测度及其位移导数必须进入热残量和一致 Jacobian。
- 内部计时使用单调时钟，至少区分问题构造、求解器设置、非线性求解、残量
  回调和 Jacobian 回调，并记录回调次数及 PETSc 工作区构造次数。

## 架构边界

- `fuelsim_elements`：位于 `elements/`，独立配置、构建和测试的局部单元计算库；
  持有 CAX4T/CAX4RT、CAX8T/CAX8RT、C3D8T/C3D8RT、C3D20T/C3D20RT
  体单元、局部边界积分与接触计算，以及材料积分、材料函数、坐标和局部几何，仅依赖 ADlite。
  局部入口接收坐标、材料、节点状态与已接受历史，返回残量、切线矩阵和试探结果；
  体单元按八种型号分别建立同名头文件、实现和测试，不建立 `detail` 目录；
  专用算法留在型号文件。CAX8T 和 C3D20T 分别持有对应二次单元的完整算法，
  CAX8RT 和 C3D20RT 通过全积分型号头文件中的显式积分规则入口调用；只有这两组
  减缩积分型号到全积分型号的调用被允许，其他型号之间不得互相调用。
  其余实际共用实现集中在平铺的 `src/cax_common.*` 和 `src/c3d_common.*`。
  边界、接触按表面拓扑分别命名。
  不持有全局网格、不执行全局装配或状态提交，也不得反向依赖 fuelsim_core。
- `fuelsim_core`：全局网格、自由度、接触候选搜索与归属、全局装配和问题定义，调用
  `fuelsim::elements`；所有时间步历史的接受与失败恢复仍由 fuelsim 负责。
- `fuelsim_io`：严格解析带版本号的 `.fsi` 输入卡，并使用 Exodus API 在
  `.e` 文件和 fuelsim 自有非结构 Quad4、QUAD8、HEX8 或 HEX20 网格及结果之间转换；
  保留元素块、节点集和边集的 ID 与名称，不使用 DMPlex，不暴露 Exodus 类型，
  也不实现对象工厂、表达式求值或兼容别名。
- `fuelsim_solver`：PETSc 会话、稀疏装配、SNES 求解、稳态加载和瞬态时间
  推进。
- `NonlinearProblem` 只作为求解器端口；不得扩张成 MOOSE 式对象工厂。
- 同一 Exodus 源节点被相邻匹配区域共同使用时，必须映射为同一个全局场节点；
  非匹配区域和通过接触耦合的区域必须保留独立源节点。轴对称默认 PCMI 算例的
  包壳高度比芯块高 `20 um`，界面通过轴向投影耦合。
- 不复制 MOOSE 的对象工厂、继承层次或输入参数系统。
- 不复制 jax_fuel 的运行时声明式 Kernel 注册系统。
- 新物理先形成具体、可验证的局部残量，再考虑通用化。
- 生产问题类型只保留 `SteadyProblem` 和 `TransientProblem`；M0/M1/M2
  只作为路线与回归名称。旧的专用问题类只能留在 `tests/support` 中支撑
  已有回归，不得重新进入公共头文件或生产库。
- `SteadyProblem` 和 `TransientProblem` 从与输入几何一致的
  `UnstructuredQuad4Mesh`、`UnstructuredQuad8Mesh`、`UnstructuredHex8Mesh` 或 `UnstructuredHex20Mesh`
  选择任意数量的命名块；一个输入算例只使用一种体单元拓扑。每个块独立建立
  区域自由度、材料和 `small|finite` 应变形式。
- Contact 输入只接受 `primary` 和 `secondary` 边集名，不接受主/从 block；
  所属区域必须由 Exodus 边集相邻单元解析。每个接触对可独立启用热接触、
  机械接触或两者。
- 瞬态问题使用 `TransientProblem` 和自由函数 `solve_transient`，再按几何
  进入 `Quad4RzTransientKernel` 或三维 Cartesian 装配路径；不得把时间状态
  职责塞入 PETSc 回调。
- 不增加材料对象工厂或标量泛型层。允许使用五个类型安全的材料函数注册表，
  分别注册热物性、弹性、本征应变、等效蠕变速率和塑性流动应力函数；注册
  函数必须使用具体 `adlite::Scalar`、严格具名参数和现有统一状态事务。
- 除非用户明确要求，不增加旧 API 别名、适配器或兼容层。
- 用户运行入口固定为 `fuelsim -i <case.fsi>`。输入 v3 只接受一个 Exodus
  文件，并要求显式选择 `axisymmetric_rz` 或 `cartesian_3d` 几何；使用 SI
  单位和严格字段集合，不提供 include、宏、表达式、单位换算、旧键别名或隐式
  默认问题。材料在 `[Materials]` 中由已注册函数组合，区域只用 `material`
  引用；网格几何与离散规模必须来自 Exodus 文件。

## 测试入口与输入卡约束

- 所有端到端测试例题（与外部程序moose或abaqus或基准解对比的例题）必须运行实际生产程序
  `fuelsim -i <case.fsi>`。
- 不同端到端测试例题必须由输入卡区分。每张测试输入卡都必须作为完整、可直接
  阅读和手工运行的 `.fsi` 文件显式提交；材料、区域、边界条件、接触、执行器、
  求解器和输出设置都写在输入卡中。不得由 C++ 测试代码、CMake 脚本、Shell
  脚本、模板替换、宏或运行时文本拼接来生成或修改问题定义。CTest 可以把完整
  输入卡逐字复制到隔离的构建目录，但复制前后的内容必须完全一致。
- 端到端结果检查程序或脚本只允许读取生产程序的返回码、标准输出、Exodus
  结果、工程历史和检查点，并与解析解或受追踪的外部求解器参考比较；不得链接
  求解器、重新求解问题或读取仅存在于生产库内部的状态。参考结果、误差门槛和
  测试调度元数据不得写入生产输入卡。
- 局部自动微分 Jacobian、中心差分方向导数、残量重复调用、历史事务、接触候选
  唯一归属以及无法从生产结果独立证明的内部数值契约，继续由构建期内部测试
  覆盖。这些测试不是用户测试例题，不用端到端终态结果替代。

## 必须执行的验收

开发使用原 `moose` Conda PETSc，并直接链接独立的串行 Exodus I/O
库。PETSc 不需要启用 Exodus；不得使用 DMPlex 或 PETSc Exodus viewer。
所有并行编译命令最多使用 4 个作业；不得使用没有显式作业数的 `--parallel`，
以免 Release 链接时优化同时启动过多链接进程并耗尽 WSL 内存。
所有 CTest 运行统一使用最多 4 个并发测试，不得使用串行完整回归；命令中必须
显式指定 `-j4` 或更小的并发数。
先安装 ADlite 和 Exodus，然后配置 fuelsim：

```bash
fuelsim_dependency_root="$(cd .. && pwd)/fuelsim-dependencies"

env \
  PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  PKG_CONFIG_PATH=/home/cooper/miniforge/envs/moose/lib/pkgconfig \
  cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH="${fuelsim_dependency_root}/adlite-0.2.3" \
  -DSEACASExodus_DIR="${fuelsim_dependency_root}/exodus-2024-06-27/lib/cmake/SEACASExodus" \
  -DFUELSIM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel 4

# 统一回归：所有已注册测试都属于同一个轻量测试套件
ctest --test-dir build -j4 --output-on-failure
```

CTest 不再区分日常回归和完整回归，也不使用 `qualification` 标签筛选扩展参数
扫描。网格、时间步、材料参数和接触参数的重复扫描保留为手动验证资料；自动测试
只保留能够增加独立物理分支、离散路径、事务行为或外部求解器对标覆盖的轻量例题。

涉及 CMake、PETSc 求解层或依赖配置的修改必须完成上述统一回归入口。

PETSc/MPICH 测试在受限沙盒内可能出现 `OFI EP enable failed`。遇到该错误应在
沙盒外重跑，不能归因于 fuelsim 数值实现。

每次物理修改至少检查：

1. 对应局部 AD Jacobian 与中心差分方向导数；
2. 相关解析解；
3. 三维和轴对称均以 Abaqus 为首要参考。轴对称采用 CAX4T 四节点温度—位移
   耦合单元，B7.0—B7.9 覆盖小应变、有限应变、热机械接触、蠕变、塑性以及
   蠕变与塑性同时作用。所有已定义的非零场相对 L2、相对绝对峰值和最大逐点
   相对误差必须小于 0.5%；零参考量单独检查绝对误差，不设置相对误差分母下限。
   具体指标和适用边界见 `verification/abaqus/b7_rz_validation.md`。
   B7.2/B9.5/B10.6/B12.6 按用户要求覆盖屈服前蠕变与随后塑性、蠕变同时
   激活。Fuelsim 全程后向欧拉，Abaqus 屈服前显式、塑性激活后隐式，各自
   与对应分段离散解析解按小于 0.1% 验收，且逐步检查两种机制同时产生增量。
   跨程序位移与蠕变差异仅作明确标注的诊断，不声称其彼此达标。
   该范围只适用于这四个应力上升与保持例题，不推广到其他验收。
   原有 MOOSE 特殊加载、非共轴历史和内部数值契约作为补充覆盖保留。

性能修改还必须检查：

1. 参考程序如moose或abaqus使用相同网格、物理、载荷步、直接求解器并关闭文件输出。
2. 并行修改必须通过 1-rank/2-rank 逐自由度等价性和 2-rank 非重复贡献区间
   检查；不得仅以多进程能够启动作为并行验收。

构建成功不等于数值验收通过。只有相关 CTest、指标全部满足门槛后，才能声称功能完成。
