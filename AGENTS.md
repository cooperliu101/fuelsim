# fuelsim Agent Guide

本文件适用于整个仓库。

## 沟通与描述方式

面向用户的所有说明、审查结论和进度汇报必须使用完整、平实的句子，不得使用
两字缩略语或只有上下文里才懂的内部简称。具体要求：

- 每个术语第一次出现时解释含义（例如不说"拒步"，而说"拒绝这个时间步、
  缩小步长重试"）。
- 结论与建议分开陈述；未能亲自验证的内容必须说明验证方式和边界。

## 目标与当前范围

`fuelsim` 使用 C++17 开发核燃料性能有限元程序。当前生产入口从一个 Exodus
文件自由组合任意数量的 2D 轴对称 RZ 区域，并选择稳态或瞬态求解：

```text
体单元/界面固定 12 DOF -> ADlite 局部 Jacobian
                       -> PETSc 统一稀疏装配 -> SNES Newton
```

接触使用 secondary-side STS 热接触、secondary 节点到 primary 线段的唯一
NTS 机械接触和 field-major 全局自由度：

```text
[T(:), ur(:), uz(:)]
```

M2.1 在相同空间离散上增加 Backward Euler 一致热容、物理时间步和
committed/trial/commit/rollback。M2.2 增加通用 J2 Norton 蠕变、J2
线性硬化塑性及两者在同一材料点的全隐式耦合；当前不包含真实燃料或包壳
经验模型。每个区域可独立选择 `small` 或 `finite` 应变；有限应变采用与
MOOSE 默认一致的轴对称增量 Taylor 应变、Rashid 转动、历史张量客观旋转和
当前构形力学弱式，已在非匹配网格 PCMI 中与 MOOSE 对比。有限应变
follower pressure 另有独立 MOOSE 对比；非共轴耦合塑性—蠕变路径另以
畸变四单元、100 个时间步和超过 25 度的转动逐步对比 MOOSE，并包含当前
构形压力和分量牵引。

## 依赖与 C++ 约束

- 仅使用 C++17。
- 不定义项目自己的 C++ 类模板、函数模板、表达式模板或标量泛型层。
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
- 不使用 `FetchContent` 或构建时网络下载。
- 类成员变量统一采用 MOOSE 风格的前缀下划线，如 `_parameters`；不得使用
  `parameters_` 后缀。

## 数值契约

- 局部自由度顺序固定为
  `[T0..T3, ur0..ur3, uz0..uz3]`。
- ADlite 只按 12 个体单元或界面局部自由度播种，禁止按全局自由度播种。
- RZ 积分测度为完整的 `2*pi*r*detJ*w`。
- 应变和应力分量顺序为 `[rr, zz, hoop, rz]`，`rz` 是张量剪应变。
- 小应变区域在参考构形装配力学；有限应变区域从轴对称变形梯度形成
  `Fhat=F_new*inverse(F_old)`，使用 MOOSE 默认 Taylor 应变增量和 Rashid
  增量转动，并用 Cauchy 应力、当前构形形函数梯度和
  `2*pi*r_current*detJ_current*w` 装配内力。
- 有限应变试探态必须保持面内变形 Jacobian、`F_hoop` 和当前半径为正；
  非法态必须作为 domain error 进入线搜索或拒步，不得夹持。
- 热传导和 M2 热容继续在参考构形积分。有限应变区域的 pressure 必须使用
  当前半径、当前法向和当前表面测度形成 follower load；小应变 pressure 与
  所有分量 traction 使用参考构形。
- 有限应变材料必须在中间构形更新，并在步末以增量转动客观旋转应力以及
  弹性、塑性和蠕变张量历史；等效塑性和等效蠕变标量不得旋转。
- MOOSE 默认 `ADComputeMultipleInelasticStress` 只客观旋转应力、弹性应变和
  `combined_inelastic_strain`；其具体模型暴露的 `plastic_strain` 与
  `creep_strain` 不旋转，不能作为 fuelsim 分机制客观张量的逐分量参考。
  非共轴对标必须比较 MOOSE 的总非弹性张量、应力、弹性张量和两个等效标量；
  fuelsim 分机制张量另由局部客观性测试约束。
- 历史变量使用 `double` 保存；只有 trial state 使用 ADlite。
- M2 热容使用参考构形一致质量矩阵：
  `N_i*rho*cp*(T_new-T_old)/dt`，不包含位移惯性。
- M2 的所有 Newton、线搜索和失败重试必须从同一 committed 积分点状态
  重算 trial；只能在最终收敛解上重算一次并提交。
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
- 接触间隙统一使用当前轴对称 RZ 几何：`g=(x_primary-x_secondary)·n_current`，
  其中 `n_current` 是由当前 primary 线段切向构造且从 secondary 指向 primary
  的单位法向；开放为正、穿透为负。圆柱侧面、水平端面和斜面不得分设不同的
  间隙、法向、切向或轴对称面积代码路径；竖直圆柱面的
  `g=(Rp+urp)-(Rs+urs)` 只能作为通用公式的退化结果。
- 接触允许零参考间隙（初始贴合）。secondary 节点恰好骑在 primary 线段
  上时，法向朝向由单元材料侧拓扑确定：以 secondary 边父单元质心相对
  primary 线段沿基准法向 `(tangent_z, -tangent_r)/length` 的有符号距离
  为提示，取法向背离 secondary 材料、指向 primary 一侧，使节点向
  secondary 材料内部移动时间隙增大。该约定与把同一几何的间隙打开无穷
  小量后按间隙符号得到的朝向完全一致。两侧父单元质心位于线段同侧
  （材料重叠）或质心恰好落在线上（退化单元）的畸形几何仍必须明确
  报错，不得隐式夹持。
- 气隙导热为 `h=k_gap/max(g,g_min)`。
- 法向压力为 `p=penalty*max(-g,0)`，罚参数单位为 `Pa/m`。
- 当启用 Coulomb 摩擦时，粘着切向牵引使用与法向相同的接触罚刚度；滑移
  牵引上限为 `mu*p`，切向方向由当前构形投影确定。摩擦粘滑状态必须进入
  committed/trial/commit/rollback 事务和检查点；`mu=0` 必须保持无摩擦路径
  的逐位结果。
- 法向增广拉格朗日接触使用非负法向乘子和互补更新；自动罚刚度按两侧法向
  柔度串联及界面网格尺度计算，显式输入优先于自动值。
- 热接触在构造期按参考 secondary-to-primary STS 重叠分片生成 secondary 侧
  积分点，并为每个积分点与完整 primary 开放链的每条线段预留固定 12 自由度
  稀疏候选。每次状态验证必须在当前构形为每个积分点选择唯一有效 primary 段；
  内部顶点使用半开区间归属，只有整条链的首端和末端可以保留所属端点。积分点
  滑出完整 primary 链时必须明确拒绝当前状态，不得夹持到链端、使用负形函数
  外插或继续使用陈旧候选。只装配唯一活动候选，在当前 primary 法向上计算通用
  RZ 有符号间隙，在 secondary 当前轴对称表面测度上积分，并将严格相反的热流
  投影到 primary 节点。
- 机械接触采用唯一 NTS 投影；secondary 节点反力按当前半边面积集总，并按
  primary 线段形函数分配相反反力。
- 机械 NTS 必须在当前构形上计算轴向或一般法向投影；每次状态验证都必须按
  当前几何重建完整 primary 链的候选段，并保持每个 secondary 节点至多一个
  有效 primary 线段。首次装配必须为所有潜在候选无条件预留稀疏零块，不得退回
  jax_fuel 当前用于特定 MOOSE 对标的参考构形固定机械投影。
- 热接触积分点和机械接触节点都采用半开区间确定内部 primary 顶点的唯一归属；
  只有整条 primary 链的首端和末端可以保留所属端点。不得让内部相邻线段重复
  装配，也不得让陈旧候选伪装成有效投影。机械接触可保留参考链首尾节点的
  物理端点支承；该机械专用端点夹持不得用于热接触积分点。
- 已激活机械接触中任一 secondary 节点若从其全部候选线段失去投影，必须通过
  `validate_state` 明确拒绝当前 Newton 状态并进入线搜索或拒步恢复，不得静默
  置零接触力或继续使用端点力。
- 热接触与机械接触都必须离散守恒。
- 有限应变 pressure 使用当前法向与当前表面测度；traction 默认使用参考构形，
  `configuration = current` 时方向仍固定为全局 R/Z 分量，但当前周长和边长
  必须进入 ADlite 几何切线。
- 一个 M1 载荷路径只能构造一次问题几何，并在所有载荷步复用同一组
  SNES、Vec、Mat、非零结构和回调缓冲区；载荷步只更新具体热源参数。
- 内部计时使用单调时钟，至少区分问题构造、求解器设置、非线性求解、残量
  回调和 Jacobian 回调，并记录回调次数及 PETSc 工作区构造次数。

## 架构边界

- `fuelsim_core`：网格、自由度、材料、Quad4 RZ 核和问题定义，仅依赖
  ADlite。
- `fuelsim_io`：严格解析带版本号的 `.fsi` 输入卡，并使用 Exodus API 在
  `.e` 文件和 fuelsim 自有非结构 Quad4 网格及结果之间转换；保留元素块、
  节点集和边集的 ID 与名称，不使用 DMPlex，不暴露 Exodus 类型，也不实现
  对象工厂、表达式求值或兼容别名。
- `fuelsim_solver`：PETSc 会话、稀疏装配、SNES 求解、稳态加载和瞬态时间
  推进。
- `NonlinearProblem` 只作为求解器端口；不得扩张成 MOOSE 式对象工厂。
- 所有区域节点必须保持独立；默认 PCMI 包壳高度比芯块高 `20 um`，界面
  通过轴向投影耦合。
- 不复制 MOOSE 的对象工厂、继承层次或输入参数系统。
- 不复制 jax_fuel 的运行时声明式 Kernel 注册系统。
- 新物理先形成具体、可验证的局部残量，再考虑通用化。
- 生产问题类型只保留 `SteadyProblem` 和 `TransientProblem`；M0/M1/M2
  只作为路线与回归名称。旧的专用问题类只能留在 `tests/support` 中支撑
  已有回归，不得重新进入公共头文件或生产库。
- `SteadyProblem` 和 `TransientProblem` 从一个 `UnstructuredQuad4Mesh`
  选择任意数量的命名块；每个块独立建立区域自由度、材料和
  `small|finite` 应变形式。
- Contact 输入只接受 `primary` 和 `secondary` 边集名，不接受主/从 block；
  所属区域必须由 Exodus 边集相邻单元解析。每个接触对可独立启用热接触、
  机械接触或两者。
- 瞬态问题使用 `TransientProblem`、自由函数 `solve_transient` 和
  `Quad4RzTransientKernel`；不得把时间状态职责塞入 PETSc 回调。
- 不增加材料工厂或标量泛型层；真实燃料/包壳关联应在通用状态事务稳定后
  作为单独里程碑实现。
- 除非用户明确要求，不增加旧 API 别名、适配器或兼容层。
- 用户运行入口固定为 `fuelsim -i <case.fsi>`。输入 v1 只接受一个 Exodus
  文件，使用 SI 单位和严格字段集合，不提供 include、宏、表达式、单位换算、
  旧键别名或隐式默认问题；网格几何与离散规模必须来自 Exodus 文件。

## 必须执行的验收

后续开发使用原 `moose` Conda PETSc，并直接链接独立的串行 Exodus I/O
库。PETSc 不需要启用 Exodus；不得使用 DMPlex 或 PETSc Exodus viewer。
先安装 ADlite 和 Exodus，然后配置 fuelsim：

```bash
fuelsim_dependency_root="$(cd .. && pwd)/fuelsim-dependencies"

env \
  PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  PKG_CONFIG_PATH=/home/cooper/miniforge/envs/moose/lib/pkgconfig \
  cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH="${fuelsim_dependency_root}/adlite-a3778d2" \
  -DSEACASExodus_DIR="${fuelsim_dependency_root}/exodus-2024-06-27/lib/cmake/SEACASExodus" \
  -DFUELSIM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

涉及 CMake、PETSc 求解层或依赖配置的修改必须完成上述入口的完整回归。

PETSc/MPICH 测试在受限沙盒内可能出现 `OFI EP enable failed`。遇到该错误应在
沙盒外重跑，不能归因于 fuelsim 数值实现。

每次物理修改至少检查：

1. 对应局部 AD Jacobian 与中心差分方向导数；
2. 相关解析解；
3. 默认端到端求解；
4. 匹配物理、罚参数、加载路径和网格设置的 MOOSE 对标量。

接触搜索或投影修改还必须检查：内部 primary 顶点参考态只有一个所有者，
secondary 节点或热接触积分点滑过该顶点后所有权唯一转移且不双计，以及滑出
完整 primary 链时失投影守卫能够拒绝该状态；参考态竖直但当前态倾斜的侧面、
水平端面和一般斜面都必须覆盖，并证明它们走同一通用 RZ 公式。热接触还必须
检查同一参考分片的两个积分点可分别选择不同 primary 段、只有唯一候选计热且
两侧热残量严格守恒。

所有 fuelsim-to-MOOSE 对比必须读取 `verification/moose/` 下由对应 MOOSE
输入生成并追踪的 `*_mesh.e`，不得在对比测试内使用 `make_annulus` 或硬编码
坐标重建等价网格。后续新增或修改的场量误差对比统一采用相对 L2 误差、
相对绝对峰值误差和最大逐点相对误差，三项都必须满足该算例的验收门槛；
相对最大范数不再作为默认对比指标。最大逐点相对误差只对非零参考值定义；
零参考点的数量和最大绝对差必须单独输出，不得静默增加分母下限。
误差门槛不得在根因未明时放宽；任何 `verified -> qualified` 或 qualified
例外扩大都必须先记录物理/离散/参考量口径原因、实际三项误差和适用边界。

纯解析、局部 Jacobian 和事务测试不受此限制。M0、M2.1
和 M2.2 使用单区域块 0；M1、M2.3 使用命名的 `fuel`/`clad` 块。当前可从
一般非结构 Quad4 文件读取元数据，生产求解网格必须保留每个选中块的原始
节点坐标和 Quad4 连接关系，不得重建张量积网格。M1 必须同时运行
`fuelsim_m1_exodus_moose_tests` 和 `fuelsim_m1_unstructured_moose_tests`；
后者逐节点比较温度、径向位移和轴向位移，并逐接触节点比较压力；全部三项误差
指标必须小于 `1%`。通用当前法向和当前热投影下，结构化和非张量 M1 的轴向
位移最大逐点相对误差实测为 `0.76044%` 和 `0.76783%`，绝对差 `1.092 nm` 和
`1.095 nm`，均发生在约 `0.143 um` 的非零参考值处。MOOSE 的双侧 quadrature
`GapHeatTransfer` 在该结构网格最终步两侧热率相差 `0.142038 W`，约占
`106.7 W` 界面热率的 `0.133%`；fuelsim 仍必须保持严格离散守恒，不得通过复制
这一参考不平衡来降低场误差，也不得增加逐点相对误差的分母下限。

M2 还必须检查：

1. residual/Jacobian 重复调用不修改 committed 历史；
2. 接受步只提交一次，拒绝步完整回滚；
3. 活跃蠕变、塑性和耦合分支的 AD 切线分别通过中心差分；
4. 多时间步只构造一次 PETSc 工作区；
5. MOOSE 瞬态温度、应力、位移、等效塑性应变和等效蠕变应变均小于
   `0.1%`；M2.3 的节点场、接触压力、总力、平均状态和 40 个积分点三指标
   也必须统一小于 `0.1%`。

有限应变修改还必须检查：

1. 均匀轴对称伸长的 Taylor 应变增量、非零 committed 构形和当前体积测度
   解析解；
2. 有限应变弹性及活跃塑性—蠕变分支的局部 AD 切线；
3. 分别隔离非正平面 Jacobian 与非正当前半径的 domain-error 路径；
4. `fuelsim_m41_finite_strain_pcmi_tests` 的非匹配网格 PCMI 全场对比；
5. 温度、径向/轴向位移、接触压力、等效应力、等效塑性应变和等效蠕变
   应变的相对 L2、相对绝对峰值和最大逐点相对误差均小于 `0.5%`。
6. 弹性、塑性和蠕变张量的独立客观旋转测试，以及 follower pressure 的
   局部 AD Jacobian、四类边界父单元法向、当前构形合力和
   `fuelsim_m42_follower_pressure_moose_tests`；端到端必须至少覆盖 `left`、
   `top` 和 `right` 压力。
7. `fuelsim_m43_noncoaxial_finite_strain_moose_tests` 必须覆盖畸变四单元上先
   拉伸、再剪切、轴向反向、剪切反向的 100 步路径，转动超过 25 度，并同时
   施加当前构形 pressure 与分量 traction。积分点历史必须先按每个单元的
   四个积分点参考 RZ 体积平均，再按时间和 Exodus 单元 ID 逐单元比较，禁止
   用跨单元平均代替。节点场、等效历史及张量 L2/峰值误差均小于 `0.5%`；
   倾斜压力面的 MOOSE 参考必须同时向径向和轴向方程施加完整当前法向，不能
   用单个 `ADPressure` 分量代替。换向低应力点的应力和弹性应变最大逐点相对
   误差按验证矩阵中的显式 `0.6%` qualified 门槛验收，总非弹性应变保持
   `0.5%`；不得增加分母下限，且必须输出最大绝对误差及其时间、单元和分量
   位置。还必须分别运行弹性、塑性、蠕变、耦合材料、pressure/traction 载荷
   隔离和单单元全位移控制材料对标；纯蠕变隔离路径的应力过零点只允许验证
   矩阵记录的 `4%` 逐点 qualified 门槛，其余聚合误差仍保持 `0.5%`。不能在
   未定位误差来源时修改生产代码。纯蠕变隔离还必须读取 MOOSE 全节点、全
   时间步历史，在同一位移路径上重放生产材料事务并检查自由力学残量；共享
   状态的材料三项误差门槛为聚合量 `1e-6%`、最大逐点 `0.001%`，自由残量
   相对局部内力尺度必须小于 `0.002%`。
   全部材料点的塑性和蠕变累计迹漂移必须分别输出并小于 `7e-6` 和 `2e-7`。
8. 有限应变 checkpoint/restart 必须在非零变形且塑性、蠕变剪切历史均活跃
   的 committed 状态保存，续算终态与不间断路径逐分量一致。

多接触组合的收敛测试必须逐接触对断言 `active_contact_nodes > 0`，不能只证明
候选面可投影或开放间隙下能够收敛。

当前 M4.1 使用 MOOSE 默认 Taylor 分解和开启的有限应变历史旋转；M4.3 已
鉴定其规定的畸变四单元、最大约 25.5 度非共轴多步路径，但不得外推为任意
转角、任意路径或任意网格的一般大转动鉴定。
稳态有限应变只从参考构形对当前总变形做一次 Taylor 更新，载荷步不累计
材料历史；不得将稳态非共轴加载路径声称为 MOOSE 增量材料路径等价。

性能修改还必须检查：

1. 默认 1,584 DOF 算例与修改前提交的同机配对计时；
2. 固定 CPU 且 OMP/OpenBLAS/MKL/NUMEXPR 均为 1 线程；
3. `benchmarks/` 中 23,010 DOF、20 步算例至少完成一次；
4. MOOSE 使用相同网格、物理、载荷步、直接求解器并关闭文件输出。
5. 并行修改必须通过 1-rank/2-rank 逐自由度等价性和 2-rank 非重复贡献区间
   检查；不得仅以多进程能够启动作为并行验收。

构建成功不等于数值验收通过。只有相关 CTest、解析指标和 MOOSE 指标全部
满足门槛后，才能声称功能完成。
