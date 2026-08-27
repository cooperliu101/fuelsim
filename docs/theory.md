# fuelsim 理论、离散与验证手册

## 1. 文档目的与适用范围

本文档统一说明 fuelsim 当前已经实现的数学约定、有限元离散、材料更新、
接触算法、时间状态事务和 PETSc 求解结构，并给出验证证据的入口。它描述的是
仓库当前行为，不是计划书，也不是核安全软件鉴定声明。

fuelsim 当前求解 2D 轴对称 RZ 热—准静态力学问题。体网格是非结构 Quad4，
边界和界面是 Line2。一个 Exodus 文件可以包含任意数量的命名区域；区域节点
保持独立，区域间只通过明确声明的界面耦合。生产问题只有
`SteadyProblem` 和 `TransientProblem`，用户入口固定为：

```text
fuelsim -i <case.fsi>
```

本手册使用以下证据层次：

- 实现位置说明程序实际执行的离散；
- 局部自动微分与中心差分检查约束一致切线；
- 解析解和制造解约束独立可计算的极限；
- 端到端 CTest 约束输入、网格、装配、求解和输出的组合；
- MOOSE 对标约束相同网格、物理和加载路径下的场量结果；
- `measured` 和 `limitation` 证据分别记录受控观测与尚未鉴定的边界。

机器可读的证据索引是
[`verification/verification_matrix.tsv`](../verification/verification_matrix.tsv)，
证据等级、误差口径和当前边界见
[`docs/verification.md`](verification.md)。MOOSE 输入、网格、结果和校验和的
溯源说明见 [`verification/moose/README.md`](../verification/moose/README.md)
及 [`verification/moose/SHA256SUMS`](../verification/moose/SHA256SUMS)。

## 2. 软件边界与状态流

实现按依赖方向分成三个库：

| 层 | 当前职责 | 直接依赖 |
| --- | --- | --- |
| `fuelsim_core` | 网格、自由度、材料、Quad4 RZ 内核、界面和问题定义 | ADlite |
| `fuelsim_io` | 严格解析版本化 `.fsi` 输入，并直接读写 Exodus 网格、结果、检查点和元数据 | Exodus、`fuelsim_core` |
| `fuelsim_solver` | 分布式向量、稀疏装配、SNES Newton、KSP 线性求解和稳态或瞬态推进 | PETSc、`fuelsim_core` |

输入 v2 只接受一个 Exodus 文件、国际单位制数值和固定字段集合。程序不实现
对象工厂、表达式求值、单位换算、旧键别名或运行时 Kernel 注册。ADlite 之外，
PETSc 是唯一直接数值依赖，Exodus 是唯一直接网格 I/O 依赖。

一次残量和 Jacobian 计算的数据流是：

```text
Exodus Quad4/Line2 几何
  -> 命名区域、边界和接触贡献
  -> 每个贡献固定 12 个局部自由度
  -> ADlite 局部残量与 12 x 12 Jacobian
  -> PETSc 分布式残量和 AIJ 稀疏矩阵
  -> SNES Newton 与 KSP 线性求解
```

`NonlinearProblem` 只暴露自由度数、局部贡献、约束和状态验证等求解器端口。
材料、时间状态和接触历史不会由 PETSc 回调提交。

主要实现位置如下：

| 内容 | 实现位置 |
| --- | --- |
| Quad4、外边界和接触界面的局部残量及 Jacobian | `src/core/rz_kernels.cpp` |
| 热弹性、Norton、J2 及耦合材料更新 | `src/core/material.cpp` |
| 区域布局、边界条件、接触搜索和贡献装配 | `src/core/rz_assembly.cpp` |
| 稳态、瞬态及 committed/trial/commit/rollback | `src/core/spatial_problem.cpp` |
| PETSc 非线性求解器与稳态、瞬态工作流 | `src/solver/solver.cpp`、`src/solver/solve_workflows.cpp` |
| 输入卡解析 | `src/io/case_input.cpp` |
| Exodus 网格和结果、检查点 | `src/io/results_io.cpp`、`src/io/checkpoint.cpp` |
| 案例运行、输出调度和程序入口 | `src/app/main.cpp` |

## 3. 坐标、场量和张量约定

### 3.1 轴对称坐标与自由度

径向坐标记为 `r`，轴向坐标记为 `z`，周向为 `hoop`。轴对称假设表示场量
不随周向角变化，位移只有径向分量 `ur` 和轴向分量 `uz`。温度记为 `T`。

全局自由度采用 field-major 排列，即先存全部温度，再存全部径向位移，最后
存全部轴向位移：

```text
[T(:), ur(:), uz(:)]
```

每个体单元和界面贡献都固定为 12 个局部自由度。Quad4 体单元的顺序是：

```text
[T0, T1, T2, T3,
 ur0, ur1, ur2, ur3,
 uz0, uz1, uz2, uz3]
```

热界面的四个节点按 `[secondary0, secondary1, primary0, primary1]` 排列，
但三个物理场仍保持相同的 field-major 顺序。ADlite 不按全局自由度数播种；
体单元 Jacobian 采用两级窄播种：运动学链按 4 个面内位移梯度分量、积分点
径向位移和积分点温度（宽度 6）播种，本构关系按 4 个应变分量加温度
（宽度 5）播种，经 `adlite::compose` 挂回运动学链后，12×12 单元 Jacobian
由参考形函数梯度与形函数的线性闭式链组装。边界和界面贡献仍按各自局部
自由度恒等播种。

### 3.2 应变、应力和内积

轴对称应变与应力的固定分量顺序是：

```text
[rr, zz, hoop, rz]
```

`rz` 是张量剪应变，不是工程剪应变。因此双重内积和 J2 等效应力使用：

```text
a:b = arr*brr + azz*bzz + ahoop*bhoop + 2*arz*brz
q   = sqrt(3/2 * s:s)
```

其中 `s` 是偏应力，`q` 是 J2 等效应力。塑性和蠕变的流动方向为
`3/2*s/q`，两种非弹性应变增量均保持无迹。积分点历史用 `double` 保存；
只有当前局部试探更新使用 `adlite::Scalar`。

### 3.3 轴对称积分测度

RZ 平面上的面积元旋转一周后形成三维体积，因此参考构形体积分测度是：

```text
dV0 = 2*pi*r0*det(J0)*w_xi*w_eta
```

Quad4 使用 `xi,eta = +/-1/sqrt(3)` 的 2 x 2 Gauss 积分，两个方向的 Gauss
权重均为 1。构造几何时，每个积分点都要求 `det(J0)>0` 且 `r0>0`，非法网格
会被明确拒绝。

Line2 的轴对称表面积测度是：

```text
dA = 2*pi*r*J_line*w
```

构形由具体边界或界面定律决定。热接触和机械接触使用当前 secondary 表面；
pressure 和 traction 通过 `configuration` 选择参考表面或当前表面。省略该字段时，
小应变自动使用参考构形、有限应变自动使用当前构形；显式指定非推荐组合仍可执行并给出提示。
对流使用参考表面。

## 4. Quad4 体单元弱式

### 4.1 插值和小应变运动学

温度和两个位移分量均使用相同的双线性 Quad4 形函数 `N_i`。小应变分量为：

```text
epsilon_rr   = d(ur)/dr
epsilon_zz   = d(uz)/dz
epsilon_hoop = ur/r
epsilon_rz   = 1/2 * [d(ur)/dz + d(uz)/dr]
```

小应变区域的力学残量在参考构形装配。

### 4.2 稳态热方程

fuelsim 的热残量采用传导项减体热源的号约定。对节点 `i`：

```text
R_T_i = integral_V0 [
    k(T) * grad(N_i).grad(T) - N_i*Q
] dV0
```

其中 `k(T)` 是导热率，`Q` 是体积热源。稳态解满足所有自由温度方程的
`R_T_i=0`。

### 4.3 Backward Euler 瞬态热方程

瞬态热方程在同一参考构形上增加一致热容矩阵：

```text
R_T_i = integral_V0 [
    N_i*rho*cp*(T_new-T_old)/dt
  + k(T_new)*grad(N_i).grad(T_new)
  - N_i*Q_new
] dV0
```

`rho*cp` 是体积热容。时间离散是一阶 Backward Euler，位移没有惯性项，
因此力学在每个物理时刻保持准静态。界面没有单独的热容。

### 4.4 小应变力学弱式

内部力残量的径向和轴向分量分别是：

```text
R_ur_i = integral_V0 [
    sigma_rr*d(N_i)/dr
  + sigma_rz*d(N_i)/dz
  + sigma_hoop*N_i/r
] dV0

R_uz_i = integral_V0 [
    sigma_rz*d(N_i)/dr
  + sigma_zz*d(N_i)/dz
] dV0
```

周向应力项 `sigma_hoop*N_i/r` 是轴对称径向平衡中不可省略的部分。边界外载
按“内部力减外力”加入同一残量。

### 4.5 有限应变运动学与当前构形弱式

有限应变区域先由当前总位移形成轴对称变形梯度：

```text
F_rr   = 1 + d(ur)/dr
F_rz   = d(ur)/dz
F_zr   = d(uz)/dr
F_zz   = 1 + d(uz)/dz
F_hoop = 1 + ur/r0
```

瞬态步内以 committed 构形为旧构形，增量变形梯度是：

```text
Fhat = F_new * inverse(F_old)
```

程序对 `Fhat` 使用与当前 MOOSE 默认一致的 Taylor 应变增量和 Rashid 增量
转动。应力和弹性、塑性、蠕变张量历史在步末按该转动客观旋转；等效塑性和
等效蠕变标量不旋转。

有限应变力学残量仍具有上一节的轴对称形式，但使用 Cauchy 应力、当前构形
形函数梯度、当前半径和当前体积测度：

```text
dV = dV0 * det(F_rz) * F_hoop
```

轴对称 RZ 与 HEX20 的热传导、所有体热源和 Backward Euler 热容继续使用参考
构形。有限应变 HEX8 的热传导按照 Abaqus C3D8T 使用当前构形温度梯度和当前
体积测度；它因此具有非零的热残量对位移导数。当前试探态必须满足：

```text
det(F_rz) > 0
F_hoop   > 0
r_current > 0
```

任何一项失败都会作为物理域错误交给 SNES 线搜索或时间步拒绝机制，不会用
夹持值继续计算。

稳态有限应变没有 committed 材料历史。每个稳态载荷延续步都从 `F_old=I`
对当前总变形做一次 Taylor 更新，所以稳态载荷步不构成增量材料路径。需要
非共轴路径历史时必须使用 `TransientProblem`。

## 5. 材料模型与局部更新

材料不是一个固定参数结构，而是由已注册的热物性、弹性、本征应变、蠕变和
塑性函数组合。每个函数拥有严格的具名参数表，活跃输入和输出使用
`adlite::Scalar`。热物性函数返回导热率、密度和比热；弹性函数返回弹性模量
与泊松比；多个本征应变函数的轴对称张量结果相加；蠕变函数返回等效蠕变
速率；塑性函数返回当前流动应力。

局部积分器仍由 fuelsim 统一管理。它从同一个 committed 状态调用纯函数，
执行 J2 关联流动、Backward Euler 更新、塑性—蠕变全隐式耦合、trial 状态
提取和有限应变客观旋转。注册函数不能自行提交历史，也不能访问全局解或
PETSc 对象。

### 5.1 热弹性

各区域独立引用组合材料。内置热弹性函数采用各向同性参数，内置热膨胀
本征应变为：

```text
epsilon_thermal = alpha(T) * (T - reference_temperature)
```

弹性模量、泊松比、热膨胀系数、导热率，以及非弹性参数可选择相对参考温度
的线性斜率。活跃温度参数直接以 `adlite::Scalar` 求值，因此材料参数的温度
链式导数进入局部 Jacobian。线性关系越过物理定义域时程序报告物理域错误，
不隐式夹持。

三维 HEX8 体单元的热膨胀温度遵循 Abaqus 一阶耦合单元口径：每个单元只用
八个节点温度的算术平均值计算本征应变。

```text
T_expansion = (T0+T1+...+T7)/8
```

导热、热容、弹性参数、塑性参数和蠕变参数仍按各积分点插值温度求值。因此该
规则只改变热膨胀本征应变及其温度—力学 Jacobian 链，不把整个材料温度场改成
单元常量。二维轴对称 Quad4 和三维 HEX20 仍保留各自现有的积分点温度热膨胀
离散。

HEX8 小应变力学还采用 Abaqus C3D8T 的选择性减缩体积积分。先从八个积分点形成
参考体积加权的单元平均应变迹：

```text
trace_average = sum_q(w_q*detJ_q*trace(strain_q)) / sum_q(w_q*detJ_q)
```

每个积分点保留原偏应变，只把三个正应变同时增加
`(trace_average-trace(strain_q))/3`。因此偏应变仍在八点求值，而体积响应使用一个
单元平均值。平均迹对节点位移的体积平均形函数梯度链以闭式形式加入 Jacobian，
不扩大运动学宽度 10 或本构宽度 7 的自动微分。

HEX8 瞬态热容同样遵循 Abaqus 一阶热单元的节点积分规则，但它与热膨胀的
单元平均温度是两个独立口径。八个热容积分站位于自然坐标角点，与八个温度
节点一一对应：

```text
R_capacity_i = detJ_corner_i * rho(T_i) * cp(T_i) * (T_i_new-T_i_old)/dt
```

因此常数或温度相关体积热容都只产生对角的温度—温度热容 Jacobian。畸变单元
使用各角点自己的参考构形 `detJ`，不把总单元体积简单平均八份。体热源仍使用
参考构形八点 Gauss 积分；RZ Quad4 与混合阶 HEX20 仍使用一致热容。

HEX8 导热的构形随应变形式变化。小应变使用参考构形；有限应变使用当前构形：

```text
R_conduction_i = integral_Vcurrent k(T)*grad_current(N_i).grad_current(T) dVcurrent
```

当前形函数梯度和当前体积测度都从试探变形梯度计算，所以有限应变热残量对位移
具有非零几何 Jacobian。角点节点热容与体热源目前仍在参考构形积分；对流边界也
仍使用参考表面。这些尚未由 Abaqus 探针识别的项不随本次导热构形修改扩张。

### 5.2 Norton 蠕变

Norton 等温幂律是：

```text
creep_rate = A * (q/q_ref)^n
```

Backward Euler 局部应力方程是：

```text
q + 3*G*dt*A*(q/q_ref)^n - q_trial = 0
```

`n=1` 使用闭式解；一般指数在对数域求有界单调根。等效应力用嵌套
`adlite::hypot` 计算，不直接构造可能上溢的平方和或幂律系数。隐式根导数
通过 `adlite::compose` 接回单元切线。

### 5.3 J2 线性硬化塑性

当前屈服应力为：

```text
yield = yield_stress + H*equivalent_plastic_strain_old
```

当 `q_trial<=yield` 时为弹性卸载或再加载，历史不变。活跃屈服采用闭式径向
返回：

```text
delta_p = (q_trial-yield)/(3*G+H)
q_new   = q_trial - 3*G*delta_p
```

### 5.4 塑性—蠕变全隐式耦合

耦合分支不是一次性算子分裂。程序先求 creep-only 松弛应力 `q_creep`；只有
`q_creep` 仍超过当前屈服应力时才激活塑性。活跃时三个关系同时成立：

```text
q_trial = q + 3*G*(delta_p + delta_c)
delta_c = dt*A*(q/q_ref)^n
q       = yield_old + H*delta_p
```

塑性和蠕变共用最终 J2 方向。所有 trial 状态都从同一 committed 材料历史
重算，Newton 残量、Jacobian 和线搜索不会修改 committed 历史。

这些本构参数用于算法验证，不是二氧化铀、混合氧化物燃料或锆合金的工程
标定关联。

## 6. 边界条件

### 6.1 Dirichlet 约束

给定自由度 `x_i=g_i` 时，PETSc 残量行和 Jacobian 行被替换为：

```text
F_i = x_i - g_i
J_ii = 1
```

实现只清约束行，不清列。约束值可以由严格分段线性时间表驱动；时间推进会
精确命中表节点。

### 6.2 Pressure

pressure 表示大小非负的压缩压力。若 `n` 是父 Quad4 逆时针边界顺序确定的
外法向，则外部牵引是 `-p*n`。因为总残量采用“内部力减外力”，pressure
对残量的贡献为 `+integral(N_i*p*n)dA`。

pressure 在 `configuration = reference` 时使用参考半径、参考法向和参考表面测度；
在 `configuration = current` 时使用当前端点计算当前半径、当前法向和当前 Line2
测度，因而是 follower load；ADlite 同时生成方向、周长和边长变化产生的几何刚度。

### 6.3 分量 traction 与对流

traction 的方向固定为用户选择的全局 R 或 Z 分量。未指定构型时，小应变使用参考表面、
有限应变使用当前表面；显式 `configuration = current` 时使用当前半径和当前边长，但方向仍不随法向
转动。

对流热流是：

```text
q_convection = h_c * (T - T_ambient)
```

它在参考 Line2 表面上积分，正值表示热量从求解区域流向环境。换热系数和
环境温度可以分别由时间表驱动。

## 7. 热接触离散

每个接触对由 `primary` 和 `secondary` 边集定义，所属区域由边集相邻 Quad4
解析。热接触采用 secondary-side segment-to-segment 离散，简称 STS；即在
secondary 线段的积分点计算热流，并投影到 primary 线段。

热接触在构造期按参考投影重叠区间切分 secondary 积分分片并生成积分点。每个
积分点独立编号，并与完整 primary 开放链的每条线段建立固定 12 自由度稀疏候选；
参考分片只定义 secondary 侧积分区域，不冻结 primary 段所有权。每次状态验证
都在当前构形搜索完整 primary 链，为每个积分点选择唯一有效线段。内部顶点采用
半开参数区间 `[0,1)`，只有整条链的最后一段采用 `[0,1]`，因此内部顶点唯一
归下一段。弯折链存在多个有效正交投影时选择绝对间隙最小的线段，平局按
primary 线段序号确定。积分点投影越出整条链时明确拒绝当前状态，不做有限距离
端点夹持。

残量和 Jacobian 计算只装配唯一活动候选，其余预留候选为零。积分点的正交投影
分数、两侧坐标、primary 法向和 secondary 表面测度都在当前轴对称 RZ 几何中
求值。所有边方向统一使用：

```text
g = (x_primary_mapped - x_secondary).n_primary_current
```

其中当前单位法向由当前 primary 线段切向构造并从 secondary 指向 primary。
竖直圆柱面只是该式的退化情况，此时
`g=(R_primary+ur_primary)-(R_secondary+ur_secondary)`；实现没有独立圆柱分支。
`g>0` 表示开放，`g<0` 表示穿透。构造时参考间隙恰好为零的贴合界面与机械
接触共用第 8.1 节的材料侧拓扑法向。默认气隙导热定律是：

```text
h_gap = k_gap / max(g, g_min)
q_gap = h_gap * (T_secondary - T_primary)
```

为对齐 Abaqus 的间隙导热表局部响应，也可以使用线性仿射定律：

```text
T_average = (T_secondary + T_primary) / 2
p = penalty * max(-g, 0)
h_gap = h_reference + h_g*g + h_p*p + h_T*(T_average - T_reference)
q_gap = h_gap * (T_secondary - T_primary)
```

间隙、罚函数压力、两侧平均温度和表面测度都保留在 ADlite 自动微分链中，因而
热残量对温度和位移的耦合块由同一个局部残量得到。非负且有限的 `h_gap` 是物理
域条件，非法试探态会进入线搜索或时间步缩小重试。B5.20 逐列识别了 Abaqus 的
这些局部导数；当前定律只代表表格中的一个线性单元，不包括任意表格的分段插值、
外推和截断语义。

`q_gap>0` 表示热量由 secondary 流向 primary。对 secondary 节点和 primary
节点，残量分别加入：

```text
R_T_secondary += integral_A_secondary N_secondary*q_gap dA
R_T_primary   -= integral_A_secondary N_primary*q_gap dA
```

因此同一积分点的界面热量严格等量反向。非匹配 Line2 网格在构造时仍检查参考
投影覆盖空洞和重叠，但一个 secondary 分片中的不同积分点可以在当前构形选择
不同 primary 段。半开区间保证内部相邻段不会重复计热，完整链候选预留保证大
滑移转移不改变 PETSc 稀疏结构，失投影守卫则阻止陈旧候选或链端夹持继续计热。

## 8. 机械接触离散

### 8.1 唯一节点到线段投影

机械接触采用 secondary 节点到 primary Line2 的唯一投影，简称 NTS。每个
secondary 节点在当前构形选择一条有效 primary 线段。一般斜面上的间隙为：

```text
g = (x_primary_projection - x_secondary).n
```

其中 `n` 是当前 primary 线段法向，从 secondary 指向 primary，仍以 `g>0`
为开放。法向的符号在构造时按参考构形确定：参考间隙非零时取间隙符号；
参考间隙恰好为零（初始贴合）
且节点投影落在 primary 线段上时，`n` 由单元材料侧拓扑确定——取
secondary 边父单元质心相对 primary 线段沿基准法向
`(tangent_z, -tangent_r)/length` 的有符号位置，`n` 背离 secondary 材料
一侧，与把间隙打开无穷小量后的符号约定一致。求值时该符号乘以当前线段
法向，因此参考时竖直的圆柱面变形后倾斜时也会同时进入径向和轴向方程，
不存在圆柱专用分支。两侧父单元质心位于线段
同侧的材料重叠几何在构造时明确报错。压缩法向力使用
secondary 节点的当前半边轴对称面积集总：

```text
F_n = pressure * A_secondary_tributary
```

该节点得到 `+F_n*n` 残量，primary 两节点按投影形函数得到总和为
`-F_n*n` 的残量，因此径向和轴向反力离散守恒。

### 8.2 三维表面到面机械接触

#### 8.2.1 HEX8 小滑移节点中心约束

HEX8 必须显式选择 `discretization = surface_to_surface` 和 `sliding = small`，并且
接触两侧都必须采用小应变。Abaqus B3.8 控制探针表明，一个 C3D8 secondary 四节点面
产生四个节点位置约束。约束开度不是单个 secondary 节点的位移，而是四节点位移经过
下列固定平均矩阵后的值：

```text
A4 = (1/16) * [9 3 1 3; 3 9 3 1; 1 3 9 3; 3 1 3 9]
```

这等价于在 QUAD4 父坐标 `(±0.5,±0.5)` 计算形函数。每个约束面积为面面积的四分
之一；匹配平面上的法向罚切线精确分解为：

```text
K_contact = penalty * transpose(A4) * diag(area/4) * A4
```

对于非匹配 primary 和 secondary 网格，构造期按参考投影所有权递归分割 secondary
面，再用 3×3 积分把上述平均约束的对偶测试函数守恒地投影到所有相交 primary 面。
每个唯一 secondary 节点最终只形成一个节点中心约束，固定其参考 primary 支持和参考
切平面。该过程没有引入接触乘子，也不是经典 mortar 离散。

法向与两个切向分量使用同一个相对位移算子，所以作用力和反作用力严格守恒。Coulomb
摩擦在固定参考切平面内保存三分量弹性切向滑移向量和粘滑状态，并参与提交、回滚及
检查点事务。小滑移不允许接触任一侧采用有限应变。HEX8 省略 `discretization` 时仍采用
原有节点到面路径。

#### 8.2.2 HEX8 有限滑移表面约束

无摩擦且参考 primary 各面的中心法向互相平行时，HEX8 有限滑移使用 Abaqus B4.4
识别出的节点中心平均约束。位移约束位于父坐标 `(±0.5,±0.5)`，翘曲 secondary 面的
平均法向位于父坐标 `(±2/3,±2/3)`；公共边和顶点上的多个有效 primary 投影线性平均，
节点传力采用 C3D8 的 `9/16、3/16、1/16` 恢复矩阵。每次状态验证都按当前构形更新
secondary 面积、法向和 primary 投影系数，同时保留覆盖整个 primary 表面的固定稀疏
支持。约束完全滑出 primary 表面时自然释放，标记为无投影并贡献严格零残量。

启用摩擦或参考 primary 各面并不平行时，仍采用 B4.1 的四个当前构形节点中心积分点。
每个点在完整 primary 四节点面集合中搜索唯一最近投影，并以当前 secondary 面测度积分：

```text
R_secondary_i = sum_q(N_secondary_i(q) * (pressure*n + traction_t)_q * J_secondary(q))
R_primary_j   = -sum_q(N_primary_j(q) * (pressure*n + traction_t)_q * J_secondary(q))
```

两种离散都不是经典 mortar，也不引入接触乘子。局部积分路径在内部 primary 边使用唯一
所有权；接触点越过分片时更新主面，滑出完整主表面时拒绝试探状态。法向、两个切向分量、
最近点参数和当前面积均进入该局部路径的自动微分切线。摩擦历史用当前 primary 面的随动
正交切向基保存；提交的
非零双分量弹性滑移先客观运输到当前切平面，共同刚体转动不产生伪滑移。该路径允许
接触两侧分别使用小应变或有限应变。

#### 8.2.3 HEX20 节点中心约束和积分路径

两侧小应变且选择小滑移时，HEX20 使用八个 Abaqus 对标识别出的二次节点中心约束。
有限滑移或任一侧采用有限应变时，则在当前 secondary 二次面上使用 3×3 高斯积分。
primary 和 secondary 面的数量及分片可以不同；每个积分点独立投影到唯一 primary
二次面，间隙、当前法向、压力和摩擦牵引都在该点计算。积分路径的离散残量为：

```text
R_secondary_i = integral(N_secondary_i * (pressure*n + traction_t) dA)
R_primary_j   = -integral(N_primary_j * (pressure*n + traction_t) dA)
```

两侧形函数都满足分片统一，所以作用反作用严格守恒。线性化直接对积分点投影、当前
法向和当前面积自动微分；对于固定法向的闭合模式，罚切线具有
`penalty * B_gap^T * B_gap * dA` 的非负耦合形式。Q8 角点的等效节点力可以为负，
但它只是上述一致残量的投影结果，不是负的对角罚弹簧。

Coulomb 弹性切向滑移和粘滑状态在积分路径中各保存九份。提交、回滚和检查点都以
积分点为事务单位。输出的节点压力由一致等效节点力除以一致等效节点面积恢复，只用于结果展示。
输入显式选择 `node_to_surface` 时才使用旧的八节点投影和
`quad8_nodal_area_rule`，该路径用于算法研究与传统 MOOSE node-face 对比。

HEX20 的 `sliding = small` 固定参考构形确定的 primary 面；`sliding = finite` 则在每次
状态验证中按当前构形为每个积分点重新搜索完整 primary 面集合。候选投影以最小绝对间隙
选择，共边的相等候选按稳定面序号只保留一个；滑出完整表面会触发几何域错误并进入
线搜索或载荷步缩减。超过 64 个 primary 面时使用当前包围盒空间树，矩阵稀疏结构仍在
构造期为全部潜在面配对预留。

有限滑移的摩擦历史除两个全局切向分量外，还保存已提交的单位法向和第一切向。旧的
弹性滑移按这组基分解，再在当前接触基中重构。曲面滑移增量使用已提交接触点的切向基，
不能使用“当前投影坐标在旧主面上的切向基”，否则圆柱面滑动会错误地按 secondary
半径累计。法向、投影参数、两个切向方向和当前 secondary 面积都由 ADlite 线性化。

### 8.3 罚接触和自动缩放

显式罚形式是：

```text
pressure = max(-penalty*gap, 0)
```

`penalty` 的单位为 `Pa/m`。省略显式值时，两侧边界相邻单元的法向尺度为
参考平面面积除以边长，界面刚度按两侧法向单元刚度串联：

```text
k_interface = 1 / (h_primary/E_primary + h_secondary/E_secondary)
penalty     = penalty_factor * k_interface
```

显式 `penalty` 与 `penalty_factor` 互斥。

### 8.4 增广拉格朗日法向接触

增广拉格朗日形式为每个 secondary 接触节点保存非负法向乘子 `lambda`。
固定乘子的内层 Newton 残量使用：

```text
pressure = max(lambda - penalty*gap, 0)
```

只有内层 Newton 收敛后才更新一次：

```text
lambda_new = max(0, lambda - penalty*gap)
```

外层同时检查穿透容差和互补状态。达到最大增广迭代数仍不满足时，整个载荷步
或时间步失败，节点场、材料历史和接触乘子恢复到同一 committed 状态。

### 8.5 Coulomb 摩擦

当前切向罚刚度与法向罚刚度使用同一个 `penalty=k`。令 `t` 为当前 primary
线段单位切向，`delta_s` 为本步 secondary 位移增量减去投影 primary 位移
增量后在 `t` 上的分量，已提交弹性切向滑移为 `s_old`：

```text
s_trial   = s_old + delta_s
tau_trial = k*s_trial
tau_limit = mu*pressure
```

预测牵引在上限内时粘着：

```text
tau   = tau_trial
s_new = s_trial
```

超过上限时滑移：

```text
tau   = sign(tau_trial)*tau_limit
s_new = tau/k
```

预测值恰好位于上限时，未提交为滑移的节点保持粘着；已提交滑移节点继续沿
同向运动时保持滑移。这一分支规则避免在上限处无条件翻转粘滑状态。

开放接触的切向牵引和弹性滑移为零。切向力使用与法向力相同的 secondary
集总面积，再按 primary 形函数分配等量反向力。稳态延续只在收敛载荷增量后
提交 `s_new` 和粘滑标志；瞬态把它们纳入完整状态事务和检查点。

### 8.6 完整 primary 链动态搜索

问题构造时，每个 secondary 节点与完整 primary 开放边链的所有线段都预留
固定 12 自由度贡献，从第一次 Jacobian 装配起就包含所有潜在稀疏耦合。每次
状态验证按当前构形检查完整链并选择绝对法向距离最小的有效段。

内部 primary 顶点采用半开区间：前一段不拥有第二端点，后一段拥有第一
端点。只有整条开放链的首端和末端保留物理端点支承。等距非相邻候选按有序
段编号确定性选择。这样内部顶点只有一个所有者，不会双计接触力。

当前状态滑出某个局部邻域时可转移到完整链上的其他线段；滑出整条 primary
链时则报告物理域错误。候选选择由 committed 几何和当前状态确定性重建，
检查点不保存临时搜索窗口。热接触对每个积分点采用同样的完整链动态候选原则，
但不使用机械参考端点节点的物理端点支承；热积分点滑出完整链时直接报告物理域
错误。M5.7 已对规定网格、载荷和接触参数下同时激活 Coulomb 摩擦、机械大滑移
与热接触动态所有权的路径完成鉴定；该结果不能外推为任意组合路径。

## 9. 时间积分和状态事务

### 9.1 committed 与 trial 的职责

`TransientProblem` 拥有唯一已提交物理状态，时间控制器和检查点在此基础上
另外保存下一控制步长：

- field-major 节点温度和位移；
- 物理时间和载荷因子；
- 每个区域、单元和四个积分点的应力及弹性、塑性、蠕变历史；
- 接触弹性切向滑移、粘滑标志和法向增广乘子；
- 最后一个接受步的守恒和耗散摘要。

下一控制步长不是材料或节点物理状态，也不由 PETSc 回调修改；检查点将它与
上述 committed 状态一起保存，使重启动后的时间控制器连续。

一次时间步的事务是：

```text
begin_time_step
  -> 所有 residual/Jacobian/line-search 从同一 committed 状态重算 trial
  -> SNES 失败：rollback_time_step
  -> SNES 收敛：在最终解上重算一次全部 trial 历史
  -> commit_time_step
```

失败的 `SolveResult.state` 只用于诊断，不能作为重试初值。回滚同时恢复节点、
时间、载荷、材料历史、接触历史和热源。

### 9.2 时间步调节

时间控制器可以根据非线性迭代窗口增长、保持或缩短下一步。求解失败时按
`cutback_factor` 缩小时间步并从同一 committed 状态重试。经历过 cutback
的接受步之后，下一步不会立即放大回刚失败的尺度；至少再经历一个没有缩步
的成功步后才允许增长。

分段线性时间表的节点是事件时刻。为精确命中事件而截短的步不会被误认为
非线性失败，也不会改变名义控制步长。

### 9.3 Backward Euler step-doubling

时间误差控制由 `time_error_relative_tolerance>0` 显式启用，生产默认关闭。
每个候选步从同一完整 committed 状态分别计算一个全步和两个半步，并比较：

- 温度、径向位移和轴向位移；
- 应力；
- 弹性、塑性和蠕变张量历史；
- 等效塑性和等效蠕变历史；
- 接触弹性滑移和法向乘子。

误差满足门槛时接受两个半步的终态。两个半步的热率按物理时间平均，功、能量
变化和耗散按两个半步求和。误差过大时完整恢复候选步之前的 committed 状态
并缩小时间步。

### 9.4 检查点

检查点只保存 committed 状态，不保存 Newton trial、残量、线搜索或失败尝试。
当前严格格式检查 magic、版本、字节序、长度、校验和和问题签名。格式不兼容
或模型签名不一致会被明确拒绝，不提供跨版本迁移。文件先完整写入临时文件，
成功关闭后再替换目标，避免把半写文件安装为新检查点。

## 10. 自动微分、装配和 PETSc 求解

### 10.1 局部自动微分

每个体、边界或界面贡献执行以下固定流程：

```text
12 个 double 局部值
  -> adlite::seed_identity
  -> 12 个 ADlite 局部残量
  -> adlite::extract_jacobian
  -> 12 x 12 double 局部 Jacobian
```

本构隐式标量根的导数通过 `adlite::compose` 嵌回局部残量。全局方向导数检查
会累加所有局部 Jacobian，按 PETSc 规则替换 Dirichlet 行，再与完整约束残量
的中心差分比较，并分别报告温度、径向位移和轴向位移误差。

### 10.2 稀疏结构与贡献所有权

PETSc 使用分布式向量和 AIJ 稀疏矩阵。全部局部贡献按连续贡献索引区间分给
MPI rank，每个贡献只由一个 rank 计算；PETSc 允许向非本地残量行和矩阵行
累加。第一次 Jacobian 装配后锁定新非零位置，完整 primary 链接触的全部
潜在耦合已在构造期预留。

每个 rank 根据自己贡献的 `contribution_dofs` 建立去重影子自由度集合。
接触搜索还显式加入该 rank 涉及的 secondary 节点在完整 primary 链上的所有
候选依赖。残量和 Jacobian 回调只收集这一集合，不再把完整 Newton 试探向量
复制到每个 rank。求解结束后仍收集一次完整状态，以保持 committed 事务和
根 rank 输出。

当前网格、问题几何、committed 节点状态和材料历史仍在每个 rank 复制，
Exodus 也由每个 rank 使用串行 API 读取。因此当前是分布式贡献装配和线性
代数，不是分布式网格与分布式历史存储。

### 10.3 非线性与线性求解

PETSc SNES 使用 Newton 线搜索。生产默认是 BASIC 全步；一次 BASIC 失败时，
fuelsim 可以从该次调用的原始初值用 BT 回溯线搜索重试。BT 仍失败后，稳态
载荷延续执行载荷增量二分，瞬态执行时间步 cutback。

SNES 报告正收敛原因后，程序仍重新计算总残量和温度、径向位移、轴向位移
三个分场残量。固定物理残量尺度会对残量和 Jacobian 同行缩放，并在所有
载荷步保持不变；它与基于初始残量的自动缩放互斥。

线性求解支持直接 LU/MUMPS，以及 GMRES 配合块 Jacobi、温度—力学乘法场
分裂或 HYPRE。默认工程路径仍使用直接求解。迭代组合的当前工程规模结果是
受控测量，不构成一般扩展性结论，具体数据见
[`docs/m5.md`](m5.md) 和 [`benchmarks/README.md`](../benchmarks/README.md)。

一个稳态载荷路径和一个瞬态时间路径都只构造一次问题几何及 PETSc SNES、
Vec、Mat、非零结构和回调缓冲区。载荷步或时间步只更新具体参数和状态。

## 11. 守恒、功和输出诊断

每个接受步记录以下独立诊断：

- 体热源生成率、单元对应离散的热容储热率、边界热率和全局热平衡；
- 每个热接触对在 secondary 与 primary 两侧的热率及不平衡；
- 内力功、外载功、约束反力功、接触功和机械功平衡；
- 弹性能变化、塑性耗散和蠕变耗散；
- SNES 与 KSP 迭代数、残量和回调次数；
- PETSc 工作区构造次数和分布式影子自由度规模。

只有 rank 0 写控制台、CSV、Exodus 和检查点。Exodus 结果按源节点和源单元
标识映射回输入网格；未选实体和无法定义的接触量写为 `NaN`，不伪造成零。

## 12. 验证方法和证据索引

### 12.1 通用误差口径

fuelsim-to-MOOSE 场量比较统一报告：

```text
relative_L2          = ||x-x_ref||_2 / ||x_ref||_2
relative_abs_peak    = |max(|x|)-max(|x_ref|)| / max(|x_ref|)
max_pointwise_relative = max_i |x_i-x_ref_i|/|x_ref_i|
```

最大逐点相对误差只在非零参考值上定义。参考值为零的点单独报告数量和最大
绝对差，不增加分母下限。门槛、qualified 例外和实测值以验证矩阵及对应 CTest
为准。

### 12.2 自动证据入口

以下索引按理论组成指向当前验证能力。完整 CTest 名称、参考文件和门槛以
验证矩阵为唯一机器可读来源。

| 理论或工程组成 | 验证矩阵标识 | 主要证据 |
| --- | --- | --- |
| 严格输入和问题构造 | `input.v3`、`io.exodus` | 输入拒绝测试、具名材料参数、Exodus 元数据回读、严格重启动 |
| 稳态 RZ 体弱式 | `m0.steady` | 实心圆柱温度、自由热膨胀、厚壁圆筒和 MOOSE 全场 |
| 无摩擦热—力接触 | `m1.contact`、`m33.contact` | 非匹配 STS/NTS、斜面、端面、多区域和 MOOSE 全场 |
| HEX8 表面到面接触 | `b38.hex8_sts_identification`、`b39.hex8_sts_multicase`、`b40.hex8_sts_friction`、`b41.hex8_sts_finite_sliding`、`b42.hex8_sts_friction_objectivity`、`b43.hex8_sts_finite_strain` | Abaqus 约束识别、匹配与非匹配场量、倾斜初始间隙、双切向摩擦、有限滑移跨面、真实当前面积和法向演化、累计滑移、逐增量法向旋转、反向再粘着、客观历史、切线、事务、重启动和 MPI 等价 |
| HEX8 热力体算子 | `b49.hex8_c3d8t_operator`、`b50.hex8_c3d8t_capacity`、`b51.hex8_c3d8t_finite_heat`、`b519.hex8_c3d8t_finite_selective` | Abaqus C3D8T 的 32 自由度切线、八点非仿射应变及热流、逐列瞬态热容、有限变形后导热构形和畸变有限应变选择性体积积分识别；热膨胀八节点平均温度、角点节点热容和有限应变当前构形导热均由生产内核复现；体热源、变形后热容及对流构形仍保持单独识别边界 |
| HEX8 热接触 | `b52.hex8_c3d8t_thermal_contact`、`b520.hex8_c3d8t_gap_conductance`、`b521.hex8_c3d8t_thermal_contact_path`、`b522.hex8_c3d8t_faceted_thermal_contact` | Abaqus 固定间隙串联热阻，导热系数对间隙、压力和平均温度的局部导数，非匹配面的闭合、跨面滑移、开放和再接触，事务回滚、重启动、分布式等价，以及三分片曲面上的压力相关热流；B5.22 的机械矢量合力差异作为明确的 qualified 边界保留 |
| Coulomb 摩擦 | `m51.friction` | 粘着、滑移、反向再粘着、局部切线、守恒和 MOOSE |
| 完整链大滑移搜索 | `m52.large_sliding` | 跨多段所有权、力连续、MPI 等价、重启动和 MOOSE |
| 自动罚刚度和增广法 | `m54.augmented_contact` | 串联刚度、乘子事务、穿透门槛和约束极限 MOOSE 对比 |
| Backward Euler 热瞬态 | `m21.transient` | 一致热容、绝热升温、非均匀制造解和 MOOSE |
| J2、Norton 和耦合本构 | `m22.inelastic` | 闭式根、卸载—再加载、极端尺度、局部切线和 MOOSE |
| 瞬态非匹配 PCMI | `m23.pcmi` | 节点场、压力、总力、平均量和 40 个积分点 MOOSE 对比 |
| 时间误差和全局诊断 | `m53.time_integration`、`m32.diagnostics`、`m40.foundation` | 100 秒时间步研究、step-doubling、分场 Jacobian 和守恒 |
| 当前构形有限应变 | `m41.finite_strain` | Taylor/Rashid 局部解析与非匹配 PCMI MOOSE 对比 |
| follower pressure | `m42.follower_pressure` | 四边法向、当前合力、几何刚度和三类边界 MOOSE 对比 |
| 非共轴有限转动 | `m43.noncoaxial_finite_strain` | 四单元 100 步路径、材料分支和共享状态重放 |
| 综合瞬态大滑移热—摩擦接触 | `m57.integrated` | 完整当前法向 MOOSE 对标和热接触跨段所有权；1-rank/2-rank/4-rank 完整状态等价保留为手动 benchmark |
| 分布式装配和影子态 | `m34.parallel`、`m55.shadow_state` | 默认 CTest 中隔离路径的 1-rank/2-rank 状态与贡献区间等价 |
| 迭代求解测量 | `m56.iterative_solver`、`performance.m34` | 固定 CPU、固定线程的中型和较大网格观测 |
| 明确适用边界 | `scope.boundary` | 未实现物理和未鉴定组合的发布级清单 |

每次发布候选必须运行验证矩阵审计。该测试会拒绝缺失的必需行、未知 CTest、
重复标识、空证据和不存在的证据文件。固定依赖、Release 构建、警告作为错误
和检测器配置见 [`docs/reproducible-build.md`](reproducible-build.md)。

### 12.3 MOOSE 参考溯源

所有 fuelsim-to-MOOSE 对比读取 `verification/moose/` 下由对应 MOOSE 输入
生成并追踪的 Exodus 网格，不在比较测试中重建等价坐标。参考目录记录：

- MOOSE 输入和生成命令；
- MOOSE 源提交和可执行文件 SHA256；
- 网格、CSV 及其他参考结果的 SHA256；
- 物理、网格、罚参数、加载路径和输出量口径；
- qualified 门槛的原因、实测相对误差、绝对误差和适用边界。

参考文件变化必须先定位来源，再更新校验和和验证矩阵。构建成功或文件可解析
不等于数值验收通过。

## 13. 已验证边界之外

当前证据不能支持以下外推：

- 3D、非轴对称载荷、非 Quad4 体单元、闭合接触面、自接触或 mortar 机械
  接触；
- 二氧化铀或包壳的燃耗、辐照、裂变气体、肿胀、开裂、重定位、氧化或
  冷却剂通道工程关联；
- 任意转角、任意加载路径和任意网格上的一般有限转动鉴定；
- 超出 M5.7 规定网格、载荷路径和接触参数的大滑移与 Coulomb 摩擦组合；
- 当前构形 traction 方向随法向旋转的一般矢量 follower law；
- 代表性 100 秒 PCMI 之外的长期时间精度，或未做独立网格、时间步和接触
  参数研究的新工况；
- 迭代求解器在任意网格上的扩展性，或当前影子自由度减少等同于进程常驻
  内存同比例下降；
- 跨检查点格式版本迁移、分布式网格、分布式材料历史或并行 Exodus I/O；
- 试验数据验证、软件质量保证流程或法规鉴定。

`m43.noncoaxial_finite_strain` 保留显式 qualified 门槛：规定换向路径的应力
和弹性应变低参考尺度点使用 `0.6%` 最大逐点相对误差门槛，纯蠕变隔离路径的
应力过零点使用 `4%` 门槛；两者都不增加分母下限，聚合误差和其他量仍使用
验证矩阵规定的较紧门槛。详细实测值和误差位置见
[`docs/m4.md`](m4.md) 与 [`docs/verification.md`](verification.md)。

任何新增工程物理都必须先形成具体局部残量、ADlite 一致切线和独立参考，
再进入验证矩阵。在证据齐全前只能标为 `qualified` 或 `limitation`，不能写成
当前已验证能力。
