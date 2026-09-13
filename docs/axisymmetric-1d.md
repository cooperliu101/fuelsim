# 轴对称一维广义平面应变与 1.5 维计算

`axisymmetric_1d` 将多个径向切片放入同一个热力耦合非线性系统。生产入口仍是
`fuelsim -i <case.fsi>`，问题类型仍是 `SteadyProblem` 或 `TransientProblem`。
首版同时支持小应变和有限应变，复用现有热膨胀、J2 塑性、Norton 蠕变及
塑性与蠕变同时作用的材料更新。

## 运动学与自由度

每个轴向切片中的温度和径向位移仅随半径变化。芯块和包壳分别拥有轴向位移链，
同一实体的相邻切片共享端截面轴向位移。同一切片中的径向单元共享上下端截面
轴向位移，因此轴向应变在截面上均匀，数值由轴向合力平衡或给定端位移确定。

`cax2t_gps` 是 Fuelsim 自有型号名；`gps` 表示广义平面应变，不是 Abaqus
型号的别名。两节点径向体单元的局部自由度为：

```text
x = [T0, T1, ur0, ur1, w_lower, w_upper]
H = z_upper - z_lower
epsilon_small = [(ur1-ur0)/(R1-R0), (w_upper-w_lower)/H, (ur0+ur1)/(R0+R1), 0]
dV_reference = 2*pi*R*H*dR
```

应变分量沿用 `[rr, zz, hoop, rz]`，最后一项是张量剪应变。截面之间不装配
轴向导热和径向位移梯度产生的剪切应变。轴向位移沿切片高度线性插值；每片
均匀轴向应变可以随切片变化。完整棒长变化自然等于各片轴向伸长之和。

小应变在参考构形上以两个径向 Gauss 点积分热传导、热容、体热源及力学残量。
热膨胀使用每个径向单元的两个节点温度算术平均 `T_exp=(T0+T1)/2`，因此同一
单元的两个材料点使用同一个膨胀温度。平均范围不跨径向单元，也不跨轴向切片。
弹性模量、塑性、蠕变及热物性仍在积分点温度求值，热容采用一致矩阵。
环向机械应变按单元参考体积平均，即 `(ur0+ur1)/(R0+R1)`；其对两个径向
位移的导数均为 `1/(R0+R1)`。同一平均规则进入材料、虚应变、内力与切线。
这对应矩形 CAX4T 在广义平面应变运动下的力学约束；径向和轴向应变本来就是
常数，完整应变迹的剩余面内修正为零。不套用 CAX4T 的角点热算子或 CAX4RT
的沙漏控制，也不跨单元平均材料应力或历史。

有限应变采用三个对角伸长比及 Hughes-Winget 中间构形增量：

```text
lambda_r = 1 + du_r/dR
lambda_z = 1 + (w_upper-w_lower)/H
lambda_hoop_bar = 1 + (ur0+ur1)/(R0+R1)
delta_e_i = 2*(lambda_i_new-lambda_i_old)/(lambda_i_new+lambda_i_old)
```

环向先平均伸长比，再计算上式增量，不平均逐点增量。三个方向没有旋转。
增量与已接受的弹性、塑性、蠕变及本征应变共同构成材料
更新输入，累计材料应变不等于总伸长比的对数。当前、已接受及中间构形的
伸长比、积分点半径和体积必须为正，非法试探态进入线搜索或缩小时间步重试。
力学采用 Cauchy 应力，径向、轴向采用当前梯度，虚环向梯度为
`1/((R0+R1)*lambda_hoop_bar_current)`。力学材料点测度为参考测度乘
`lambda_r*lambda_z*lambda_hoop_bar`，即单元整体当前/参考体积比。
热传导、热容和体热源仍采用逐点真实伸长比 `(R+u_r)/R` 对应的当前几何和
体积，不用力学平均测度替换。两套测度和虚应变的全部几何导数进入切线。
本模型不表示二维局部剪切、弯曲和端部应力集中。

有限应变材料历史中的热应变更新使用 `eigenstrain(T_exp_new)-eigenstrain(T_exp_old)`，
保持前后均温一致。平均温度对两个节点温度的导数均为 1/2，同时保留积分点
温度对其他物性的导数。材料调用复用现有四应变分量加温度的五分量自动微分
接口，切线通过显式节点链
组装成 6×6 矩阵。单元接收已接受的节点状态和两个材料点历史，只返回试探历史；
不保存全局状态，不提交时间步。

截面轴向合力为 `N = sum(sigma_zz*W_mechanical)/H_current`。体单元对上下端
轴向残量分别贡献 `-N` 和 `+N`。小应变使用参考测度与参考高度；有限应变
采用前述力学测度与当前高度，等价于参考加权平均轴向应力乘当前截面积。
轴向平衡包含界面摩擦
的端点等效力，内部控制节点的平衡会同时改变相邻切片的轴向应变。弹性应变能
变化与塑性、蠕变耗散诊断同样使用力学测度；弹性能分别按当前和已接受构形积分。

## 热机械接触与摩擦

圆柱界面 `ring_gps` 接收外侧实体内表面作为 primary、内侧实体外表面作为
secondary。局部排列为：

```text
[T_primary, T_secondary, ur_primary, ur_secondary,
 w_primary_lower, w_primary_upper, w_secondary_lower, w_secondary_upper]
g = R_primary + ur_primary - R_secondary - ur_secondary
p = max(-penalty*g, 0)
```

法向接触采用显式罚刚度，不提供自动罚刚度或增广拉格朗日法。小应变两侧
使用同一参考面积 `2*pi*R_secondary*H`；有限应变两侧使用同一当前面积
`2*pi*(R_secondary+ur_secondary)*(H+w_secondary_upper-w_secondary_lower)`。
热流、径向接触力、轴向摩擦力分别严格成对守恒。
开放间隙仍可导热。导热率复用已有气隙和仿射接触导热律。

界面沿轴向使用两个 Gauss 点。每点按上下端位移插值计算芯块相对包壳的轴向
位移增量，再调用已有库仑摩擦返回映射。切向牵引按轴向形函数分配到四个轴向
端点自由度。不能只比较本切片的轴向应变差，也不能只在切片中点计算摩擦。

接触使用当前试探间隙和压力，保留压力对滑动极限的导数。材料历史、弹性切向
滑移和累计切向滑移由外层统一接受和恢复。摩擦耗散单独报告，首版不把它自动
转换为热源。

两侧参考轴向区间必须一一匹配，配对保持固定，接触输入为
`sliding = small`。这一限制可以与区域有限应变同时成立，但不支持
滑过切片边界后的重新投影、非匹配轴向分片和轴向端面接触。
`slip_tolerance` 是无量纲比例，乘该接触对 primary 切片平均参考高度得到
弹性滑移的长度阈值。界面只播种间隙、两侧温度、相对滑移增量及当前半径、
高度这六个局部量，经闭式链组装 8×8 切线，保留压力和面积的导数。

## 生产网格与边界

Exodus 使用二维 RZ 坐标和真实 `BAR2` 径向连接，在同一节点表保存独立轴向控制节点。
每个径向单元通过 `axial_lower_node`、`axial_upper_node` 元素属性引用两个
控制节点从一开始的内部编号，属性值必须能精确表示为整数。切片高度从这两个
节点的轴向坐标读取，径向节点位于参考切片中面。节点共享关系由编号给定，
不按坐标自动合并。径向节点与控制节点角色不能重叠，径向坐标严格递增。
允许内节点位于轴线上，生产装配会添加零径向位移约束并拒绝冲突的非零约束。
输入卡只声明物理设置，
不生成网格和切片几何。示例网格与离线生成脚本在
[`verification/meshes`](../verification/meshes/gps_meshes.py)。

全局场顺序为 `[T(径向节点), ur(径向节点), w(轴向控制节点)]`，三个场无需
等长。径向节点不因为轴向位置相邻而共享温度或径向位移。芯块与包壳使用独立
轴向控制节点，接触前各自的刚体平移需要明确约束。

径向端点边集使用 BAR2 的侧号 1、2。温度与径向位移约束施加在径向节点集或
端点边集，轴向位移约束施加在控制节点集。`axial_force` 对节点集中每个控制
节点施加以 N 为单位的轴向力，正号沿正轴向。圆柱压力、径向牵引、表面热流
及对流沿用现有输入；省略 `configuration` 时，小应变采用参考面积，有限
应变采用当前面积。详见[输入卡说明](input-card.md)。

## 可运行算例与验收

最终统一回归、正式误差和适用边界见[验收记录](../verification/axisymmetric_1d_validation.md)。

以下输入卡均为完整文件，可以直接交给生产程序：

| 输入卡 | 独立验收对象 |
| --- | --- |
| [uniform_small](../verification/fuelsim/transient_gps_uniform_small.fsi) | 均匀温升与轴向力，Abaqus CAX4T 小应变全场比较 |
| [uniform_finite](../verification/fuelsim/transient_gps_uniform_finite.fsi) | 有限应变温升与轴向力，同步增量的 Abaqus 全场比较 |
| [nonuniform_finite](../verification/fuelsim/transient_gps_nonuniform_finite.fsi) | 非比例径向位移与非均匀温度，Abaqus 有限应变应力、历史和反力比较 |
| [contact_small](../verification/fuelsim/transient_gps_contact_small.fsi) | 两种切片高度、径向导热及接触粘着、滑动、反向、张开和再接触 |
| [chain_small](../verification/fuelsim/transient_gps_chain_small.fsi) | 两条两切片轴向链，内部控制节点自由，解析解与原生 Abaqus 连续网格共同检验摩擦传力 |
| [chain_finite](../verification/fuelsim/transient_gps_chain_finite.fsi) | 总轴向伸长达到 10%，当前面积和分片轴向力平衡的解析检验 |
| [inelastic_finite](../verification/fuelsim/transient_gps_inelastic_finite.fsi) | 两个材料点十步同时产生塑性与蠕变，完整历史与后向 Euler 耦合解析解比较 |
| [steady_finite](../verification/fuelsim/steady_gps_uniform_finite.fsi) | 稳态有限应变自由热膨胀，五个加载步及结果写出 |

Abaqus 原生输入、提取脚本、完整参考字段及来源记录在
[`verification/abaqus/gps_1d`](../verification/abaqus/gps_1d/README.md)。

局部测试检查一致切线与中心差分方向导数、独立残量重复调用、自由均匀热膨胀、
轴向均匀应变、轴线、材料历史不变、热容和体热源守恒，以及接触粘着、滑动、
反向加载、张开再接触和两端反向滑移。

外部参考使用施加相同降维运动学约束的 Abaqus CAX4T。矩形单元在此运动下
的力学平均规则已保持一致。CAX2T_GPS 的材料参数使用插值温度，CAX4T 使用
配对角点温度，因此跨型号的局部力学投影比较限定为常数弹性模量和泊松比，
仍保留非均匀温度及随温度变化的平均热膨胀。CAX2T_GPS 自身的温变材料解析
和切线检查独立保留。热算子采用各自型号的定义，非均匀有限应变
算例不把热反力和体热流列为跨型号相等指标。原接触算例保留全部内外体的
160 个应力及弹性应变样本，并检查全部径向反力、轴向控制节点反力和截面力。
非零验收场分别记录相对
L2、相对绝对峰值与最大逐点相对误差，门槛均为 0.1%；零参考量单独检查 SI
绝对误差，不设置相对误差分母下限。

内部测试还覆盖不同场长度、贡献唯一分配、完整历史恢复、时间误差与检查点。
实际检查点续算与连续运行逐帧比较，有限应变多切片算例的单进程和双进程结果
逐字段比较，包括全部材料、接触、反力及守恒量，几何与 NaN 分布也必须一致。
生产检查器只读取实际输出，不调用求解库。统一自动验收入口为
`ctest --test-dir build -j4 --output-on-failure`。

## 资料来源与边界

[BISON 的分层一维介绍](https://mooseframework.inl.gov/bison/tutorials/layered1D_introduction.html)
说明了各片分别设置芯块和包壳广义平面应变、以径向传热为主要适用范围。
[MOOSE 广义平面应变说明](https://mooseframework.inl.gov/modules/solid_mechanics/generalized_plane_strain.html)
给出了截面均匀非零轴向应变的定义。本文的端截面位移链和摩擦等效力装配由
虚功原理推导，不表示直接复制了 BISON 的摩擦实现。
