# 生产结果字段

以下字段由生产程序 `fuelsim -i <case.fsi>` 写入 `[Outputs]` 指定的 Exodus
文件，数值采用 SI 单位。变量名称包含输入卡中的接触对名称，结果文件为名称
预留 256 个字符；名称超出文件容量时明确报错，不允许截断后混淆接触对。

## 材料与反力

瞬态结果对每个材料积分点输出 `stress_<component>_q<index>`、
`elastic_<component>_q<index>`、`plastic_<component>_q<index>`、
`creep_<component>_q<index>`、`equiv_plastic_q<index>` 和
`equiv_creep_q<index>`。轴对称分量为 `rr,zz,hoop,rz`，三维分量为
`xx,yy,zz,xy,yz,xz`；剪应变均为张量分量，不是工程剪应变。积分点下标从零
开始，采用生产单元的内部积分点顺序，不等同于其他求解器的积分点编号。
三维一点减缩积分单元目前将同一个材料状态写入八个位置编号，其余八节点单元
有八个独立材料点，二十节点单元有 27 个独立材料点。

二十节点单元的瞬态结果还输出每个材料积分点的 `current_x_q<index>`、
`current_y_q<index>` 和 `current_z_q<index>`。坐标由参考积分点加上完整
二次位移形函数插值得到，单位为米，可用于与外部求解器逐单元唯一匹配
材料积分点；它们不是温度角点坐标，也不代表积分体积。

八节点单元的瞬态结果还输出 `reference_x/y/z_q<index>`、
`current_x/y/z_q<index>`、`material_temperature_q<index>`、
`heat_flux_x/y/z_q<index>`、`logarithmic_strain_<component>_q<index>` 和
`infinitesimal_strain_<component>_q<index>`。热流采用小应变参考构形、有限应变
当前构形的温度梯度；完整积分材料温度来自对应角点，减缩积分采用体积加权
单元材料温度。对数应变由总变形的左伸长张量计算，有限应变完整积分采用
单元平均对数体积修正；线性化应变采用参考位移梯度及完整积分的平均迹修正，
不得把有限应变下的线性化应变当作材料历史。小应变模型没有限制总变形梯度，
当派生对数应变没有定义时该字段为空值，不因此拒绝有效的小应变结果。

`integration_measure_q<index>` 是与该材料点对应的离散体积：完整积分有限
应变使用参考点体积分数乘单元当前总体积，减缩积分为整个单元体积。
`current_measure_q<index>` 另外给出该位置梯度计算的局部构形测度。减缩积分
的全部派生字段与材料历史一样在八个编号中重复，不能按八点再次求和。

轴对称瞬态积分点还输出 `reference_r_q<index>`、`reference_z_q<index>` 和
`reference_measure_q<index>`。前两个字段为参考构形坐标，第三个是完整的
`2*pi*r*detJ*w` 参考积分测度；它们可用于按参考体积计算逐单元材料平均值，
不得当作有限应变的当前体积。摘要中的 `minimum_accepted_time_step`、
`maximum_accepted_time_step` 和 `maximum_accepted_time_error_estimate` 只统计
本次执行接受的时间步；续算时不包括检查点之前的步骤。

瞬态节点输出还包含 `reaction_heat_flux` 和 `reaction_force_<direction>`。
它们是接受这个时间步时、提交材料状态之前计算的未施加 Dirichlet 行替换的
物理残量；不是缩放后的求解器残量。受约束自由度上的值是热反力或机械反力，
自由自由度上的值是平衡残差。输出不重新执行本构更新，也不修改已提交历史。
初始帧使用初始保存的零残量。二十节点单元非温度自由度节点的热反力为空值，
不能用中间节点的插值温度误认为该处存在独立的温度自由度。

三维瞬态节点还输出 `dirichlet_temperature` 和 `dirichlet_displacement_x/y/z`，
用零或一明确标识对应场是否施加位移或温度约束。二十节点非温度自由度节点的
温度约束标识为空值。对标受约束反力时可以使用这些标识，不能通过反力数值大小
猜测边界条件。

瞬态全局字段 `conservation_<name>` 保存接受这个时间步时的 24 项守恒诊断，
名称与工程历史一致，包括热率、机械功增量、弹性能变化、各机制耗散与沙漏能。
它们直接读取已提交诊断，不重新装配残量或更新材料。输出间隔大于一步时，
每个输出帧只保存最近一个接受步的诊断，不能将稀疏输出帧求和当作全程累计能量。

## 接触节点

接触节点变量统一命名为 `contact_<field>_<pair>`，其中 `<pair>` 为输入卡中的
接触对名称。非该接触对节点上的值为空值。`projected` 在接触约束节点上始终
写为零或一，因而可以区分失去投影与不属于该接触对；不能把两者都当作零压力。

共同字段包括 `gap`、`pressure`、`tangential_traction`、
`elastic_tangential_slip`、`sliding`、`projected`、`tributary_area`、
`normal_force` 和 `tangential_force`。后两个字段是该节点的标量接触力，不是
压力；`sliding` 用零或一标识未滑动与滑动。

轴对称接触还输出 `primary_segment`、`tributary_length`、`current_r`、
`current_z`。`primary_segment` 是生产接触链中的零起始线段下标，仅在成功
投影时有值，可用于检查跨线段后投影归属的改变。

CAX8T（八节点轴对称温度—位移耦合单元）的 secondary 接触边，即节点约束所在侧，
另外输出 `contact_recovered_pressure_<pair>` 和 `contact_recovered_shear_<pair>`，
单位均为 Pa。前者为正压缩的恢复压力，对应 Abaqus 的 `CPRESS`；后者为沿当前
primary 边切向的物理剪切应力，对应 `CSHEAR1`。原始 `tangential_traction` 使用
残量符号，因此恢复前须反号才能与这里的物理剪切应力比较。

恢复使用二次边节点值的线性最小二乘投影，在相邻边共享节点处取算术平均，并对
恢复结果应用经 Abaqus 例题识别的极值限制。恢复只在写结果时计算；原始 `pressure`、
`tangential_traction`、节点力、接触活动状态和历史仍保持其原有物理含义。
局部接触时，恢复压力可延伸至没有节点接触力的相邻节点，不能用恢复压力乘节点
面积重建接触力，也不能用它判断该节点是否真正接触。完全分离时两项恢复值为零。
不属于 secondary 接触边的节点保持空值。适用证据见
[CAX8T 接触恢复验证](../verification/abaqus/b114_cax8t_recovery_validation.md)。

三维接触还输出 `primary_face`、`constraint_pressure`、`current_x/y/z`、
`normal_force_x/y/z`、`tangential_force_x/y/z`、`total_slip_x/y/z` 和
`elastic_slip_x/y/z`。`primary_face` 是生产搜索中的零起始面下标，不是 Exodus
边集编号。`constraint_pressure` 是原始积分约束压力，`pressure` 是节点恢复
压力；部分接触下二者可能有不同的正值节点范围，不得相互替代。

原有全局变量 `contact_force_<pair>`、`contact_tangential_force_<pair>` 和
`contact_heat_rate_<pair>` 继续给出该接触对的标量总量。三维矢量合力应从完整
节点力分量求和，不能把标量压力合力作为任意方向的矢量合力使用。
