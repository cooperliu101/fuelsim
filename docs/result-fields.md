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

三维接触还输出 `primary_face`、`constraint_pressure`、`current_x/y/z`、
`normal_force_x/y/z`、`tangential_force_x/y/z`、`total_slip_x/y/z` 和
`elastic_slip_x/y/z`。`primary_face` 是生产搜索中的零起始面下标，不是 Exodus
边集编号。`constraint_pressure` 是原始积分约束压力，`pressure` 是节点恢复
压力；部分接触下二者可能有不同的正值节点范围，不得相互替代。

原有全局变量 `contact_force_<pair>`、`contact_tangential_force_<pair>` 和
`contact_heat_rate_<pair>` 继续给出该接触对的标量总量。三维矢量合力应从完整
节点力分量求和，不能把标量压力合力作为任意方向的矢量合力使用。
