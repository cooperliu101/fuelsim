# 初始密度、热容质量与重力

输入材料的 `density` 是区域 `initial_temperature` 对应的初始密度，单位为
kg/m³。内置热物性函数不再接受 `density_temperature_coefficient`。导热率和
比热仍可随当前温度变化；自定义热物性函数的密度只在初始温度、初始时间 0
和参考位置求值。后续升温、时间推进和重启动不会把当前温度对应的密度重新
当作初始密度。

令 `J=dV_current/dV_reference` 为局部当前体积与参考体积之比，有限变形下：

```text
rho_current = rho_initial / J
dm = rho_initial * dV_reference = rho_current * dV_current
R_capacity_i = integral_reference(N_i * rho_initial * cp(T) * dT/dt dV_reference)
R_gravity_i = -integral_current(N_i * rho_current * acceleration dV_current)
            = -integral_reference(N_i * rho_initial * acceleration dV_reference)
```

这里 `dT/dt` 使用后向欧拉离散。热容保持各型号原有温度插值和积分规则，
但质量测度固定为参考构形：CAX2T_GPS、CAX4T、CAX4RT 使用参考形函数体积分
形成的节点对角质量；C3D8T 使用配对参考 Gauss 点的节点对角质量；C3D8RT
使用参考一致热容矩阵的行和；CAX8T/CAX8RT 和 C3D20T/C3D20RT 使用各自参考
积分规则的一致质量。比热在相应当前节点或积分点温度求值，温度导数完整保留。
固定质量热容对位移的导数为零。导热、体热源和表面热载荷继续遵循各型号
原有构形规则；没有把体积热源改成单位质量热源。

加速度按区域在 `[Regions]` 中设置。例如三维重力为：

```text
body_acceleration = 0 0 -9.81
```

三维依次为固定全局 x、y、z 分量，轴对称依次为径向、轴向分量：

```text
body_acceleration = 0 -9.81
```

单位均为 m/s²，省略时为零。当前入口支持区域内均匀、时间恒定的加速度；
稳态和瞬态都可使用。径向分量表示轴对称径向加速度场，不能用它表示固定
横向地球重力。程序仍求解准静态力学，没有增加位移惯性或摩擦生热。

重力使用位移形函数。二次单元保留有符号的一致节点力，不把力均分给角点；
一维广义平面应变单元的轴向力各分配一半至两个端截面控制节点。装配直接
使用上式密度与体积相消后的初始质量形式，避免两个几何导数相减造成舍入误差。
材料接口 `current_density` 提供保留自动微分导数的当前密度，并拒绝非正或
非有限体积比；原有当前、已接受和中间构形有效性检查继续执行。

工程历史记录 `body_force_work_increment`，即该时间步加速度体力所做的功。
两个半时间步被接受时，该项按两步求和，并进入机械功平衡。稳态结果也输出
节点反力，可直接检查支承反力与初始质量乘加速度的平衡。

检查点版本为 22，拒绝旧热容语义的检查点。区域加速度进入模型签名，改变
加速度后不能直接沿用旧检查点。

验证包含九种单元的小应变和有限应变局部总力、密度导数、热容与几何导数
检查，以及五种网格拓扑的完整生产输入 `steady_gravity_*.fsi`。生产重力
检查比较显著热膨胀后的约束总反力与初始质量乘加速度。原生 Abaqus 比较的
验收状态必须单独报告，不能用上述质量守恒检查替代应力或热反力比较。

旧版与新版 Abaqus 文档差异、原生热容证据及本次验收状态见
[密度模型核查](../verification/abaqus/density_mass_model_audit.md)。

Abaqus 2025 已完成原生重算。通过范围、未通过比较及非仿射热容离散差异见
[2025 验证记录](../verification/abaqus/2025/README.md)。

## Abaqus 2025 探测后的实施决定

统一采用初始质量的物理定义，保留上述九种单元的参考质量热容和重力实现。
Abaqus 2025 的热容离散不作为这一质量定义的替代：不引入 C3D8T 的整体/局部
体积比混合权重，不引入减缩积分单元随当前温度更新的热容密度，也不增加
切换到 Abaqus 热容的输入选项。完整原生证据见
[密度与热容实现探测](../verification/abaqus/2025/density_capacity_probe/README.md)。

验收分开记录：质量守恒、局部残量与切线以及生产重力反力由独立解析检查
验证；与 Abaqus 的热反力、应力和接触场仍按原门槛比较。热容离散差异可以
解释部分误差，但不能直接豁免应力、接触误差，也不能把失败比较改记为通过。
截至 2026-09-13，完整回归为 306/314 通过，八项原生比较仍未通过，逐项边界见
[2025 验证记录](../verification/abaqus/2025/README.md)。

实施决定后重新运行九种单元、五种生产重力输入、输入解析、检查点、三维求解器
和验证清单检查，20/20 通过，耗时 4.38 s。逐项结果见
[初始质量专项检查](../verification/initial_mass_policy_test_results.tsv)。本次仅补充
实施约束和验证描述，生产算法与上一轮完整回归相同。

用户随后明确授权放宽上述八项例题的比较容差。当前分项门槛见
[八项例题容差约定](../verification/abaqus/2025/tolerance_qualification.md)，
此前的八项失败及 306/314 统计保留为放宽前的历史记录。质量定义不变。
