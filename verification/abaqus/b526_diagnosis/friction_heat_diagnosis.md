# B5.26 摩擦反转热反力误差的独立审计

B5.26 原参考包含 Abaqus 默认摩擦生热，而 Fuelsim 只有接触导热和机械摩擦耗散。
这项物理定义差异已经由输入、文档、全历史能量和独立原生诊断共同确认。
只关闭 Abaqus 的摩擦生热后，同一份历史 Fuelsim 生产结果通过了既有完整比较规则，
热反力最大逐点相对误差从 75.914168697% 降至 0.071710348%，没有调整验收门槛。

此处的历史生产代码为 `e495bb821fbc039f2a1adff314d8ec0c9106582f`，使用相同 GCC 14.3.0、
ADlite 0.2.3、PETSc 与链接时优化配置运行当前完整 Fuelsim 输入。
该对比排除了固定初始质量修改引起原 28 条失败的可能性。
随后热容规则回退后的当前生产程序，也实际运行并通过了
`fuelsim_b526_hex8_c3d8t_friction_reversal_abaqus_tests`。
当前生产验证与上述历史对照分别保留证据，不能将历史指标当作当前复算的具体数值。

原生与 Fuelsim 均为两个接触长方体，每个体沿法向两层 C3D8T，合计 24 个节点、4 个单元。
主侧外端节点 1、4、7、10 固定 300 K；从侧外端升至 400 K，先压紧、正向滑动、反向滑动，
最后重新黏着。全部 20 个时间步均为 0.02 s。
热接触导通率为 `h=50+0.001*p`，摩擦系数为 0.05。
原生全历史压力约为 14359—75796 Pa，处于 0—500000 Pa 的导通率表范围内；
温度处于 300—400 K，也在本例 300—500 K 的材料表内。

原输入没有 `*GAP HEAT GENERATION`。
本地 2018 `SIMACAEKEYRefMap/simakey-r-gapheatgeneration.htm` 第 31、38、41 行明确：
默认热转换比例为 1，热量平均分给两侧；该关键字用于修改默认行为。
同版本 `SIMACAEITNRefMap/simaitn-c-thermalinteraction.htm` 第 415—433 行
明确将这一规则用于全耦合温度—位移中的摩擦滑动。
[官方关键字文档](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-gapheatgeneration.htm)
和[热接触说明](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEITNRefMap/simaitn-c-thermalinteraction.htm)
给出相同默认定义。

Fuelsim 的 [quad4_face.cpp](../../../elements/src/quad4_face.cpp:829) 将接触热流计算为
`h*(T_secondary-T_primary)`，热残量在两侧成对相反；机械摩擦耗散单独记录，未加入热残量。
因此两者在黏着且无摩擦耗散时可以吻合，在滑动产生热源后会偏离。
原生 `ALLFD` 在第 8 个增量首次非零，该步开始两侧界面温度同时明显高于 Fuelsim；
第 10 步两侧平均温度的差分别为约 0.432263 K 和 0.436817 K。
这种两侧共同升温与热源差异一致。

对归档原生输出做了独立能量检查。采用其原生 `IVOL` 积分体积、`TEMP` 材料点温度及前一步
温度，按本例已输出的角点配对关系计算：

```text
Q_storage = sum(IVOL * 100 * [1+0.001*(TEMP-300)] * [TEMP_new-TEMP_old] / 0.02)
Q_boundary = sum(RFL11)
Q_friction = [ALLFD_new-ALLFD_old] / 0.02
```

所有 480 条节点记录、640 条活跃材料点记录和 20 条能量记录均参与计算。
这里直接使用原生 `IVOL`，没有用 Fuelsim 的几何权重替换原生体积。
全 20 步 `Q_storage-Q_boundary-Q_friction` 的最大绝对残差仅为 `5.7478253e-5 W`。
例如第 10 步，储热为 2409.241540759 W，边界热输入为 1976.002535574 W，
摩擦生热功率为 433.239030838 W，剩余误差为 `-2.5653597e-5 W`。
全历史储热减边界热输入的积分为 39.0203215444 J，`ALLFD` 末值为 39.0203208923 J，
相差 `6.52e-7 J`。这确认了原参考确实把摩擦耗散计入热方程。
逐步证据见 [friction_heat_energy_audit.tsv](friction_heat_energy_audit.tsv)，
数据来源及散列见 [friction_heat_audit_provenance.json](friction_heat_audit_provenance.json)。
其中 Fuelsim 列对应审计时散列锁定的输出快照；后续生产复算已经覆盖了相同路径的结果文件，
归档表格仍保留原快照的数值。

边界热反力的提取和比较没有发现新的单位、节点或时间错误。
[extract_b524_b525.py](../extract_b524_b525.py:136) 将原生 `RFL11` 按节点原样写入
`reaction_heat_flux_w`；能量 CSV 中的边界热率就是这些值的和。
比较器按增量、节点和 `1e-7 s` 时间容差对应，没有给相对误差分母加下限。
这 28 条失败都在规定 300 K 的主侧外端，均同时超出原 1% 相对或 0.2 W 绝对限定。
最后一步节点 1，原生热反力为 `-0.7917937157378266 W`，历史 Fuelsim 为
`-0.2715601991636575 W`，绝对误差为 0.520233516574 W。
同一步邻近节点 2 的原生温度为 300.105397746 K，当前固定初始质量 Fuelsim 为
300.036163502 K；边界热流偏差对应实际内部温度差，不是单独的反力输出错误。

另发现一个独立的输出口径和覆盖问题：接触 CSV 列 `contact_heat_flux_w` 实际直接提取原生
接触 `HFL`，其单位是 W/m²，表示离开从侧表面的净热流密度。
2018 热接触文档第 508—514 行区分了 `HFL` 与乘以节点面积后的 `HFLA`。
在存在摩擦生热时，`HFL` 还不能直接等同于 Fuelsim 的纯导热总功率；例如第 10 步
原生接触 `HFL` 的节点平均值约为 `-185.31355 W/m²`，而 Fuelsim 的总导热功率为
`+31.88518 W`，其中有摩擦热分配造成的物理差异，不能只凭符号判错。
当前比较器读取接触 `HFL` 列但未用于任何误差指标，因此这个命名问题没有造成 28 条
`RFL11` 边界失败，同时也意味着接触热流尚缺独立的数值比较。
后续应显式区分热流密度和热率，并按原生 `HFLA` 或明确面积积分补充同口径比较；
本次没有改提取器字段、比较器或门槛。

[完整关闭摩擦生热输入](b526_friction_reversal_no_friction_heat.inp)仅增加以下两行：

```text
*Gap Heat Generation
0.0, 0.5
```

该原生诊断在 `2026-09-13T07:46:23.4847114Z` 完成。节点、网格、材料、初始密度、
热接触、摩擦系数和全部时间历程均与原卡相同；原卡每一行都保留。
机械摩擦与 `ALLFD` 仍存在，只关闭其向热方程的转换。
相同历史 Fuelsim 结果在修改前后参考下的主要指标如下，表中均为百分比。

| 指标 | 原默认摩擦生热参考 | 关闭摩擦生热参考 |
|---|---:|---:|
| 温度 L2 | 0.079855942 | 0.000322247 |
| 温度最大逐点 | 0.240546287 | 0.001360668 |
| 位移向量 L2 | 0.038230908 | 0.006174861 |
| 位移向量最大逐点 | 0.691266007 | 0.198688916 |
| 应力张量 L2 | 0.166007982 | 0.015330973 |
| 应力张量最大逐点 | 0.325112928 | 0.066227478 |
| 热反力 L2 | 0.025640510 | 0.024520784 |
| 热反力最大逐点 | 75.914168697 | 0.071710348 |
| 超过热反力原限定的样本数 | 28 | 0 |

新参考下既有完整比较器实际输出 `PASS`，接触状态匹配率为 100%。
完整日志保存在 [no_friction_heat_full_field_comparison.log](no_friction_heat_full_field_comparison.log)。
上述通过结论仅覆盖现有比较器真正检查的场，不能据此声称被忽略的接触 `HFL` 已经过数值验收。
当前生产程序通过同一原生验收的记录见
[current_production_focused_ctest.log](current_production_focused_ctest.log)，该次五项针对性测试全部通过。

原 active 输入及四份 CSV 已逐字归档于
[b526_friction_heat_history](../b526_friction_heat_history/)，并保存局部 `SHA256SUMS`。
当前 active 输入等同于完整诊断卡，四份 active CSV 是新原生结果的逐字复制，仅修改外部文件名。
保留原生运行记录及文件名映射，见
[active_reference_provenance.json](active_reference_provenance.json)。
`generate_b524_b525.py` 仅对 `b526_friction_reversal` 补入同一设置，未改其他算例。
没有修改 Fuelsim 生产代码、输入或任何验收门槛，也没有把这项修复推广到其他综合算例。
