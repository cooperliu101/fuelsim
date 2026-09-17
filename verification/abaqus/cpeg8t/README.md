# CPEG8T 与二维热机械接触验证

本目录包含完整 Fuelsim 生产输入、Exodus 网格、Abaqus 2025 原生输入及结果，
并保留开发时使用的独立算子探测。参考于 2026-09-17 使用单核 Abaqus 2025 生成。

## 生产比较

| Fuelsim 输入 | Abaqus 参考 | 覆盖范围 | 最大非零逐点相对误差 |
|---|---|---|---|
| `prescribed.fsi` | `small_prescribed.inp` | 小应变受约束热膨胀、伸长和两个弯曲参数 | 4.75e-8 |
| `free_controls.fsi` | `free_controls.inp` | 非均匀温度、自由伸长和自由弯曲、导热反力 | 3.57e-8 |
| `finite_prescribed.fsi` | `finite_reference_node_probe.inp` | 有限应变、厚度变化、弯曲、瞬态热容 | 3.98e-8 |
| `contact.fsi` | `contact.inp` | 小应变、无摩擦面到面接触、间隙传热 | 1.67e-8 |
| `finite_contact.fsi` | `finite_contact.inp` | 有限应变热机械接触、截面反力 | 1.70e-6 |

误差是无量纲比例，较大的参考控制量误差主要受 Abaqus 历史输出精度影响。
检查范围包含全部有效节点的温度、两个位移、面内反力和热反力，每单元九个材料点的
四个应力分量，每个 section 的三个控制量及三个反力，以及接触间隙、压力、合力和传热量。
零场根据约束和矩形对称性预先指定，单独检查绝对误差；不调整相对误差分母。
所有非零场的相对二范数、相对绝对峰值和最大逐点相对误差均要求小于 0.5%。

接触模型由两个 0.02 m × 0.01 m 的矩形组成，初始间隙 0.0001 m，
厚度均为 0.1 m，弹性模量 1 MPa，泊松比 0.25。底边固定，顶边下压 0.0003 m，
全部水平位移及截面控制量固定。上下温度为 400 K 和 300 K，导热系数为
10 W/(m K)，接触导热系数为 1000 W/(m² K)，罚刚度为 1e9 Pa/m。
小应变结果为接触压力 11320.7547 Pa、合力 22.6415094 N、传热量 66.6666667 W。
有限应变例题采用一个相同载荷增量，并关闭热容项以比较准静态热平衡。

必须在 Abaqus `*SURFACE INTERACTION` 的数据行显式写入厚度 `0.1`。
其默认值为 `1.0`，不会自动继承本例体单元厚度；遗漏会改变接触面积及传热量。
该设置依据[官方接触属性说明](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-surfaceinteraction.htm)。

## 自动测试与运行

Python 比较脚本需要 NumPy 和 netCDF4。配置时可通过
`-DFUELSIM_TEST_PYTHON=/path/to/python3` 指定已安装这两个包的解释器；
它只用于测试，不是生产程序的数值依赖。

```bash
build/fuelsim -i verification/abaqus/cpeg8t/contact.fsi
ctest --test-dir build -j8 -R 'cpeg8t|plane_transaction' --output-on-failure
```

`compare.py` 将完整输入卡及网格逐字复制到隔离目录，运行生产程序，只读取生产结果
与已保存的 Abaqus 字段，不重新求解参考问题。每个比较目录保存 `production.log`
和逐字段的 `comparison.json`。日常 CTest 不调用 Abaqus。
本次五个算例的完整逐字段记录汇总在 `validation.json`。
输入已经使用扁平的 `[GeneralizedPlaneStrain]` 多 section 语法，网格不含参考节点或
自定义元素属性。原生 Abaqus 仍使用其参考节点定义；比较时将位移及弯矩换算到
Fuelsim section 的初始面积形心。有限应变接触的最大误差来自这一弯矩换算，
原生单精度历史值相减后相对误差放大，但仍远小于 0.5% 门槛。

`sections.fsi` 由真实生产程序运行两个时间步：下块厚度 0.1 m，三个约束按线性函数变化；
上块厚度 0.2 m，伸长自由，两个转角在第二个时间步开始变化。
`check_sections.py` 检查两个 section 的全部控制量、截面力和弯矩，以及十八个
材料点的四个应力分量，并与小应变弹性解析解比较。它不使用 Abaqus 参考。
内部测试另外覆盖多 block 共用 section、初始面积形心、重复或遗漏归属、未知 block、
非正厚度及共用场节点时禁止不同 section。

此前经用户明确授权，仅将 `finite_contact` 双进程比较中九个材料点的理论零剪应力
`stress_xy` 绝对门槛统一为 1e-7 Pa；两侧零值及两侧差值均需满足该门槛。
节点自由度、其他材料历史、非零场和接触量的门槛保持不变。
独立提交版本的统一回归为 **335/335 通过**，用时 151.90 秒，其中包含九项 CPEG8T 相关测试。
该版本仅包含本功能及必要的并行状态提交支持；此前包含其他未提交工作的工作区为
336/336 通过，两次验证的范围和指标均保存在 `validation.json`。
`reference.sha256` 固定五个生产比较的原生输入与字段文件校验值，运行比较前先检查完整性。
两个 MPI 进程的有限应变接触测试同时与 Abaqus 和单进程结果比较，逐项检查场值、
全部材料历史和接触输出；内部测试检查贡献区间不重复及所需影子自由度的完整性。
内部测试另检查材料与边界切线、曲边接触切线及守恒、负厚度拒绝、检查点与完整状态回滚。

Abaqus 参考重算命令为：

```powershell
.\run.ps1 -SourceDirectory <本目录的绝对路径> -JobName contact
```

完整实现边界见[二维广义平面应变说明](../../../docs/generalized-plane-strain.md)。
上述表格记录首轮五个例题。后续补充了有限应变非弹性历史、非匹配曲边反力、
跨单元滑移分段积分、热算子鉴定、检查点、自适应时间步和收敛性验证。
具体通过范围及尚未通过的 Abaqus 比较见 [补充验证记录](coverage.md)。
`extended_reference.sha256` 固定新增原生证据，原有校验文件和误差门槛保持不变。
逐字段补充结果，包括未通过项，保存在 `extended_validation.json`。

## 独立原生探测

`reference_node_probe.inp` 使用一个 0.02 m × 0.01 m 的矩形 CPEG8T，
初始厚度 0.1 m，参考节点位于截面原点。全部面内位移固定，温度从
300 K 上升到 400 K。参考节点的厚度伸长为 0.0001 m，两个转动参数为
0.002 和 -0.003。弹性模量 1 MPa，泊松比 0.25，热膨胀系数 1e-5/K。
小应变下厚度方向应变为 `0.001 + 0.02*y + 0.03*x`。
原生参考节点反力为 `RF3 = -0.16 N`、`RM1 = 4e-5 N m`、`RM2 = -2.4e-4 N m`。
此探测用于确定参考节点约定；生产比较另使用表中相同载荷及热求解类型的参考。

有限应变单步探测识别的厚度应变增量为 `2*(t-t_old)/(t+t_old)`。
面内虚轴向梯度使用参考转动参数本身，参考控制量内力则除以当前厚度。
`finite_thickness_probe.inp` 在初始厚度 0.2 m 下独立复核该区别。

程序在 Windows 临时目录运行单核 Abaqus，复制文本日志和原生字段提取结果到本目录。
二维结果数据库的 `RF` 字段只包含面内两个分量，参考节点的第三方向反力
从原生 `RF3` 历史输出提取，不用面内场重构。

理论依据为[广义平面应变理论](https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-genplanestrain.htm)
与[二维单元库及自由度定义](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEELMRefMap/simaelm-r-2delem.htm)。
