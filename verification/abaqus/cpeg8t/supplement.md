# CPEG8T 组合工况补充验证

本轮增加十四张完整的静态生产输入卡，覆盖可变形曲面滑移、非恒定接触传热、
热膨胀与塑性蠕变共同作用、失败恢复与检查点续算，以及接触网格和罚刚度敏感性。
没有修改生产算法，也没有放宽任何验收门槛。完整分项数据保存在
[`supplement_validation.json`](supplement_validation.json)，原生输入和参考字段由
[`supplement_reference.sha256`](supplement_reference.sha256) 固定。

四项未通过比较的后续独立研究见[逐项研究记录](failure_diagnosis/README.md)。
该记录进一步区分了热接触平均方式、投影区间限制、进程数相关舍入以及
压力表在较小正间隙下仍传热的行为；尚未修改生产算法或改变这里的通过状态。

## 已通过的比较

### 间隙与温度相关传热

`thermal_clearance.fsi` 使用两个宽 0.02 m、高 0.01 m、初始厚度 0.1 m 的矩形。
下侧固定，上侧刚体平移，使间隙从 0.0001 m 增至 0.0004 m，再恢复。
下边界保持 300 K，上边界在 1 s 内升至 400 K；界面温度由生产求解器计算。
传热系数为 `1000 - 1e6*g + 2*(T_average-300)`，体导热系数为 10 W/(m K)。

每个非初始时刻的全部角点温度、热反力、位移、机械反力和九点应力通过 Abaqus 比较。
同时使用独立串联热阻解：界面平均温度为 `300+50*t`，接触面积为 0.002 m²，
两块材料的总热阻为 1 K/W，因此总传热率为 `100*t/(1+1/(0.002*h))`。
机械应力为零，接触处于张开状态。该例不是瞬态热容验证。

### 非弹性热机械接触完整历史

`coupled_contact.fsi` 使用相同的两块矩形。每个截面的伸长在 1 s 内从零增加到
0.003 m，两个转动固定。上边界竖向位移在 0、0.25、0.5、0.75、1 s 分别为
0、-0.0003、-0.0005、0、-0.0007 m；上边界温度从 300 K 升到 400 K。
两块材料都包含热膨胀、线性硬化塑性和 Norton 蠕变，并实际计算瞬态热容。
采用小应变，避免把有限应变下两程序不同的质量规则混入本次组合验证。

17 个输出时刻、18 个材料点的应力、弹性/塑性/蠕变张量和两个等效标量，
全部有效节点场、两个截面控制量及反力、接触合力通过原有门槛。
两个程序的全部 18 个材料点都存在同一时间步内塑性与蠕变同时增长，接触经历
闭合、张开、再次闭合。最大逐点相对误差约为 **0.038533%**。

原点处参考点的历史力矩转换到形心时，两个单精度历史量相减曾造成一个约
`5.8e-9 N m` 的力矩误差达到 0.5376%。`coupled_contact_centered.inp` 把参考点
设在初始形心，直接输出该力矩，消除了相减误差。两个转动始终为零，因此
这不会改变厚度或任何物理约束；原点和形心参考输入的全部物理节点场、材料张量、
等效标量在每个时刻完全相同。两份原生输入和输出都保留，不能把该修正称为
生产算法修复。参考点坐标与截面伸长关系见
[Abaqus 单元维数说明](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEELMRefMap/simaelm-c-dimension.htm)。

### 失败后的恢复与串行续算

`coupled_contact_split.fsi` 在 0.5 s 写入包含全部材料历史的检查点。
`coupled_contact_failure.fsi` 从该状态计算到 0.5625 s，故意把非线性迭代次数
限制为一次，使进入接触张开过程的试探计算失败。生产输出仍停留在 0.5 s，
所有节点场、材料历史和接触输出与完整计算在该时刻一致，输入检查点字节不变。
随后 `coupled_contact_restart.fsi` 续算到 1 s，逐字段通过 `rtol=atol=1e-12`
的原有续算门槛。这证明了失败返回及磁盘续算，不替代自适应时间步缩小后的
同进程重试验证。

## 已执行但不能标为全部通过的比较

### 两侧同时变形的曲面滑移

`deforming_contact.fsi` 包含两块各有两个单元的物体。初始从侧界面为抛物线，
只驱动外边界，从侧上边界先压入，随后水平滑动 0.012 m，跨越主侧单元分界。
体单元使用有限应变，界面节点的位移和温度由求解器决定。

Abaqus 完成 16 个时间步。Fuelsim 在首个 0.0625 s 时间步内拒绝状态，原因是
`CPEG8T primary projection intervals overlap ambiguously`。相邻二次边变形后
不保证切线连续，按主侧最近点投影得到的热积分区间可能重叠。
当前明确拒绝有歧义的投影，没有夹持、重复积分或放宽检查。
这说明当前投影边界限制会影响普通可变形多单元表面，不能仅当作罕见折角问题。

`deforming_fixed_primary.fsi` 固定主侧全部位移，保留从侧曲面变形与滑移，
两程序均完成计算。位移和材料应力比较一致，但节点热反力的最大逐点相对误差为
**0.883628%**，超过 0.5%，最大绝对差为 **0.185953 W**。
完整比较还记录了极小节点反力的相对误差问题；没有按幅值自动改用绝对门槛。
本例不能标为热机械接触全场比较通过。

### 压力相关传热的张开行为

`thermal_pressure.fsi` 的仿射定律为 `1000 + 0.002*p + 2*(T_average-300)`。
间隙经历张开、闭合、再次张开。整个过程通过独立热阻解析解，闭合阶段通过
Abaqus 比较。张开阶段两种定律不同：Fuelsim 在压力归零后仍保留基础传热系数，
本次 Abaqus 压力表参考停止传热。完整原生比较因此失败，分别输出解析、闭合
原生比较和张开原生比较的状态。没有擅自添加开关或修改现有仿射定律。
本例界面压力均匀，尚不能鉴定非均匀平均力学压力与逐点热传导定律之间的对应关系。

### 四进程续算

单进程写入、四进程读取同一检查点并完成续算。最终输入卡的独立检查中，
约 `0.0024 Pa` 的非零应力出现约 `9.83e-12 Pa` 差异，超过既定
`rtol=atol=1e-12` 的逐字段门槛。该项是未通过的数值比较，不是启动失败，
也不是此前仅针对理论零剪应力获批的例外。

## 接触网格和罚刚度检查

`contact_refine_2/4/8/16.fsi` 固定平直主侧，对可变形抛物线从侧施加 0.0005 m
压入，分别沿接触方向划分 2、4、8、16 个单元。接触合力依次为
73.775081、73.727389、73.656314、73.613496 N，取自主侧全部节点的竖向反力总和。
8 到 16 单元的合力变化约 **0.05817%**，峰值压力从 43136.72 Pa 变为
43043.57 Pa，变化约 **0.2164%**。各网格整体竖向反力守恒。
前三级变化不单调，且这里只加密接触方向，因此不宣称统一收敛阶或连续体误差界。

`contact_penalty_low.fsi`、`contact_refine_8.fsi`、`contact_penalty_high.fsi`
分别使用 `5e8`、`1e9`、`2e9 Pa/m`。合力为 66.897206、73.656314、77.627598 N，
最大穿透量为 77.9106、43.1367、22.7939 微米。增大罚刚度确实减少穿透，但反力
仍有明显敏感性，不能把 `1e9 Pa/m` 当作与罚刚度无关的硬接触结果。
这些重复参数扫描仅作为手动研究，不加入自动测试套件。

## 复现入口

网格和原生输入由 `create_deforming_contact.py`、`create_contact_laws.py`、
`create_contact_convergence.py` 生成；它们不生成或改写任何 `.fsi`。
生产输入卡全部显式保存，可直接执行 `fuelsim -i <case.fsi>`。
原生参考通过既有 `run.ps1` 单核计算，并由 `extract.py` 读取原始字段。

```bash
python check_contact_laws.py --case thermal_clearance --executable /path/to/fuelsim --work /tmp/law
python check_contact_laws.py --case thermal_pressure --executable /path/to/fuelsim --work /tmp/pressure
python check_coupled_contact.py --executable /path/to/fuelsim --work /tmp/coupled
python check_contact_restart.py --executable /path/to/fuelsim --work /tmp/restart
python check_contact_restart.py --executable /path/to/fuelsim --work /tmp/restart4 --mpiexec /path/to/mpiexec
python check_contact_supplement.py --case deforming_contact --executable /path/to/fuelsim --work /tmp/deforming
python check_contact_supplement.py --case deforming_fixed_primary --executable /path/to/fuelsim --work /tmp/fixed
python check_contact_convergence.py --executable /path/to/fuelsim --work /tmp/refinement
```

已知未通过的命令继续返回失败，不将诊断结果伪装成通过。自动测试增加间隙传热
原生比较、非弹性组合工况原生比较和串行失败恢复三项。

## 回归范围

从已提交版本 `ae68002` 建立独立源码目录，只加入本轮验证文件和三项 CTest 注册，
没有包含共享工作区中的其他未提交改动。Release 构建最多使用四个编译作业，
`ctest --test-dir /tmp/cpeg8t-supplement-build -j8 --output-on-failure`
为 **354/354 通过**，用时 **172.14 s**。初次配置误选不含 `netCDF4` 的 Python，
改为 `-DFUELSIM_TEST_PYTHON=/home/cooper/miniforge/bin/python3` 后重新运行完整套件。
组合工况后续增加的初始状态检查另行通过定向测试。上述自动测试通过不改变手动
研究中明确记录的未通过状态。
