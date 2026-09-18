# 小间隙传热与首步塑性验收

用户指定压力相关传热在较小距离内验证，并明确以首步即塑性例题通过作为
非弹性验证通过的依据。本次没有修改生产算法或任何误差门槛。

## 压力相关传热

新增完整生产输入卡 `thermal_pressure_small_gap.fsi`。沿用两个宽 0.02 m、
高 0.01 m、厚 0.1 m 的矩形 QUAD8，初始间隙从原例的 100 μm 改为 1 μm。
两块均为小应变，上侧作规定的竖向刚体平移；两侧温度角点仍由热方程求解，
仅底边保持 300 K、顶边在 1 s 内从 300 K 升至 400 K。没有热膨胀。

间隙在 0—0.25 s 保持 1 μm，0.5 s 达到 -10 μm，0.75—1 s 再保持 1 μm。
负间隙表示罚函数接触允许的穿透。罚刚度为 `1e9 Pa/m`，最大接触压力为
10000 Pa。传热定律不变：`h = 1000 + 0.002*p + 2*(T_average-300)`。
体导热系数为 10 W/(m K)，关闭热容项，逐时刻求准静态热平衡。

五个张开时刻和三个闭合时刻分别检查全部有效节点温度、热反力、位移、机械
反力及两单元各九个材料点的应力。非零场的三种相对指标均沿用 0.5% 门槛；
零量沿用原有绝对门槛。完整原生比较通过，最大逐点相对误差约为 **0.005261%**，
对应 0.625 s 热反力，最大绝对差约 **0.001121 W**。

另以串联热阻作独立检查。两块的总热阻为 1 K/W，接触面积为 0.002 m²，
热端温升为 `100*t`，因此总热流为 `100*t/(1+1/(0.002*h))`。
解析检查、张开阶段原生比较、闭合阶段原生比较均通过。8 个时刻的间隙、压力、
传热系数、热流及全部误差记录在 `small_gap_acceptance_validation.json`。

这一结果仅覆盖本例的最大正间隙 1 μm，以及规定的几何和接触参数，不能据此
推断任意更大间隙或不同网格下的 Abaqus 启停距离。原 `thermal_pressure` 的
较大张开工况仍保留为定律差异诊断，不改写原始失败结果，也不再作为本次
用户限定范围内的未通过验收项。

## 首步塑性

重新运行 `inelastic_initial_plastic.fsi` 并启用 `--require-native`，要求解析
检查及 Abaqus 比较均通过。全部 13 个时刻、27 个材料点、分机制应变张量、
等效标量、节点场及截面量通过；原有输入、参考和门槛均未改变。

按用户决定，小应变非弹性验证记为通过。原 `inelastic_history` 仍用于观察
Abaqus 显式/隐式切换与 Fuelsim 后向 Euler 的积分差异，其原生比较失败记录
保留，不把该诊断结果伪装成通过。

## 重现

`create_pressure_small_gap.py` 只生成网格和 Abaqus 输入，不生成或修改 `.fsi`。
新网格、完整输入卡、原生输入及原生字段的摘要由 `small_gap_reference.sha256`
固定。原生任务使用一个 CPU，日常 CTest 只读取已保存的参考。

```bash
python verification/abaqus/cpeg8t/check_contact_laws.py \
  --case thermal_pressure_small_gap --executable /path/to/fuelsim --work /tmp/small-gap
python verification/abaqus/cpeg8t/check_inelastic.py \
  --initial-plastic --require-native --executable /path/to/fuelsim --work /tmp/initial-plastic
```

小间隙例题加入统一 CTest；首步塑性例题已在原有套件内。

以已提交版本 `5c287c7` 加本次验证改动作隔离 Release 构建，统一回归
**379/379 通过**，并发数为 8，总耗时 295.13 s。首步塑性原生比较的最大
逐点相对误差为 `4.76697e-8`，约 `0.000004767%`。构建基线、可执行文件
摘要和完整比较结果一并保存在上述 JSON 记录中。
