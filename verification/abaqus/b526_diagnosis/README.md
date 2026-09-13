# B5.26 摩擦反转的摩擦生热诊断

[完整诊断输入](b526_friction_reversal_no_friction_heat.inp)与原 B5.26 摩擦反转卡逐行相同，
唯一新增内容是接触属性下的两行：

```text
*Gap Heat Generation
0.0, 0.5
```

第一项将摩擦耗散转为热的比例设为零，第二项保留两侧均分的默认分配值。
摩擦系数、机械摩擦耗散、热接触导热、初始密度和全部载荷历程保持原定义。
这是热方程物理定义对齐，没有拟合材料参数或修改误差门槛。

[run_no_friction_heat.ps1](run_no_friction_heat.ps1)逐字复制完整输入并调用 Abaqus，
复用上级 `extract_b524_b525.py` 提取全部 20 个增量，写入独立文件名。
本次实际原生运行已在 `2026-09-13T07:46:23.4847114Z` 完成；工作目录、输入和提取器散列
以及四份原生 CSV 散列见
[原生运行记录](b526_friction_reversal_no_friction_heat_provenance.txt)。

原基线输入和四份 CSV 已逐字保存在
[b526_friction_heat_history](../b526_friction_heat_history/)。
当前基线四份 CSV 是本目录新原生 CSV 的逐字复制，仅去掉文件名中的诊断后缀，
映射记录见 [active_reference_provenance.json](active_reference_provenance.json)。

历史 `e495bb8` 完整结果比较和随后热容规则回退后的当前生产验收均实际通过。
分别保留完整比较器日志与
[当前生产 CTest 记录](current_production_focused_ctest.log)，原有门槛没有改变。

原基线热源、接触热流口径、逐点热反力与关闭摩擦生热后的实际比较结果，见
[friction_heat_diagnosis.md](friction_heat_diagnosis.md)。
