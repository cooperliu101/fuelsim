# B5.23 原摩擦生热参考归档

这里逐字保留 2026-09-13 正式参考更新之前的 B5.23 完整输入和四份原生 CSV。
[SHA256SUMS](SHA256SUMS)记录这五个文件的原始散列。

这张输入没有 `*Gap Heat Generation`，因此 Abaqus 使用默认的摩擦耗散转热比例 1，
并把热量在界面两侧平均分配。Fuelsim 的该算例包含摩擦力与接触导热，但不含摩擦生热，
所以这组旧结果存在物理定义差异。它仍作为历史证据保留，不再是正式对比参考。

上级正式输入仅增加 `*Gap Heat Generation` 和 `0.0, 0.5`，其余定义不变。
四份正式 CSV 来自实际完成的无摩擦生热原生任务，并保持原始数值和字节内容。
同时修正了 Fuelsim C3D8T 非仿射有限变形下的热容与导热体积权重，修正后的实际
生产结果已经通过原有完整比较器；没有修改误差门槛。

详见[诊断报告](../b523_diagnosis/thermal_weight_diagnosis.md)、
[正式参考来源映射](../b523_diagnosis/active_reference_provenance.json)和
[完整通过日志](../b523_diagnosis/final_no_friction_heat_comparison.log)。
