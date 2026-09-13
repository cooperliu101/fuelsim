# B5.23 热反力差异诊断

本目录保留两个独立的原生诊断。关闭摩擦生热后，B5.23 的原有六个热反力样本仍同时超过
0.5% 相对误差和 0.05 W 绝对误差限定。主要原因是 C3D8T 非仿射有限变形下的热容和
导热积分权重，而不是摩擦生热。完整证据见[诊断报告](thermal_weight_diagnosis.md)。

- [完整无摩擦生热输入](b523_integrated_no_friction_heat.inp)只比原输入增加
  `*Gap Heat Generation` 和 `0.0, 0.5`。摩擦力学仍然保留。
  [运行脚本](run_no_friction_heat.ps1)逐字暂存输入，并复用上级 `extract_b523.py`。
- [完整非仿射容量输入](c3d8t_nonaffine_capacity_probe.inp)包含规则和畸变两种几何、
  每种八个单节点升温基，共 16 个独立单元。其
  [运行脚本](run_nonaffine_capacity.ps1)和[提取器](extract_nonaffine_capacity.py)
  不生成或改写输入。
- [综合例题代数分析](analyze_thermal_weights.py)只读取已有节点、材料点和能量输出。
  它比较候选权重，不链接生产库、不求解问题、不替代全场验收。
- [独立容量分析](analyze_nonaffine_capacity.py)从原生热反力恢复节点热容权重，比较
  参考角点、参考配对 Gauss 点、当前角点、当前配对 Gauss 点、当前行和以及总体积候选。

两个原生任务均已实际完成，原生输入与输出散列分别记录在同名前缀的
`_provenance.txt` 中。分析使用的完整文件散列见
[综合分析来源](thermal_audit_provenance.json)和
[容量分析来源](nonaffine_capacity_analysis_provenance.json)。
目录中保存的 Fuelsim 差异是散列对应的修正前输出快照，后续生产程序重新计算可能替换
`build/blackbox/b523` 下的结果，不能把这里的数值当作后续二进制的结果。

当前无摩擦生热数据仍为诊断资料，尚未替换上级 B5.23 原生参考。生产修改能否通过原有
完整场门槛，应以修正后的实际生产计算和完整检查器结果判断。
