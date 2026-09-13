# B5.23 与 B5.44 正式参考来源同步

2026-09-13 运行统一测试时发现，两个正式 Abaqus 输入及其生成脚本已经更新，
但 SHA256 清单和用于比较的正式 CSV 仍对应此前输入。本次同步使用已有的原生
Abaqus 结果，没有重新运行 Abaqus，没有修改生产输入、数值实现或验收门槛。

| 正式参考前缀 | 原生运行前缀 | 输入 SHA256 |
|---|---|---|
| `b523_hex8_c3d8t_integrated_path` | `b523_diagnosis/b523_integrated_no_friction_heat` | `a19d5c87517177511192bf123d6628b81aba6346aaafd915be46291c158a81f0` |
| `b544_hex8_c3d8rt_distorted_bending` | `b544_diagnosis/b544_extended_tables_fixed_hourglass` | `fec3052af7d4e4a5611759fc3961a14ad759a1a3556201ce821b3c2bdacc337a` |

B5.23 的当前输入保留摩擦，但明确禁止将摩擦耗散转化为热量，与 Fuelsim 的
物理范围一致。B5.44 的当前输入保留初始单值密度，将原有材料直线延伸至
280 K 和 800 K，并显式采用初温对应的沙漏刚度 400000 Pa。

同步前逐项验证了当前输入与原生运行记录的 SHA256，四份原生 CSV 的 SHA256，
以及在临时目录中运行两个生成脚本所得输入与正式输入逐字一致。节点、材料点、
接触和能量 CSV 均从匹配的原生运行逐字复制；没有修改数值或筛选数据行。
正式前缀下的 `*_provenance.txt` 保留原始运行记录，包括原生模型名、原始输出
文件名及其散列。因此其中的 `diagnostic only` 描述原运行目的；本文件记录将
这两组已匹配输入的结果用于正式比较的变更。

在当前生产输出上，两个匹配参考均通过各自现有的完整场比较。B5.44 的
[初始密度报告](b544_initial_density_validation.md)记录此前材料表及默认沙漏
设置的历史比较，其误差表不代表当前正式输入。该历史输入与原生结果保存在
`b544_initial_density_history/`；本次同步没有改写历史资料。

运行生产程序并按原有门槛比较：

```bash
ctest --test-dir build -j4 --output-on-failure \
  -R '^fuelsim_b(523_hex8_c3d8t_integrated|544_hex8_c3d8rt_distorted_bending)_abaqus_tests$'
ctest --test-dir build -j4 --output-on-failure -R '^fuelsim_abaqus_reference_sha256$'
```
