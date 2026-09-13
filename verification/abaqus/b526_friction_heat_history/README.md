# B5.26 原默认摩擦生热参考归档

本目录逐字保存 2026-09-13 修正物理定义前的 B5.26 摩擦反转输入和四份原生 CSV。
五个文件的原始 SHA256 见 [SHA256SUMS](SHA256SUMS)。这些文件未重写数值、节点编号或时间。

原输入同时启用热机械接触和摩擦，省略 `*GAP HEAT GENERATION`。
按照 Abaqus 2018 的默认规则，全部摩擦耗散被转为热，并平均分给接触两侧。
Fuelsim 此例没有摩擦生热，因而原参考与 Fuelsim 的热方程不一致。
直接读取原生 `IVOL`、`TEMP`、`RFL11` 和 `ALLFD` 可验证该热源：全历史积分得到的
储热减边界热输入为 39.0203215444 J，摩擦耗散末值为 39.0203208923 J。

用修正后的逐点比较器检查历史代码 `e495bb8` 时，原参考存在 28 条热反力超限记录。
因此本目录是历史证据，不再作为当前同物理定义的验收参考。
原有误差门槛未改变；详细结果见
[摩擦生热诊断报告](../b526_diagnosis/friction_heat_diagnosis.md)。

上级目录中的当前输入仅增加 `*GAP HEAT GENERATION` 与 `0.0, 0.5` 两行。
当前四份 CSV 来自独立的完整诊断输入
[b526_friction_reversal_no_friction_heat.inp](../b526_diagnosis/b526_friction_reversal_no_friction_heat.inp)
的实际原生计算。原生输出到当前文件名的逐字映射与散列见
[active_reference_provenance.json](../b526_diagnosis/active_reference_provenance.json)。
