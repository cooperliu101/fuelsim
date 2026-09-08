# C3D20RT 无摩擦算例精简

完整有摩擦算例 `../fuelsim/transient_c3d20rt_finite_contact_safeguarded.fsi`
保持七秒、280 步，作为完整往复加载的手动外部验证。现有小应变有摩擦
自动测试也保留，它验证另一种应变形式。

无摩擦生产输入 `../fuelsim/transient_c3d20rt_finite_contact_frictionless.fsi`
缩短为两秒、80 步。0 至 1 秒压紧至 ux=-0.02 m，同时 uy 增至 0.005 m；
1 至 2 秒保持压紧位移，uy 继续增至 0.01 m。步长仍为 0.025 s，几何、
材料、边界、罚刚度、零摩擦系数和热接触参数均不变。步数减少 71.43%，
不据此声称相同比例的运行时间提升。

原七秒输入逐字保存在 `c3d20rt_frictionless_full_history.fsi`，仍可手动运行。
原 280 时刻 Abaqus CSV 和历史结果全部保留；三个 `c3d20rt_finite_contact_frictionless_short_*.csv`
文件分别为原节点、接触和材料积分点参考的前 80 个时刻，保留原始字节和时间值，
未重新生成参考值。比较范围为 0.025 至 2 s，无需重新运行 Abaqus。

自动测试 `fuelsim_c3d20rt_finite_frictionless_abaqus_tests` 使用逐字复制的完整
输入卡运行实际 `fuelsim -i`，并检查全部 80 步、88 个节点、13 个接触节点和
每单元八个材料积分点。它单独检查摩擦系数为零时切向力严格为零，并要求
1 秒之后在受压接触下存在几何滑移，以及非零法向接触力和接触传热。
所有原有温度、位移、反力、间隙、压力、接触力、滑移、热流和应力对标仍保留。

无摩擦容差沿用原有批准范围：反力、滑移最大逐点相对误差小于 1%，
零参考反力绝对误差小于 0.1 N；整体相对误差及其他非零场仍小于 0.1%，
其他零参考绝对门槛不变。有摩擦专用容差不用于此测试。

复现：

```bash
cmake --build build --parallel 4
ctest --test-dir build -j4 --output-on-failure -R '^fuelsim_c3d20rt_finite_frictionless_abaqus_tests$'
ctest --test-dir build -j4 --output-on-failure
```

验证结果：新增无摩擦生产测试通过，完整有摩擦已有结果重新检查也通过。
短算例切向接触力全程严格为零，最大逐点相对误差：反力 0.00144982%，
滑移 0.00264361%，压力 0.0000342637%，接触总热流 0.00000201246%。
零参考反力最大绝对差为 0.000418130 N。证据见
`c3d20rt_frictionless_short_metrics.txt`、`c3d20rt_frictionless_short_ctest.txt`
和 `c3d20rt_friction_retained_metrics.txt`。

统一回归运行了 304 项测试，用时 188.69 s；新增无摩擦测试在其中通过。
首次结果为 299/304：四个 `fuelsim_b13_small_cax4t/cax4rt/cax8t/cax8rt_abaqus_tests`
未通过，另一个失败是本次修改文档后尚未刷新 SHA256 清单。文档清单随后更新并
单独复查通过，见 `c3d20rt_frictionless_short_manifest_recheck.txt`。四个 B13 小应变失败仍保留，不声称全套回归通过；本次没有修改
这些算例，也未建立干净基线来重新归因其失败。完整日志见
`c3d20rt_frictionless_short_full_ctest.txt`。
