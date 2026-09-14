# 自动测试的 Abaqus 2025 参考更新

2026-09-14 已完成四项失败的后续独立诊断，见 [诊断报告](../four_case_diagnosis.md)。
下文“原因尚未鉴定”是首次更新时的状态；诊断阶段未改变正式输入和门槛。随后已按用户授权调整 B5.23/B11.4，
当前状态见 [两项分项容差记录](../two_case_tolerance_qualification.md)。

本次将自动测试仍在读取的 Abaqus 2018 数值结果改为 Windows Abaqus 2025
RELr427 的原生计算结果。145 项运行任务全部成功；一项任务可以包含多个原生作业。
原有输入的物理参数、网格和时间历史保持原样，输入字节同一性见
[input_identity.tsv](input_identity.tsv)。原先已更新为 2025 的数值结果继续保留。

日常 `ctest --test-dir build -j8 --output-on-failure` 不启动 Abaqus。
它运行 Fuelsim 生产程序或内部自检，再读取保存的参考文件。CTest 使用 8 个
并发额度，双进程 MPI 测试占两个额度，各进程的数学库线程数为 1。
编译继续使用 `cmake --build build --parallel 4`。

## 原生参考的生成与追溯

[native_tasks.tsv](native_tasks.tsv) 保存全部成功任务的运行入口、参数、开始和结束时间、
对应日志；[logs](logs/) 中保留原生版本与完成记录。最初的批量脚本把 Abaqus 的
许可证提示误当作 PowerShell 异常，另有 B7 运行命令引号错误，均修复后重新运行。
三个成功任务的状态因并发收集时序遗漏，修正收集逻辑后也重新运行并验证。
这些脚本问题与数值比较失败分别记录。

在 Windows PowerShell 中复现本批参考：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\verification\abaqus\2025\suite_refresh\run_all.ps1 -MaximumJobs 8
python .\verification\abaqus\2025\suite_refresh\select_frictionless_history.py
```

运行入口固定为 `C:\SIMULIA\Commands\abq2025.bat`。批量程序最多启动八个
独立单核原生作业，总计算额度为八核。此时应停止其他测试，以避免超过该额度。
C3D20RT 无摩擦自动例题仍只比较前 80 个时刻；选择程序从新生成的完整历史中
逐行保留 0.025 至 2 秒的数据，不修改结果数值，也不改变原有比较范围。

旧目录中的静态 Exodus 网格、验证规则表与历史 SHA256 清单继续使用。
这些文件不是 Abaqus 2018 的数值输出。当前参考路径见 [test_reference_paths.tsv](test_reference_paths.tsv)。
旧目录中的 2018 数值结果保留用于追溯，自动数值比较不再读取它们。

## 数值验证结果与边界

首次完整运行耗时 112.68 秒，309/314 项通过。四项为数值比较失败，另一项是
更新运行脚本后尚未刷新的参考文件 SHA256 清单检查；首次完整输出保存在
[initial_full_ctest.log](initial_full_ctest.log)。

更新清单及验证矩阵后再次完整运行，最终 **310/314 项通过，耗时 113.00 秒**。
参考 SHA256 清单检查与验证矩阵检查均通过，仍失败的只有下表四项。
最终完整日志见 [full_ctest.log](full_ctest.log)，逐项状态见
[test_results.tsv](test_results.tsv)。

本次新增 397 份 CSV 参考文件；183 项测试的命令参数显式使用 2025 数值参考。
本次没有修改生产算法、输入物理或误差门槛。此前用户批准的八项分项容差继续生效。
新参考暴露的四项差异不自动套用这些容差，验证矩阵相应标记为 `unverified`：

| 例题 | 未通过的量 | 当前证据 |
|---|---|---|
| B5.23，C3D8T 有限应变综合热机械接触 | 节点热反力 | 最大逐点相对误差约 0.72854%，原条件为小于 0.5% 或绝对误差小于 0.05 W；有五个样本同时超过这两个条件。其他正式验收量通过。 |
| B11.4，CAX8T 接触压力恢复 | 规定状态下的节点热反力绝对误差 | 1.43067424e-11 W，原门槛为 1e-11 W；接触压力、剪应力和接触力检查通过。 |
| B14，CAX4T 节点到面接触 | 界面节点热反力 | 最大逐点相对误差约 6.29443454e5（比值）。一个对应样本为 Fuelsim 4.21429818 W、Abaqus 6.69526621e-6 W。既有比较范围内的机械量通过。 |
| B14，CAX8T 节点到面接触 | 界面节点热反力 | 最大逐点相对误差约 6.29443454e5（比值）。一个对应样本为 Fuelsim 4.21429833 W、Abaqus 6.69526644e-6 W。既有比较范围内的机械量通过。 |

B14 两项的输入均与旧输入逐字相同，提取程序仍读取原生 RFL11 节点热反力。
本次确认了新旧原生响应不同，尚未鉴定造成这一差异的 Abaqus 内部机制。
B11.4 的热反力在极小绝对量上超过既有门槛，也仍如实保留失败状态。
这些结果不能表述为完整套件数值验收通过。
