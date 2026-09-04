# Fuelsim 测试入口

测试分成两个边界清楚的类别。

用户可运行的端到端算例只执行生产程序：

```bash
fuelsim -i verification/fuelsim/<case>.fsi
```

每张 `.fsi` 输入卡都是版本库中完整、可直接阅读和手工运行的文件。CTest 只把
输入卡逐字复制到 `build/blackbox/<case>/verification/fuelsim/`，并为每个算例
提供独立输出目录，因此这些测试可以通过 `ctest -j4` 并行运行。生产程序退出后，
`fuelsim_production_result_tests` 只读取 Exodus、CSV 和外部参考结果；它不链接
`fuelsim_core`、`fuelsim_io` 或 `fuelsim_solver`，也不会再次求解。

局部自动微分 Jacobian、材料状态提交与回滚、接触候选所有权、网格映射以及
PETSc 并行贡献所有权等不能仅由最终结果文件证明的契约，继续由普通 C++ 内部
测试覆盖。这些程序不是用户算例入口，也不得代替生产程序的端到端测试。

新增端到端算例时，应当提交显式输入卡，为输出指定算例专用文件名，并通过
`run_fuelsim_case.cmake` 调用生产程序和只读结果检查器。不得在 C++、CMake 或
Shell 中拼装输入卡正文。
