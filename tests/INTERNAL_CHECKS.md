# 内部检查的合并入口

迁移时拆出的内部检查与二十节点非匹配接触检查，共 30 项，现在归入三个
CTest 入口：

| 入口 | 原检查数量 |
| --- | ---: |
| `fuelsim_internal_rz_tests` | 4 |
| `fuelsim_internal_hex8_tests` | 16 |
| `fuelsim_internal_hex20_tests` | 10 |

这里合并的是测试调度入口，不是删减断言，也不是减少生产输入卡。
各检查仍使用原有内部辅助程序、命令参数和独立进程，避免共享全局初始化状态。
同一组内顺序执行，某项失败或超时后仍继续执行后续检查，最后汇总失败名称。
不同组可以与其他测试并行调度；每组继承所有成员所需的生产结果依赖。

配置时生成的 `build/internal_checks/membership.tsv` 给出完整的组与检查对应关系，
同目录下的脚本列出实际执行命令。这些文件仅包含调度信息，不生成输入卡。
原有数值核综合测试暂不改变；本次合并针对迁移造成的细碎检查项。

统一回归仍使用 `ctest --test-dir build -j4 --output-on-failure`。
仅检查这些内部组时，可以使用下面的命令；CTest 会自动补入所需的生产测试：

```bash
ctest --test-dir build -j4 -R '^fuelsim_internal_' --output-on-failure
```

新增同类检查应使用 `fuelsim_add_internal_check` 加入已有组，不再单独注册 CTest。
检查源代码中的失败信息应指明物理分支或断言，方便从组日志定位问题。
