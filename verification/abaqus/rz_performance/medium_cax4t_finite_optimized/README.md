# CAX4T 代码优化与精度、速度验证

本次仅减少单元调用中未使用的计算，没有修改物理输入、材料公式、积分规则、
20 个载荷增量、MUMPS、求解容差或预测设置。原版为提交 `011636f` 对应的
生产可执行文件，在修改代码前直接保存；它的 SHA256 与此前正式性能比较一致。
正式测量前已移除临时计时代码。

## 结果

两种程序均采用单进程、单线程、逻辑 CPU 0 和关闭结果输出的完整输入。
Fuelsim 两个版本分别预热一次，再依次执行“原版、优化版、优化版、原版”。
完成有限应变和小应变计时后，才运行 Abaqus 的一次预热和两次正式测量。
没有同时运行其他求解、编译或回归任务。

| CAX4T 执行方式 | 原版正式两次外部耗时 | 原版均值 | 优化版正式两次外部耗时 | 优化版均值 | 耗时减少 |
|---|---:|---:|---:|---:|---:|
| 有限应变，逐增量提交历史 | 45.3800 / 45.0806 s | 45.2303 s | 35.8703 / 36.0972 s | 35.9838 s | **20.44%** |
| 小应变，稳态加载 | 21.9705 / 21.6687 s | 21.8196 s | 21.7320 / 21.6228 s | 21.6774 s | 0.65% |

小应变差异较小，视为基本持平，不作统计显著性结论。
全部 Fuelsim 运行均完成 20 个增量或载荷步、46 次非线性迭代，有限应变没有
失败重试。两种版本的最终区域及接触工程输出逐项一致。

同批 Abaqus 有限应变正式外部耗时为 43.7553、43.6932 s，均值 **43.7243 s**；
优化版 Fuelsim 用时少 **17.70%**，Abaqus/Fuelsim 耗时比约 **1.215**。
Abaqus 均为 20 个增量和 53 次包含接触不连续迭代的总迭代。
其 `JOB TIME SUMMARY` 分析时间两次均为 40 s。
Fuelsim 的内部 `total_seconds` 与各次外部计时独立保存在 `timing.json`，
不以内部时间代替外部程序总耗时。

这些是同一台 i9-13980HX 主机上 WSL2 Fuelsim 与 Windows Abaqus 的观测。
虚拟机内外逻辑 CPU 0 不足以证明物理核映射相同，两个正式重复也不支持推广成
所有硬件及所有 CAX4T 工作负载的通用速度比。

## 定位与代码变化

一次原版内部定位运行记录如下，单位为秒。这些值用于定位，不作为正式外部
计时数据；未请求切线的调用包含残量、诊断和历史更新等路径。

| 单元调用 | 调用次数 | 几何和平均量 | 本构、应力及其导数处理 | 完整历史输出 | 热力残量和切线装配 |
|---|---:|---:|---:|---:|---:|
| 未请求切线 | 786,944 | 3.0717 | 6.1866 | 5.2581 | 2.0659 |
| 请求切线 | 341,504 | 4.7885 | 10.1924 | 2.5018 | 4.8999 |

原实现即使最终丢弃历史输出，也会额外调用一次材料响应并旋转、整理试探历史。
完整历史输出阶段合计约 7.76 s，其中接受增量时需要的计算必须保留。
另有用于跨积分点闭式导数的系数张量，其转动结果只读取数值，却传播了随后
丢弃的自动微分导数。

修改集中在以下位置：

1. `src/core/cax4_evaluation.hpp` 向局部库传递已有的 `ElementRequest`。
   求残量和切线时不请求历史；应力输出单独请求应力；接受增量仍请求完整历史。
   `src/core/spatial_problem.cpp` 的接受路径明确保留历史请求。
2. `elements/src/cax4t.cpp` 仅在请求历史时执行额外的完整历史输出计算。
   实际材料应力和本构切线计算仍然保留，试探历史只在全局接受路径提交。
3. 跨积分点链中的应力修正系数只使用转动的普通数值；完整应力仍使用带导数
   的转动，因此几何切线保留。仅求残量时跳过这些切线系数的转动与累加。
4. 同一积分点、同一次调用中的已接受热应变直接复用，避免相同坐标、温度和
   时间上下文的重复材料函数调用。没有增加跨试探态或跨增量缓存。

仍使用宽度 6 的运动学、宽度 5 的本构自动微分和闭式节点链，没有扩大播种
宽度，也没有跳过必须的应力或历史事务计算。公共局部单元接口保持不变。

## 精度与回归

- 有限应变：全部 21 个输出时刻的 **228 个数值数组及 7 个名称等文本数组**
  与原版完全一致；只有创建信息和来源文字记录不参与比较。
- 小应变：最终输出的 **89 个数值数组及 7 个文本数组**与原版完全一致。
- 有限应变对 Abaqus 的全部最终场指标仍通过 **0.01%** 门槛，最大逐点
  相对误差仍为 **0.00000430214%**。小应变全场比较也通过。
- `elements/tests/cax4t_tests.cpp` 在已有小应变和有限应变塑性、蠕变共同激活
  检查中，增加输出请求组合的逐项等价性、空历史返回、独立应力输出和热量
  诊断检查；已有中心差分方向导数和失败后重复调用检查继续保留。
- Release 完整构建成功，`ctest --test-dir build -j4 --output-on-failure`
  **281/281 通过**，用时 172.48 s，包含历史事务、失败恢复和 MPI 等价性。
  没有新增中等规模自动测试，也没有改变任何已有误差门槛。

完整数值数组比较由 `check_equivalence.py` 读取两个生产 Exodus 结果执行。
相对 Abaqus 的误差仍由上级目录 `compare.py` 计算，并保留全部节点、29,696 个
活跃材料积分点和 65 个接触节点，零参考量单独按绝对误差检查。

## 重复运行

保留或按相同编译器、ADlite、PETSc、Exodus 和 Release LTO 配置构建提交
`011636f` 的原版可执行文件，再构建当前版本。使用全新的输出目录执行：

```bash
python verification/abaqus/rz_performance/medium_cax4t_finite_optimized/run_comparison.py \
  --baseline /absolute/path/to/baseline/fuelsim \
  --candidate "$PWD/build/fuelsim" --root "$PWD" --results /tmp/cax4t-comparison-repeat
```

该脚本只选择已提交的完整 `.fsi` 输入，不生成或修改物理定义。Abaqus 使用
上级目录 `run.ps1 -Size medium -Element cax4t -Strain finite -Timing -Runs 3`。
所有求解应依次执行，并避免与编译或回归同时运行。

数组等价性检查示例：

```bash
env NETCDF_LIBRARY=/home/cooper/miniforge/envs/moose/lib/libnetcdf.so \
  /home/cooper/miniforge/envs/moose/bin/python \
  verification/abaqus/rz_performance/medium_cax4t_finite_optimized/check_equivalence.py \
  /absolute/path/to/baseline_results.e /absolute/path/to/candidate_results.e \
  --report /tmp/cax4t-equivalence.json
```

`summary.json` 保存汇总，`timing.json` 保存原版与优化版的逐次时间，
`provenance.json` 保存源代码、输入、可执行文件和完整结果的 SHA256；
`comparison_provenance.json` 保存实际计时路径及摘要。
`profile.json` 和压缩日志保存定位证据，`ctest.log.gz` 保存完整回归记录。
精度和等价性报告分别为 `finite_accuracy.json`、`small_accuracy.json`、
`finite_equivalence.json`、`small_equivalence.json`。
