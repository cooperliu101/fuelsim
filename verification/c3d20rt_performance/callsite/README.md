# 形函数调用边界整理与复验

日期：2026-09-16。本次基线为已推送的 `ccc0ca0`，已包含上一轮接触性能优化。

## 修改内容

`quad8_shape` 只计算自动微分形函数，`double_quad8_shape` 只计算双精度形函数。
原来在 `quad8_shape` 内将六组双精度数组转写成 `adlite::Scalar` 数组的分支已删除。

投影迭代调用处根据两个局部坐标是否携带导数选择计算入口。
双精度形函数及一、二阶导数直接传给接受 `std::array<double, 8>` 的具体插值函数。
局部坐标没有导数时，面节点坐标仍可能携带导数，因此插值保留节点的自动微分类型，
让双精度权重直接乘节点坐标，继续传播节点导数。

投影结果需要保存自动微分形函数权重的调用处直接使用自动微分入口。
没有增加模板、通用标量层或新的存储转换函数。

## 同批次性能比较

使用上一轮的六单元有限应变接触首增量例题，固定 CPU 0、一个求解进程和单线程，
保留输入卡原有输出与 MUMPS 直接求解器。每版预热一次，再交替测量五次。
两个版本的编译器、依赖和 Release 优化配置相同，测量期间不运行编译或测试。

| 指标 | 整理前 | 整理后 |
|---|---:|---:|
| 外部总耗时中位数，秒 | 4.337302 | 4.195413 |
| 程序内部总耗时中位数，秒 | 3.352598 | 3.199533 |
| 接受的时间步数 | 1 | 1 |
| 非线性迭代次数 | 14 | 14 |
| 残量计算次数 | 22 | 22 |
| 导数矩阵计算次数 | 14 | 14 |

本批次外部总耗时减少约 3.27%。这与上一轮不同批次的时间分别报告，不能将两个百分比直接相加。

```bash
python benchmarks/run_c3d20rt_comparison.py \
  /path/to/ccc0ca0/fuelsim build/fuelsim /tmp/c3d20rt-callsite-comparison \
  --case finite_contact_first_increment --cpu 0 --repeats 5
```

证据：[逐次测量](runs.csv)、[输入校验值](inputs.sha256)、[二进制校验值](executables.txt)。

读取两个版本的完整 Exodus 输出，比较时间、节点及单元字段的 1608 个数组，
所有有限数值差异均为零，数组形状及 NaN 位置相同。
具体结果见[逐字段比较](field-comparison.csv)。

局部形函数与曲面接触测试通过。完整 Release 构建后，运行统一回归：

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  taskset -c 0-7 ctest --test-dir build -j8 --output-on-failure
```

326 项测试全部通过，总耗时 150.91 秒；原有误差门槛和输入卡保持不变。
