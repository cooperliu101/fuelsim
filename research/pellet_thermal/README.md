# 三维芯块热传导代理超单元研究

本研究的实现、算例、训练脚本、模型和验证记录统一保存在本目录。生产目录只保留
输入解析、纯热区域装配、自由度分配和结果输出所需的接入修改。

- `include/`、`src/`：精确凝聚、原生 C++ 网络推理及完整有限元采样程序。
- `cases/`：受版本控制的 Exodus 网格和完整生产输入卡。
- `python/`：网格生成、PyTorch 训练和权重导出。
- `models/`：可复现回归所需的网络权重及训练统计。
- `tests/`：四项核心验证，端到端检查只运行生产程序和读取生产输出。
- `docs/`：方法、适用范围和验证结果。
- `runs/`、`datasets/`：忽略的本地运行产物与训练数据。

实现细节和边界见 [方法说明](docs/method.md)，具体验证结果见 [验证记录](docs/validation.md)。
新增的 50 μm 偏心气隙工况见 [偏心验证](docs/eccentric.md)，使用原有权重比较八个周向区间的非均匀热率。

切线直接比较、物理恒等式检查和八组参数组合见 [补充验证](docs/qualification.md)。

## 构建与回归

从 fuelsim 仓库根目录运行。沿用已经配置好的 `build` 目录以及原 `moose` Conda PETSc、
ADlite 和独立串行 Exodus。Python 需要 NumPy、netCDF4 和 PyTorch；只有重新生成网格需要
SciPy。C++ 推理不依赖 Python。

```bash
cmake --build build --target fuelsim fuelsim_pellet_exact_tests fuelsim_pellet_response_tests fuelsim_pellet_infer fuelsim_pellet_sample --parallel 8
ctest --test-dir build -R '^fuelsim_pellet_' --output-on-failure -j 2
```

CTest 使用已导出的模型，不在每次回归时重新训练。`FUELSIM_SECTION_PYTHON` 指向具有上述
Python 软件包的解释器。本次实测为 `/home/cooper/miniforge/bin/python`，PyTorch 运算全部在 CPU。
运行 PETSc 的环境必须允许其 MPI 本地通信。

## 重新采样和训练

```bash
mkdir -p research/pellet_thermal/datasets
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/research/pellet_thermal/fuelsim_pellet_sample research/pellet_thermal/cases/pellet_full.fsi research/pellet_thermal/datasets/samples.txt 1024
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/cooper/miniforge/bin/python research/pellet_thermal/python/train.py --data research/pellet_thermal/datasets/samples.txt --model research/pellet_thermal/models/pellet_mlp.txt
ctest --test-dir build -R '^fuelsim_pellet_' --output-on-failure -j 2
```

采样程序逐个求解完整三维芯块；没有使用精确凝聚或网络产生训练标签。网格和材料签名、
数据 SHA-256、随机种子以及训练和保留样本误差保存在模型及训练记录中。

## 手动运行包壳耦合算例

```bash
mkdir -p research/pellet_thermal/runs/manual
cp -R research/pellet_thermal/cases research/pellet_thermal/models research/pellet_thermal/runs/manual/
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/fuelsim -i research/pellet_thermal/runs/manual/cases/coupled_full.fsi
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/fuelsim -i research/pellet_thermal/runs/manual/cases/coupled_exact.fsi
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/fuelsim -i research/pellet_thermal/runs/manual/cases/coupled_surrogate.fsi
```

三个输入卡的芯块网格、材料、热源、间隙和包壳相同，只改变芯块响应选择及输出文件名。
研究的网络准确性只覆盖所提交的小网格和常数导热材料，尚未鉴定其他网格、温度相关导热
或整体运行时间收益。
