# 工程验证矩阵与适用边界

fuelsim 的“已验证”只表示仓库中存在可重复执行的测试、明确的判据和追踪的
证据文件，不表示完成核安全软件鉴定。机器可读的唯一索引是
`verification/verification_matrix.tsv`，`fuelsim_verification_matrix_tests`
会检查矩阵模式、必需能力行、CTest 名称以及每个证据文件是否存在且非空。

## 证据等级

- `verified`：有自动 CTest，并且验收判据和参考文件均进入 Git；
- `qualified`：有自动 CTest，但某项门槛尚未达到项目目标；必须明确记录例外，
  不能写成完全验证；
- `measured`：按受控环境运行的性能观测，不作为每次 CTest，也不外推为缩放
  结论；
- `limitation`：明确未实现或未鉴定的边界，是发布证据的一部分，不能省略。

所有 fuelsim-to-MOOSE 全场误差使用相对 L2、相对绝对峰值和最大逐点相对
误差。参考值为零的点不构造相对误差，而是另外报告数量和最大绝对差。矩阵中
M0、M1 和 M3.3 的门槛为 `1%`；M2.1、M2.2 和 M3.1 的门槛为 `0.1%`。
M2.3 当前是 `qualified`：温度、轴向位移和平均非弹性状态使用 `0.1%`，径向
位移 L2、峰值和逐点门槛分别为 `0.095%`、`0.125%` 和 `0.27%`，接触压力
为 `0.11%`；积分点应力、塑性和蠕变应变按各自当前误差采用
`0.015%`--`0.305%` 的分项
门槛。具体实际数值由对应 CTest 输出，不能只凭构建成功宣称通过。
M4.1 有限应变 PCMI 使用统一 `0.5%` 门槛；温度、径向/轴向位移、接触压力，
以及逐积分点等效应力、等效塑性应变和等效蠕变应变的三项误差必须全部通过。

## 发布验收

规定环境中的候选提交至少执行：

```bash
env \
  PATH=/home/cooper/miniforge/envs/moose/bin:/usr/local/bin:/usr/bin:/bin \
  PKG_CONFIG_PATH=/home/cooper/miniforge/envs/moose/lib/pkgconfig \
  cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/home/cooper/miniforge/envs/moose/bin/c++ \
  -DCMAKE_PREFIX_PATH=/tmp/adlite-fuelsim-install \
  -DSEACASExodus_DIR=/home/cooper/.local/exodus-2024-06-27/lib/cmake/SEACASExodus
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

涉及性能的候选还必须按 `benchmarks/README.md` 固定 CPU 和全部线程环境变量，
完成 1,584 DOF 配对计时以及 23,010 DOF、20 步工况。MPI 沙盒初始化失败需在
沙盒外复现；不能跳过并行测试后仍声明发布验收通过。

## 当前工程适用边界

当前版本是燃料棒性能分析的基础有限元求解器，适合继续承载经过独立验证的
物理模型，但不是完整的燃料性能产品。边界包括：

- 仅支持 2D 轴对称 RZ、Quad4 体单元和 Line2 边；区域可选小应变或
  Eulerian Hencky 有限应变，没有 3D 或非轴对称载荷；
- 材料只有线性温度斜率热弹性、通用 J2 线性硬化塑性、Norton 蠕变及其
  全隐式耦合；斜率采用 AD 活跃求值，但尚无 UO2/包壳温度、燃耗、辐照
  相关的工程关联；
- 当前有限应变塑性—蠕变证据仅覆盖关闭有限应变历史旋转的近共轴 PCMI；
  历史张量没有随有限转动客观旋转，压力和牵引仍为参考构形 dead load，
  因而一般非共轴大转动与 follower load 尚未实现或鉴定；
- 尚无裂变气体释放、棒内气压与气体组分演化、芯块致密化/肿胀/开裂/重定位、
  包壳辐照生长/肿胀、氧化腐蚀、冷却剂通道或功率—燃耗历史模型；
- 接触为无摩擦罚接触与气隙导热；没有摩擦、粘结、脱粘、接触热阻经验关联或
  压力相关导热、辐射或燃料碎块之间的三维运动；当前候选面只适合小滑移，
  越界会拒步而不是动态搜索全部主面；
- 已对代表性稳态热问题、Backward Euler 热瞬态和罚接触分别自动检查空间
  网格、时间步和罚参数收敛；J2 独立证据仍以单调加载和卸载为主，每个新增
  工程工况仍必须执行自己的收敛研究；
- 时间积分只有一阶 Backward Euler；可选 step-doubling 控制温度/位移场的
  局部时间误差，但多年蠕变和强历史依赖工况仍须独立做全局时间步收敛；
- 分场残量有独立诊断和复核；基于初始残量的自动行缩放仅为可选功能且默认
  关闭，尚无由用户指定且跨载荷步固定的物理特征尺度；
- 多 rank 已分布贡献装配和 PETSc 线性代数，但网格、问题几何和回调完整状态
  仍复制，Exodus 仍由各 rank 串行读取，因此内存扩展性尚未完成；
- 检查点只保证当前严格格式和相同问题签名，不提供跨版本迁移；
- 当前验证工况是小型或中型基准，不能替代独立代码审查、网格/时间步收敛、
  试验数据验证、软件质量保证和法规鉴定。

新增工程物理时，必须先在矩阵中增加 `verified` 行或扩展现有行，绑定输入卡、
参考网格/结果、局部 AD 切线、解析或独立代码对标和明确门槛；在证据齐全前应
标为 `qualified` 或 `limitation`，不能用“计划支持”描述为当前能力。
