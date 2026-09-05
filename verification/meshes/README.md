# 显式生产测试网格

- `b510_hex8_unit_cube.e` 保留原 B5.10 系列的八个源节点、单位立方体单元、
  `solid` 区域以及 `left`、`right`、`y0`、`z0` 边集，供十一条非弹性路径使用。

这里保存从原内部测试几何导出的 Exodus 网格。输入卡直接引用这些受版本管理
的文件，CTest 不生成网格或输入卡。网格导出只用于维护资料，不执行物理求解。

- b55_hex8_two_elements.e 保留原 B5.5 的十二个源节点、两个八节点六面体、
  solid 元素块和六组命名外表面；源节点编号与 Abaqus 参考一致。

- `b544_distorted_bending.e` 保留原 B5.44 的三十个节点、八个畸变单元和
  `left`、`right`、`y_high`、`z_high` 边集。
- `b59_coarse.e`、`b59_refined.e`、`b59_distorted.e` 保留原 B5.9 的三种
  两层板件网格。粗网格有十八个节点、四个单元；加密和畸变网格各有七十五个节点、
  三十二个单元。`lower` 和 `upper` 两块在材料界面共用源节点。

上述网格由 `benchmarks/test_mesh_export.cpp` 中对应的 `b544`、
`b59_coarse`、`b59_refined`、`b59_distorted` 分支手动导出。导出工具
保留原测试的坐标和连接关系，输入卡正文不由工具生成。
