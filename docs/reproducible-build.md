# 可复现构建与持续集成

本文给出在新的 Linux 机器上重建 fuelsim 验收环境的完整步骤。所有依赖安装到
显式、持久的前缀；构建目录与安装目录分离。命令假定当前目录是 fuelsim 仓库
根目录。

## 固定版本

```text
MOOSE development environment: 2026.06.16, MPICH build
PETSc supplied by MOOSE:        3.25.2
ADlite version:                 0.2.3
ADlite commit:                  fd319e00234e18280319d141f17d9fa015c2501b
SEACAS Exodus tag:              v2024-06-27
```

`dependencies/moose-2026.06.16-linux-64.yml` 固定 MOOSE 开发包及直接使用的
构建工具版本。该环境文件只适用于 Linux x86-64。先安装 Conda 或兼容实现，
然后执行：

```bash
conda env create --file dependencies/moose-2026.06.16-linux-64.yml
conda activate fuelsim-moose-2026.06.16
```

选择持久目录。下面的示例把依赖源、构建树和安装前缀放在 fuelsim 仓库的同级
目录，不依赖系统临时目录：

```bash
fuelsim_source_root="$(cd .. && pwd)/fuelsim-dependency-sources"
fuelsim_build_root="$(cd .. && pwd)/fuelsim-dependency-builds"
fuelsim_dependency_root="$(cd .. && pwd)/fuelsim-dependencies"
fuelsim_toolchain_prefix="${CONDA_PREFIX}"
mkdir -p "${fuelsim_source_root}" "${fuelsim_build_root}" \
  "${fuelsim_dependency_root}"
```

## 安装 ADlite

从 ADlite 仓库取得源码，并检出下面的固定提交。安装脚本要求工作树没有修改，
并严格拒绝不是固定提交的源码：

```bash
git clone git@github.com:cooperliu101/ADlite.git "${fuelsim_source_root}/ADlite"
git -C "${fuelsim_source_root}/ADlite" checkout \
  fd319e00234e18280319d141f17d9fa015c2501b

./scripts/install_adlite.sh \
  "${fuelsim_source_root}/ADlite" \
  "${fuelsim_dependency_root}/adlite-0.2.3" \
  "${fuelsim_build_root}/adlite-0.2.3" \
  "${fuelsim_toolchain_prefix}"
```

脚本配置 Release ADlite 并启用链接期跨翻译单元优化，关闭示例、编译并运行
ADlite CTest，然后安装 CMake 软件包配置。Fuelsim 的 Release 构建默认使用同一
优化；两者必须同时启用才能复现 B5.46 的单核计时。提交不匹配、源码不干净、
测试失败或安装配置缺失都会返回非零。

## 安装串行 Exodus

PETSc 只负责求解，Exodus 作为独立串行输入输出库构建。脚本要求 SEACAS
源码精确位于 `v2024-06-27` 标签且工作树干净：

```bash
git clone --branch v2024-06-27 --depth 1 \
  https://github.com/gsjaardema/seacas.git \
  "${fuelsim_source_root}/seacas-2024-06-27"

./scripts/build_exodus.sh \
  "${fuelsim_source_root}/seacas-2024-06-27" \
  "${fuelsim_dependency_root}/exodus-2024-06-27" \
  "${fuelsim_build_root}/exodus-2024-06-27" \
  "${fuelsim_toolchain_prefix}"
```

该配置只启用 SEACASExodus，关闭 MPI、测试和其他 SEACAS 软件包，并复用固定
MOOSE 环境的 NetCDF。

## Release 验收

```bash
env \
  PATH="${fuelsim_toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
  PKG_CONFIG_PATH="${fuelsim_toolchain_prefix}/lib/pkgconfig" \
  cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="${fuelsim_toolchain_prefix}/bin/c++" \
    -DCMAKE_PREFIX_PATH="${fuelsim_dependency_root}/adlite-0.2.3" \
    -DSEACASExodus_DIR="${fuelsim_dependency_root}/exodus-2024-06-27/lib/cmake/SEACASExodus" \
    -DFUELSIM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel 4
ctest --test-dir build -j8 --output-on-failure
```

已有构建目录升级依赖时，还需在配置命令中显式指定
`-Dadlite_DIR="${fuelsim_dependency_root}/adlite-0.2.3/lib/cmake/adlite"`，
避免 CMake 沿用缓存中的旧版软件包路径。

Fuelsim 的 Release 配置默认启用链接期跨翻译单元优化；可用
`-DFUELSIM_ENABLE_RELEASE_IPO=OFF` 显式关闭，但这种构建不用于性能验收。只有
CTest 全部通过才构成 Release 验收；仅配置或编译成功不构成数值验收。

## 检测器验收

检测器配置同时启用地址检测器和未定义行为检测器：

```bash
env \
  PATH="${fuelsim_toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
  PKG_CONFIG_PATH="${fuelsim_toolchain_prefix}/lib/pkgconfig" \
  cmake -S . -B build-sanitize \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER="${fuelsim_toolchain_prefix}/bin/c++" \
    -DCMAKE_PREFIX_PATH="${fuelsim_dependency_root}/adlite-0.2.3" \
    -DSEACASExodus_DIR="${fuelsim_dependency_root}/exodus-2024-06-27/lib/cmake/SEACASExodus" \
    -DFUELSIM_ENABLE_SANITIZERS=ON \
    -DFUELSIM_WARNINGS_AS_ERRORS=ON
cmake --build build-sanitize --parallel 4
ASAN_OPTIONS=detect_leaks=0 \
MPIR_CVAR_CH4_NETMOD=ofi \
FI_PROVIDER=tcp \
  ctest --test-dir build-sanitize -j8 --output-on-failure
```

泄漏检查关闭是显式边界：当前 MOOSE PETSc 初始化会加载系统 CUDA 驱动，
`libcuda.so.1` 在进程退出时保留分配，无法归因于 fuelsim。该入口仍完整执行
地址越界和未定义行为检查，不能据此声称 LeakSanitizer 已通过。检测器任务和
Release 任务运行同一个轻量 CTest 套件，不再通过标签形成两套回归。检测器任务
还固定 MPICH 使用 OFI 网络模块和 TCP provider；默认 UCX 网络模块在
`PetscInitialize` 进入 fuelsim 代码前的地址交换中会触发外部库越界读取。
Release 任务仍使用环境默认网络模块，且所有双进程等价性测试都必须通过。

## 持续集成执行器

`.github/workflows/source-audit.yml` 保留 GitHub 托管执行器上的源码和参考文件
哈希审计，并增加两个专用 Linux 自托管任务。执行器必须带有
`self-hosted, linux, x64, fuelsim` 标签，并配置两个仓库变量：

```text
FUELSIM_TOOLCHAIN_PREFIX = 固定 MOOSE Conda 环境的绝对路径
FUELSIM_DEPENDENCY_ROOT  = 包含 adlite-0.2.3 和 exodus-2024-06-27 的绝对路径
```

Release 与检测器任务分别从空的配置目录重新运行 CMake、全量编译和全部 CTest；
两者都启用编译警告即错误。检测器任务采用上面相同的泄漏和 MPI 网络模块边界。
自托管任务不运行 MOOSE 生成参考结果；它只消费已经追踪且通过 SHA256 审计
的快照。
