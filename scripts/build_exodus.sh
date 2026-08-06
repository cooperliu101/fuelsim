#!/usr/bin/env bash

set -euo pipefail

expected_tag=v2024-06-27

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <seacas-source-dir> <install-prefix> <build-dir> <toolchain-prefix>" >&2
    exit 2
fi

seacas_source_dir=$1
exodus_install_prefix=$2
build_dir=$3
toolchain_prefix=$4

if [[ ! -f "${seacas_source_dir}/CMakeLists.txt" ]]; then
    echo "SEACAS CMake project not found: ${seacas_source_dir}" >&2
    exit 2
fi
if [[ ! -x "${toolchain_prefix}/bin/cc" ||
      ! -x "${toolchain_prefix}/bin/c++" ]]; then
    echo "C/C++ compilers not found below: ${toolchain_prefix}/bin" >&2
    exit 2
fi
if ! git -C "${seacas_source_dir}" rev-parse --is-inside-work-tree \
    >/dev/null 2>&1; then
    echo "SEACAS source must be a Git worktree: ${seacas_source_dir}" >&2
    exit 2
fi

actual_tag=$(git -C "${seacas_source_dir}" describe --tags --exact-match HEAD \
    2>/dev/null || true)
if [[ "${actual_tag}" != "${expected_tag}" ]]; then
    echo "SEACAS tag mismatch: expected ${expected_tag}, got ${actual_tag:-none}" >&2
    exit 2
fi
if [[ -n $(git -C "${seacas_source_dir}" status --porcelain) ]]; then
    echo "SEACAS source worktree must be clean: ${seacas_source_dir}" >&2
    exit 2
fi

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake -S "${seacas_source_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${exodus_install_prefix}" \
    -DCMAKE_C_COMPILER="${toolchain_prefix}/bin/cc" \
    -DCMAKE_CXX_COMPILER="${toolchain_prefix}/bin/c++" \
    -DBUILD_SHARED_LIBS=ON \
    -DSeacas_ENABLE_ALL_PACKAGES=OFF \
    -DSeacas_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
    -DSeacas_ENABLE_SEACASExodus=ON \
    -DSeacas_ENABLE_TESTS=OFF \
    -DTPL_ENABLE_DLlib=OFF \
    -DTPL_ENABLE_MPI=OFF \
    -DTPL_ENABLE_Netcdf=ON \
    -DNetCDF_ROOT="${toolchain_prefix}" \
    -DSEACASExodus_ENABLE_MPI=OFF \
    -DSEACASExodus_ENABLE_Pnetcdf=OFF

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake --build "${build_dir}" --parallel 4

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake --install "${build_dir}"

config_file="${exodus_install_prefix}/lib/cmake/SEACASExodus/SEACASExodusConfig.cmake"
if [[ ! -f "${config_file}" ]]; then
    echo "Exodus package configuration was not installed: ${config_file}" >&2
    exit 1
fi
