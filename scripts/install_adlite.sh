#!/usr/bin/env bash

set -euo pipefail

expected_commit=a0e75a887017135d5e520062bb330fcb194b4399

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <adlite-source-dir> <install-prefix> <build-dir> <toolchain-prefix>" >&2
    exit 2
fi

adlite_source_dir=$1
adlite_install_prefix=$2
build_dir=$3
toolchain_prefix=$4

if [[ ! -f "${adlite_source_dir}/CMakeLists.txt" ]]; then
    echo "ADlite CMake project not found: ${adlite_source_dir}" >&2
    exit 2
fi
if [[ ! -x "${toolchain_prefix}/bin/c++" ]]; then
    echo "C++ compiler not found: ${toolchain_prefix}/bin/c++" >&2
    exit 2
fi
if ! git -C "${adlite_source_dir}" rev-parse --is-inside-work-tree \
    >/dev/null 2>&1; then
    echo "ADlite source must be a Git worktree: ${adlite_source_dir}" >&2
    exit 2
fi

actual_commit=$(git -C "${adlite_source_dir}" rev-parse HEAD)
if [[ "${actual_commit}" != "${expected_commit}" ]]; then
    echo "ADlite commit mismatch: expected ${expected_commit}, got ${actual_commit}" >&2
    exit 2
fi
if [[ -n $(git -C "${adlite_source_dir}" status --porcelain) ]]; then
    echo "ADlite source worktree must be clean: ${adlite_source_dir}" >&2
    exit 2
fi

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake -S "${adlite_source_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="${toolchain_prefix}/bin/c++" \
    -DCMAKE_INSTALL_PREFIX="${adlite_install_prefix}" \
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
    -DBUILD_TESTING=ON \
    -DADLITE_BUILD_EXAMPLES=OFF

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake --build "${build_dir}" --parallel 4

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    ctest --test-dir "${build_dir}" --output-on-failure

env PATH="${toolchain_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake --install "${build_dir}"

config_file="${adlite_install_prefix}/lib/cmake/adlite/adliteConfig.cmake"
if [[ ! -f "${config_file}" ]]; then
    echo "ADlite package configuration was not installed: ${config_file}" >&2
    exit 1
fi
