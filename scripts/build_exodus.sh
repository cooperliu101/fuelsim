#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <seacas-source-dir> <install-prefix> <dependency-prefix>" >&2
    exit 2
fi

seacas_source_dir=$1
exodus_install_prefix=$2
dependency_prefix=$3
build_dir="${seacas_source_dir}/build-fuelsim"

if [[ ! -f "${seacas_source_dir}/CMakeLists.txt" ]]; then
    echo "SEACAS CMake project not found: ${seacas_source_dir}" >&2
    exit 2
fi

env PATH="${dependency_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake -S "${seacas_source_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${exodus_install_prefix}" \
    -DCMAKE_C_COMPILER="${dependency_prefix}/bin/cc" \
    -DCMAKE_CXX_COMPILER="${dependency_prefix}/bin/c++" \
    -DBUILD_SHARED_LIBS=ON \
    -DSeacas_ENABLE_ALL_PACKAGES=OFF \
    -DSeacas_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
    -DSeacas_ENABLE_SEACASExodus=ON \
    -DSeacas_ENABLE_TESTS=OFF \
    -DTPL_ENABLE_DLlib=OFF \
    -DTPL_ENABLE_MPI=OFF \
    -DTPL_ENABLE_Netcdf=ON \
    -DNetCDF_ROOT="${dependency_prefix}" \
    -DSEACASExodus_ENABLE_MPI=OFF \
    -DSEACASExodus_ENABLE_Pnetcdf=OFF

env PATH="${dependency_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake --build "${build_dir}" --parallel 4

env PATH="${dependency_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    cmake --install "${build_dir}"
