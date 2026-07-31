#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <petsc-source-dir> <install-prefix> <dependency-prefix>" >&2
    exit 2
fi

petsc_source_dir=$1
petsc_install_prefix=$2
petsc_dependency_prefix=$3
petsc_arch=arch-fuelsim-exodus-opt

if [[ ! -x "${petsc_source_dir}/configure" ]]; then
    echo "PETSc configure script not found: ${petsc_source_dir}/configure" >&2
    exit 2
fi

env PATH="${petsc_dependency_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    "${petsc_source_dir}/configure" \
    PETSC_ARCH="${petsc_arch}" \
    --prefix="${petsc_install_prefix}" \
    --with-cc="${petsc_dependency_prefix}/bin/mpicc" \
    --with-cxx="${petsc_dependency_prefix}/bin/mpicxx" \
    --with-fc=0 \
    --with-cuda=0 \
    --with-x=0 \
    --with-debugging=0 \
    COPTFLAGS=-O3 \
    CXXOPTFLAGS=-O3 \
    --with-shared-libraries=1 \
    --with-blaslapack-dir="${petsc_dependency_prefix}" \
    --with-hdf5-dir="${petsc_dependency_prefix}" \
    --with-netcdf-dir="${petsc_dependency_prefix}" \
    --with-zlib-dir="${petsc_dependency_prefix}" \
    --download-pnetcdf \
    --download-exodusii

env PATH="${petsc_dependency_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    make -C "${petsc_source_dir}" --jobs=4 \
    PETSC_DIR="${petsc_source_dir}" \
    PETSC_ARCH="${petsc_arch}" \
    all

env PATH="${petsc_dependency_prefix}/bin:/usr/local/bin:/usr/bin:/bin" \
    make -C "${petsc_source_dir}" \
    PETSC_DIR="${petsc_source_dir}" \
    PETSC_ARCH="${petsc_arch}" \
    install
