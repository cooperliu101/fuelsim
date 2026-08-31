#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <fuelsim-source-root>" >&2
    exit 2
fi

source_root=$1

require_text() {
    local path=$1
    local expected=$2
    if ! grep -Fq -- "${expected}" "${path}"; then
        echo "[FAIL] ${path} does not contain required text: ${expected}" >&2
        exit 1
    fi
}

bash -n "${source_root}/scripts/install_adlite.sh"
bash -n "${source_root}/scripts/build_exodus.sh"

if [[ ! -x "${source_root}/scripts/install_adlite.sh" ||
      ! -x "${source_root}/scripts/build_exodus.sh" ]]; then
    echo "[FAIL] dependency installation scripts must be executable" >&2
    exit 1
fi

require_text "${source_root}/scripts/install_adlite.sh" \
    "a0e75a887017135d5e520062bb330fcb194b4399"
require_text "${source_root}/scripts/build_exodus.sh" "v2024-06-27"
require_text \
    "${source_root}/dependencies/moose-2026.06.16-linux-64.yml" \
    "moose-dev=2026.06.16=mpich"
require_text \
    "${source_root}/dependencies/moose-2026.06.16-linux-64.yml" \
    "numpy=2.5.1"
require_text \
    "${source_root}/dependencies/moose-2026.06.16-linux-64.yml" \
    "netcdf4=1.7.3"

workflow="${source_root}/.github/workflows/source-audit.yml"
require_text "${workflow}" "runs-on: [self-hosted, linux, x64, fuelsim]"
require_text "${workflow}" "-DFUELSIM_WARNINGS_AS_ERRORS=ON"
require_text "${workflow}" "-DFUELSIM_ENABLE_SANITIZERS="
require_text "${workflow}" \
    "ctest --test-dir build-ci-release -j4 --output-on-failure"
require_text "${workflow}" \
    "ctest --test-dir build-ci-sanitizer -j4 --output-on-failure"
require_text "${workflow}" "MPIR_CVAR_CH4_NETMOD=ofi"
require_text "${workflow}" "FI_PROVIDER=tcp"

for path in \
    "${source_root}/AGENTS.md" \
    "${source_root}/README.md" \
    "${source_root}/docs/reproducible-build.md" \
    "${source_root}/docs/verification.md"; do
    if grep -Fq -- "/tmp/adlite-fuelsim-install" "${path}"; then
        echo "[FAIL] ${path} still uses the ephemeral ADlite prefix" >&2
        exit 1
    fi
done

# The dependency scripts must reject unpinned or dirty source trees with
# exit code 2 before any configure, build, or install action runs. The fake
# toolchain compilers only satisfy the existence checks; they fail loudly if
# a script ever reaches a compiler invocation, and every case below also
# asserts that no build or install directory was created.
work_root=$(mktemp -d)
trap 'rm -rf "${work_root}"' EXIT

fake_toolchain="${work_root}/toolchain"
mkdir -p "${fake_toolchain}/bin"
for tool in cc c++; do
    cat > "${fake_toolchain}/bin/${tool}" <<'EOF'
#!/usr/bin/env bash
echo "fake toolchain compiler must never run" >&2
exit 1
EOF
    chmod +x "${fake_toolchain}/bin/${tool}"
done

make_git_tree() {
    local path=$1
    mkdir -p "${path}"
    git -c init.defaultBranch=main -C "${path}" init --quiet
    echo "cmake_minimum_required(VERSION 3.16)" > "${path}/CMakeLists.txt"
    git -C "${path}" add CMakeLists.txt
    git -C "${path}" -c user.name=fuelsim-test \
        -c user.email=fuelsim-test@localhost commit --quiet -m "fake tree"
}

expect_rejection() {
    local description=$1
    local script=$2
    local source_dir=$3
    local build_dir="${work_root}/build"
    local install_prefix="${work_root}/install"
    rm -rf "${build_dir}" "${install_prefix}"
    local status=0
    bash "${script}" "${source_dir}" "${install_prefix}" "${build_dir}" \
        "${fake_toolchain}" > /dev/null 2>&1 || status=$?
    if [[ ${status} -ne 2 ]]; then
        echo "[FAIL] ${description}: expected exit code 2, got ${status}" >&2
        exit 1
    fi
    if [[ -e "${build_dir}" || -e "${install_prefix}" ]]; then
        echo "[FAIL] ${description}: rejection happened only after build or install started" >&2
        exit 1
    fi
    echo "[PASS] ${description}"
}

adlite_tree="${work_root}/adlite"
make_git_tree "${adlite_tree}"
expect_rejection "install_adlite.sh rejects an unpinned ADlite commit" \
    "${source_root}/scripts/install_adlite.sh" "${adlite_tree}"

plain_tree="${work_root}/adlite-not-git"
mkdir -p "${plain_tree}"
echo "cmake_minimum_required(VERSION 3.16)" > "${plain_tree}/CMakeLists.txt"
expect_rejection "install_adlite.sh rejects a non-Git ADlite source" \
    "${source_root}/scripts/install_adlite.sh" "${plain_tree}"

seacas_tree="${work_root}/seacas"
make_git_tree "${seacas_tree}"
git -C "${seacas_tree}" tag v1999-01-01
expect_rejection "build_exodus.sh rejects an unpinned SEACAS tag" \
    "${source_root}/scripts/build_exodus.sh" "${seacas_tree}"

dirty_seacas_tree="${work_root}/seacas-dirty"
make_git_tree "${dirty_seacas_tree}"
git -C "${dirty_seacas_tree}" tag v2024-06-27
echo "uncommitted change" >> "${dirty_seacas_tree}/CMakeLists.txt"
expect_rejection "build_exodus.sh rejects a dirty SEACAS worktree" \
    "${source_root}/scripts/build_exodus.sh" "${dirty_seacas_tree}"

echo "[PASS] pinned dependency and numerical CI source audit"
