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
    "a3778d2b2fd87cd70da59d1cc32ed3b856020a3e"
require_text "${source_root}/scripts/build_exodus.sh" "v2024-06-27"
require_text \
    "${source_root}/dependencies/moose-2026.06.16-linux-64.yml" \
    "moose-dev=2026.06.16=mpich"

workflow="${source_root}/.github/workflows/source-audit.yml"
require_text "${workflow}" "runs-on: [self-hosted, linux, x64, fuelsim]"
require_text "${workflow}" "-DFUELSIM_WARNINGS_AS_ERRORS=ON"
require_text "${workflow}" "-DFUELSIM_ENABLE_SANITIZERS="
require_text "${workflow}" \
    "ctest --test-dir build-ci-release --output-on-failure"
require_text "${workflow}" \
    "ctest --test-dir build-ci-sanitizer --output-on-failure"
require_text "${workflow}" "MPIR_CVAR_CH4_NETMOD=ofi"
require_text "${workflow}" "FI_PROVIDER=tcp"

for path in \
    "${source_root}/README.md" \
    "${source_root}/docs/reproducible-build.md" \
    "${source_root}/docs/verification.md"; do
    if grep -Fq -- "/tmp/adlite-fuelsim-install" "${path}"; then
        echo "[FAIL] ${path} still uses the ephemeral ADlite prefix" >&2
        exit 1
    fi
done

echo "[PASS] pinned dependency and numerical CI source audit"
