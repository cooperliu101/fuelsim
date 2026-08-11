#!/usr/bin/env bash

set -euo pipefail

repository_root="$(git rev-parse --show-toplevel)"
cd "${repository_root}"

formatter="${1:-${CLANG_FORMAT:-clang-format}}"
if [[ "${formatter}" == */* ]]; then
    if [[ ! -x "${formatter}" ]]; then
        echo "clang-format is not executable: ${formatter}" >&2
        exit 1
    fi
    formatter_path="$(cd "$(dirname "${formatter}")" && pwd)/$(basename "${formatter}")"
else
    formatter_path="$(command -v "${formatter}" || true)"
    if [[ -z "${formatter_path}" ]]; then
        echo "clang-format was not found: ${formatter}" >&2
        echo "Usage: scripts/install_git_hooks.sh /path/to/clang-format" >&2
        exit 1
    fi
fi

"${formatter_path}" --version
git config --local core.hooksPath .githooks
git config --local hooks.clangFormat "${formatter_path}"

echo "Installed fuelsim Git hooks from .githooks."
echo "The pre-commit hook will use ${formatter_path}."
