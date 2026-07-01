#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:-/tmp/fylat_build_check}"

if ! python -c "import pybind11" >/dev/null 2>&1; then
  echo "pybind11 is not importable. Activate the cloudmask conda environment first." >&2
  exit 2
fi

pybind11_dir="$(python -c "import pybind11; print(pybind11.get_cmake_dir())")"

echo "Repository : ${repo_root}"
echo "Build dir  : ${build_dir}"
echo "Python     : $(command -v python)"
echo "pybind11   : ${pybind11_dir}"

cmake -S "${repo_root}" -B "${build_dir}" -Dpybind11_DIR="${pybind11_dir}"
cmake --build "${build_dir}" -j"${FYLAT_BUILD_JOBS:-2}"
ctest --test-dir "${build_dir}" --output-on-failure
