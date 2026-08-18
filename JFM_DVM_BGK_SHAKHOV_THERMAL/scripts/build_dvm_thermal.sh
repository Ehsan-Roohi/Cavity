#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

if ! command -v gfortran >/dev/null 2>&1; then
  if command -v module >/dev/null 2>&1; then
    module load gcc
  fi
fi
if ! command -v gfortran >/dev/null 2>&1; then
  echo "ERROR: gfortran is unavailable. On Unity, load a GCC module first." >&2
  exit 2
fi

mkdir -p bin
gfortran -O3 -march=native -funroll-loops -fopenmp \
  -ffree-line-length-none -fimplicit-none -Wall -Wextra \
  src/dvm2d_thermal_cavity.f90 \
  -o bin/dvm2d_thermal_cavity

echo "Built: $ROOT/bin/dvm2d_thermal_cavity"
gfortran --version | head -1

