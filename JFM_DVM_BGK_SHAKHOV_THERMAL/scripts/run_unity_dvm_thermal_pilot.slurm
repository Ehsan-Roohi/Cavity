#!/bin/bash
#SBATCH --job-name=jfm-dvm-pilot-v3
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-1%1
#SBATCH --output=slurm-dvm-pilot-v3-%A_%a.out
#SBATCH --error=slurm-dvm-pilot-v3-%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
if ! command -v gfortran >/dev/null 2>&1; then
  module load gcc
fi

# Build a versioned executable so an older binary cannot be reused.
mkdir -p bin
gfortran -O3 -march=native -funroll-loops -fopenmp \
  -ffree-line-length-none -fimplicit-none -Wall -Wextra \
  src/dvm2d_thermal_cavity.f90 \
  -o "bin/dvm2d_thermal_cavity_v3_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export DVM_BINARY="$PWD/bin/dvm2d_thermal_cavity_v3_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

MODELS=(bgk shakhov)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

mamba run -n dsmc-gpu bash scripts/run_dvm_thermal_case.sh \
  "$MODEL" 0.2 pilot
