#!/bin/bash
#SBATCH --job-name=jfm-dvm-thermal
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-3
#SBATCH --output=slurm-dvm-%A_%a.out
#SBATCH --error=slurm-dvm-%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Each task compiles privately to avoid an array-start race on the executable.
if ! command -v gfortran >/dev/null 2>&1; then
  module load gcc
fi
mkdir -p bin
gfortran -O3 -march=native -funroll-loops -fopenmp \
  -ffree-line-length-none -fimplicit-none -Wall -Wextra \
  src/dvm2d_thermal_cavity.f90 \
  -o "bin/dvm2d_thermal_cavity_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export DVM_BINARY="$PWD/bin/dvm2d_thermal_cavity_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

MODELS=(bgk shakhov bgk shakhov)
RTS=(0.2 0.2 0.5 0.5)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
RT="${RTS[$SLURM_ARRAY_TASK_ID]}"

mamba run -n dsmc-gpu bash scripts/run_dvm_thermal_case.sh \
  "$MODEL" "$RT" production
