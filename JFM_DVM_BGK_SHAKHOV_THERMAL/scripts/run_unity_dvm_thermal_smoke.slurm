#!/bin/bash
#SBATCH --job-name=jfm-dvm-smoke
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=0-1%1
#SBATCH --output=slurm-dvm-smoke-%A_%a.out
#SBATCH --error=slurm-dvm-smoke-%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
if ! command -v gfortran >/dev/null 2>&1; then
  module load gcc
fi
bash scripts/build_dvm_thermal.sh

MODELS=(bgk shakhov)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

mamba run -n dsmc-gpu bash scripts/run_dvm_thermal_case.sh \
  "$MODEL" 0.2 20 smoke
