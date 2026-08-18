#!/bin/bash
#SBATCH --job-name=jfm-s-full30
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-jfm-shakhov-full-kn30-%j.out
#SBATCH --error=slurm-jfm-shakhov-full-kn30-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_STACKSIZE=256M

if ! command -v gfortran >/dev/null 2>&1; then
  module load gcc
fi

mkdir -p bin
gfortran -O3 -march=native -funroll-loops -fopenmp \
  -ffree-line-length-none -fimplicit-none -Wall -Wextra \
  src/dvm2d_thermal_cavity.f90 \
  -o "bin/dvm2d_thermal_cavity_${SLURM_JOB_ID}"
export DVM_BINARY="$PWD/bin/dvm2d_thermal_cavity_${SLURM_JOB_ID}"

# Exact JFM case discussed in the revision work:
# stationary full square cavity; hot horizontal walls T=1; cold vertical
# walls T=0.2; paper Kn=30; Shakhov Pr=2/3.
mamba run -n dsmc-gpu bash scripts/run_dvm_thermal_case.sh \
  shakhov 0.2 30 production full

RESULT_REL="results_dvm_jfm/shakhov_full_Kn30_RT0p2_production_N102_V61"
METADATA="$RESULT_REL/ThermalCavity_DVM_shakhov_full_Kn30_RT0p2_metadata.json"
ARCHIVE="JFM_DVM_SHAKHOV_FULL_Kn30_RT0p2_job${SLURM_JOB_ID}.zip"

if command -v zip >/dev/null 2>&1; then
  zip -qr "$ARCHIVE" "$RESULT_REL"
else
  mamba run -n dsmc-gpu python -m zipfile -c "$ARCHIVE" "$RESULT_REL"
fi

echo "ARCHIVE=$PWD/$ARCHIVE"
echo "METADATA=$PWD/$METADATA"

if ! grep -Eq '"converged"[[:space:]]*:[[:space:]]*true' "$METADATA"; then
  echo "WARNING: run ended without satisfying the steady residual gate; keep the ZIP for diagnosis, not as a JFM production result." >&2
  exit 3
fi

echo "JFM full-cavity Shakhov production gate passed."
