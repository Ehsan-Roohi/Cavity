#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MODEL="${1:-shakhov}"
RT="${2:-0.2}"
LEVEL="${3:-production}"

if [[ "$MODEL" != "bgk" && "$MODEL" != "shakhov" ]]; then
  echo "MODEL must be bgk or shakhov" >&2
  exit 2
fi
if [[ "$RT" != "0.2" && "$RT" != "0.5" ]]; then
  echo "RT must be 0.2 or 0.5 for the prepared JFM cases" >&2
  exit 2
fi

case "$LEVEL" in
  smoke)
    NX=13; NY=13; NV=21; STEPS=20; MIN_STEPS=20; LOG_EVERY=5
    ;;
  pilot)
    NX=25; NY=25; NV=41; STEPS=5000; MIN_STEPS=3000; LOG_EVERY=100
    ;;
  production)
    NX=51; NY=51; NV=61; STEPS=120000; MIN_STEPS=20000; LOG_EVERY=200
    ;;
  fine)
    NX=81; NY=81; NV=81; STEPS=180000; MIN_STEPS=30000; LOG_EVERY=200
    ;;
  *)
    echo "LEVEL must be smoke, pilot, production, or fine" >&2
    exit 2
    ;;
esac

BIN="${DVM_BINARY:-$ROOT/bin/dvm2d_thermal_cavity}"
if [[ ! -x "$BIN" ]]; then
  bash "$ROOT/scripts/build_dvm_thermal.sh"
  BIN="$ROOT/bin/dvm2d_thermal_cavity"
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_STACKSIZE="${OMP_STACKSIZE:-256M}"

RTTAG=${RT/./p}
OUT="$ROOT/results_dvm_kn30_v3/${MODEL}_Kn30_RT${RTTAG}_${LEVEL}_N${NX}_V${NV}"
mkdir -p "$OUT"
cd "$OUT"

cat > dvm_cavity.in <<EOF
&params
 nx = ${NX},
 ny = ${NY},
 nv = ${NV},
 vmin = -5.0,
 vmax = 5.0,
 kn = 30.0,
 t_hot = 1.0,
 t_cold = ${RT},
 pr = 0.6666666666666666667,
 cfl = 0.45,
 tol = 1.0e-8,
 floor_val = 0.0,
 shakhov_weight_min = 0.0,
 shakhov_weight_max = 2.0,
 projection_max_iter = 12,
 projection_tol = 1.0e-12,
 steps = ${STEPS},
 min_steps = ${MIN_STEPS},
 log_every = ${LOG_EVERY},
 save_every = 0,
 model = '${MODEL}',
 out_prefix = 'ThermalCavity_DVM_${MODEL}_Kn30_RT${RTTAG}',
/
EOF

"$BIN" 2>&1 | tee stdout.log

PREFIX="ThermalCavity_DVM_${MODEL}_Kn30_RT${RTTAG}"
python "$ROOT/tools/plot_dvm_thermal.py" \
  --csv "${PREFIX}_full.csv" \
  --metadata "${PREFIX}_metadata.json" \
  --output-prefix "$PREFIX"

echo "DONE: $OUT"
find . -maxdepth 1 -type f -printf '%f\n' | sort
