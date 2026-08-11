# Deterministic DVM for the JFM stationary thermal cavity

This package adapts the recovered `dvm2d_cavity_shakhov.f90` lid-driven
cavity code to the stationary thermal cavity in JFM-2026-1491.

## Physics implemented

- Square cavity with **no external force** and no moving wall.
- Horizontal walls: `T_hot=1`; vertical walls: `T_cold=RT`.
- Fully diffuse, impermeable physical walls.
- Quarter domain `x,y in [0,1/2]`, with exact specular symmetry conditions
  at `x=0` and `y=0`; full-domain fields are reconstructed by parity.
- Selectable model-equation target: `bgk` or `shakhov` (`Pr=2/3`).
- The Shakhov correction uses the same recorded `[0,2]` weight limiter as
  the GPU particle implementation.  Limiting is applied before the local
  conservative projection, so discrete mass, momentum and energy remain
  conserved.
- Exact paper Knudsen-number scaling:

  ```text
  Kn = lambda/L = 1/(sqrt(2)*n0*pi*d^2*L)
  nu = rho*sqrt(pi*T)/2
  tau = Kn/nu = 2*Kn/(rho*sqrt(pi*T))
  ```

  There is no `Kn=30 -> 96` conversion and no viscosity-based hidden
  rescaling.

The solver uses the two reduced distributions `h` and `b`, retaining the
three-dimensional monatomic energy.  A cell-local positive exponential
projection uses four Lagrange multipliers to remove discrete
velocity-quadrature defects in mass, both momentum components, and total
energy before each BGK/Shakhov relaxation update.  Unlike the earlier linear
projection, this correction cannot create a negative target distribution.

This is a deterministic **model-equation DVM**, not a deterministic solver for
the full hard-sphere Boltzmann collision integral.  The hard-sphere reference
continues to be DSMC.

## Prepared production cases

The production Slurm array contains the 12 BGK/Shakhov cases needed for the
current revision campaign:

| Task | Model | Paper Kn | RT |
|---:|---|---:|---:|
| 0--5 | BGK/Shakhov | 20, 10, 5 | 0.2 |
| 6--11 | BGK/Shakhov | 20, 10, 5 | 0.5 |

The single-case runner is general in paper Kn and is called as
`run_dvm_thermal_case.sh MODEL RT KN LEVEL`.  For example:

```bash
OMP_NUM_THREADS=16 bash scripts/run_dvm_thermal_case.sh shakhov 0.2 10 production
```

Production resolution is `51 x 51` cells in the quarter domain (equivalent to
`102 x 102` reconstructed full-domain cells), a `61 x 61` velocity grid on
`[-5,5]^2`, and at most 120,000 time steps.  The solver stops after 20,000
steps if the logged residual falls below `1e-8`.

## Unity submission

DVM is CPU/OpenMP code.  Reserving a GPU would not accelerate this solver.

After checking out this repository on Unity, the one-line pilot command is:

```bash
cd /project/pi_roohie_umass_edu/github_sync/Cavity/JFM_DVM_BGK_SHAKHOV_THERMAL && sbatch scripts/run_unity_dvm_thermal_pilot.slurm
```

After both pilot cases pass the acceptance checks, the one-line four-case
production command is:

```bash
cd /project/pi_roohie_umass_edu/github_sync/Cavity/JFM_DVM_BGK_SHAKHOV_THERMAL && sbatch scripts/run_unity_dvm_thermal_cpu.slurm
```

The July 21, 2026 `25 x 25`, `41 x 41`, 5,000-step pilot compiled and ran on
Unity with zero projection failures, zero negative targets, wall mass flux at
roundoff, and relative mass drift below `5e-15`.  It did **not** satisfy the
steady residual criterion (`final_residual` was approximately `8.9e-2`), so
that pilot demonstrates executable integrity only and must not be used as a
converged JFM result.

```bash
cd /project/pi_roohie_umass_edu/JFM_revision_2026
unzip -o JFM_DVM_BGK_SHAKHOV_THERMAL.zip
cd DVM2D_THERMAL_JFM
bash -n scripts/run_unity_dvm_thermal_cpu.slurm
sbatch scripts/run_unity_dvm_thermal_cpu.slurm
```

For a compilation and 20-step smoke check before submitting production:

```bash
cd /project/pi_roohie_umass_edu/JFM_revision_2026/DVM2D_THERMAL_JFM
bash scripts/build_dvm_thermal.sh
OMP_NUM_THREADS=8 bash scripts/run_dvm_thermal_case.sh bgk 0.2 20 smoke
OMP_NUM_THREADS=8 bash scripts/run_dvm_thermal_case.sh shakhov 0.2 20 smoke
```

The same check can be submitted as a two-task CPU array:

```bash
sbatch scripts/run_unity_dvm_thermal_smoke.slurm
```

After smoke succeeds, run the supplied BGK/Shakhov pilot for the stronger
`RT=0.2` case.  It uses `25 x 25` quarter cells, a `41 x 41` velocity grid and
5,000 iterations:

```bash
sbatch scripts/run_unity_dvm_thermal_pilot.slurm
```

Do not submit production until both pilot metadata files pass the acceptance
checks below.

For one production case only (for example, task 0 = BGK, Kn=20, RT=0.2):

```bash
sbatch --array=0 scripts/run_unity_dvm_thermal_cpu.slurm
```

## Outputs

Current cases are written below `results_dvm_jfm/`, leaving the earlier pilot
results untouched for comparison. Each case contains:

- `*_full.csv`: reconstructed full-domain quantitative fields;
- `*_full.dat`: Tecplot point data;
- `*_quarter.dat`: raw quarter-domain block data;
- `*_quarter_moments.dat`: higher moments;
- `*_metadata.json`: definition, grid, convergence, conservation and wall
  diagnostics;
- `*_postcheck.json`: exact parity checks from the reconstructed fields;
- `*_fields.png`: unfiltered temperature, density and speed plots;
- `*.hst` and `stdout.log`: convergence histories.

The output heat flux is the paper convention
`q = integral(|c|^2 c f dc)`.  Internally the common half-energy convention is
used and multiplied by two when written.

## Acceptance checks

Do not use a completed case in the referee response until its metadata shows:

- `converged: true`;
- `mass_drift_relative` small (target below `1e-8`);
- `max_physical_wall_mass_flux` near roundoff;
- `max_discrete_target_conservation_error` near roundoff;
- `clipped_updated_values: 0` (required for both models); the Shakhov
  metadata also reports the limiter interval, raw weight extrema and number
  of limited target values;
- `positive_projection_failures: 0`; the positive exponential conservative
  projection must match the four collision invariants within its recorded
  tolerance;
- a stable velocity field under the optional `fine` resolution.

If a production case is smooth but has not met the residual tolerance at
120,000 steps, retain the output for diagnosis but do not label it converged.
The `fine` level (`81 x 81`, `81 x 81` velocity grid) is available through the
manual run script after the standard cases are assessed.
