# Cavity and DSMC learning notebooks

Exploratory notebooks and scripts for lid-driven cavity flow, DSMC relaxation,
shock-wave surrogates, PINNs, and supervised neural-network comparisons.

## Contents

- `Supervised_Learning_Cavity_Re=100.ipynb` — supervised cavity-flow study.
- `DSMC_Relaxation_Machine_Learning_*.ipynb` and `Relaxation_DSMC*.ipynb`
  — relaxation and PINN experiments.
- `Supervised_Shock_wave_DSMC*.ipynb` — shock-wave surrogate experiments.
- `Shock_Wave_Rotational_NN.py` — TensorFlow model for velocity and
  translational/rotational-temperature profiles; it expects local `M*.txt`
  data files that are not distributed here.
- `Cavity Re=100 CFD.ipynb` — despite its historical extension, this file is a
  plain Python script, not a valid Jupyter notebook.

Most valid notebooks include an “Open in Colab” badge.  Otherwise clone the
repository and install the imports declared by the selected experiment
(typically NumPy, Matplotlib, scikit-learn, and TensorFlow).

## Scientific status

This is a research/teaching archive, not a packaged or benchmark-qualified CFD
solver.  Filenames preserve the experiment history.  Before citing a numerical
result, record the exact commit, input data, split/seed, mesh or DSMC settings,
and an independent reference comparison.

The deterministic JFM thermal-cavity draft is tracked separately in
[`RESEARCH_STATUS.md`](RESEARCH_STATUS.md); it is not part of the default
notebook baseline.

## Data, citation, and license

Large input/output data are not versioned consistently across the historical
experiments.  No repository-wide citation metadata or license is declared yet;
choose them only after confirming authorship and the intended release scope.
