# Transmission-delay quantitative verification

This isolated experiment recovers one homogeneous target permittivity from
transmission delays and an exactly known target mask. Target material values
are used only to score the reconstruction.

## Recommended reconstruction

The validated workflow matches object/reference pulse envelopes, measures
centerline chords through the known mask, and estimates refractive-index
contrast from the maximum target chord:

```matlab
verificationSeed = 1;
forceRecompute = false;
run('inverse_problem/algorithms/td_verification/tds_2d_improved.m')
```

It reduced mean absolute percentage error from 36.64% to 2.79% across eight
deterministic cases. See `documentation.md` for every attempted estimator,
the development/holdout results, mathematical reasoning, and limitations.

## Other entry points

```matlab
% Reproduce the original 15-pixel arithmetic-mixture baseline.
run('inverse_problem/algorithms/td_verification/build_cfg.m')
run('inverse_problem/algorithms/td_verification/tds_2d.m')

% Load or acquire seeds 1-8 and compare baseline versus improved errors.
run('inverse_problem/algorithms/td_verification/run_experiment_suite.m')
```

Per-seed measurements are cached under `results/` to make estimator studies
repeatable without rerunning FDTD.
