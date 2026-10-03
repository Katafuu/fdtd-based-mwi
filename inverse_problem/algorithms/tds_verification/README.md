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
run('inverse_problem/algorithms/tds_verification/tds_2d_improved.m')
```

It reduced mean absolute percentage error from 36.64% to 2.79% across eight
deterministic cases. See `documentation.md` for every attempted estimator,
the development/holdout results, mathematical reasoning, and limitations.

## Other entry points

```matlab
% Reproduce the original 15-pixel arithmetic-mixture baseline.
addpath('inverse_problem/algorithms/tds_verification')
setup = struct('targetOpts', struct(), 'randomSeed', 1);
cfg = build_cfg(setup);
run('inverse_problem/algorithms/tds_verification/tds_2d.m')

% Load or acquire seeds 1-8 and compare baseline versus improved errors.
run('inverse_problem/algorithms/tds_verification/run_experiment_suite.m')
```

Per-seed measurements are cached under `results/` to make estimator studies
repeatable without rerunning FDTD.

## Configuring targets and material variation

`build_cfg()` creates a fixed off-center circle. Pass `setup.targetSpecs` for
exact targets, `setup.targetOpts` for random targets, or both; use
`setup.randomSeed` to repeat random generation. The eight-case measurement
runner explicitly requests its original seeded random rectangle, so its
cached study remains comparable.

The shared `fdtdmat.applyEpsrInhomogeneity(cfg, mask, opts)` function is
available after `build_cfg` adds `buildLib` to the path. It applies bounded,
correlated variation to selected epsilon-r cells and returns the varied
configuration plus diagnostic `info`. The current transmission-delay
estimator still verifies one homogeneous target at a time.

## Shape-estimation robustness

Run `shape_estimate_robust_test.m` from this folder to test target/background
inhomogeneity and receiver noise. It explicitly uses the 12-antenna
`hybrid_tr_mwi/build_cfg.m` configuration and hybrid TR helpers; the local
verification builder uses 16 antennas. Results are saved under
`figs/shape_estimate/run_XXXX/`.
