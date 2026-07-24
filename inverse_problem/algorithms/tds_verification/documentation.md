# Transmission-delay permittivity experiments

## Objective and fixed conditions

The experiment assumes that the target mask is known exactly and that the
target is homogeneous. The target permittivity is hidden from every estimator
and is used only after reconstruction to calculate error. All tests retain the
16-antenna circular array, the 40-degree off-boresight limit, the 31-by-31
rectangular target, and directed TX-to-RX measurements.

Seeds 1-5 were used to investigate methods. After choosing the final rule, its
parameters were frozen and seeds 6-8 were run as unseen holdout cases.

## Mathematical correction

The original experiment treated the delay-derived quantity as an arithmetic
average of relative permittivity. Travel time instead integrates refractive
index (or slowness):

\[
c_0\Delta t_i = \int_{L_i}\left(n(\mathbf x)-n_b\right)\,ds,
\qquad n=\sqrt{\epsilon_r}.
\]

For a homogeneous target with known intersection length \(L_{t,i}\), the
straight-ray estimate is therefore

\[
\widehat n_t=n_b+\frac{c_0\Delta t_i}{L_{t,i}},
\qquad \widehat\epsilon_{r,t}=\widehat n_t^2.
\]

This agrees with the travel-time model in
[Tepe, Schuster, and Littau (2017)](https://doi.org/10.1080/17415977.2016.1267168),
which reconstructs refractive index from path differences and known
interfaces. Radar travel-time tomography likewise inverts slowness rather than
permittivity directly
([Day-Lewis, Singha, and Binley, 2005](https://doi.org/10.1029/2004JB003569)).

## Experiments attempted

### 1. Original 15-pixel arithmetic-permittivity mixture

The original formula was retained as the baseline. Seed 1 recovered 2.2765
instead of 3.8887, an error of 41.46%. Across all eight cases its MAPE was
36.64%, with a worst error of 58.06%.

Conclusion: pixel-count mixing in the permittivity domain is not compatible
with the measured optical delay.

### 2. Refractive-index mixing with thick strips

The same 1-40 pixel strip masks were tested after changing the mixture to the
refractive-index domain. On seed 1, width 14 happened to give 0.45% error and
width 15 gave 1.05% error. This did not generalize: the best common width on
seeds 1-5 was 19 pixels with 11.73% MAPE, and one development case still had
26.13% error.

The seed-1 result was caused partly by cancellation: 32 of 54 accepted
15-pixel rays had zero prominent-peak delay, while a small number of grazing
rays produced very large estimates. Changing the width by only a few pixels
changed the result sharply; for example seed-1 width 6 exceeded 200% error.

Conclusion: correcting the material domain is necessary, but an arbitrary
finite-width strip and unweighted per-ray mean are not stable.

### 3. Straight centerline chords with the first prominent peak

Centerline intersection lengths were measured directly from the known mask.
The prominent-peak delay was bimodal. Several rays with 25-41 pixel target
chords reported only -2 to 6 delay steps, while other comparable chords
reported 50-54 steps. Clamping the early values to zero made the centerline
estimate unusable.

Interpretation: for this high contrast, the first detectable energy can be a
diffracted or target-avoiding arrival rather than the transmitted pulse whose
delay is required. Day-Lewis et al. warn that straight rays become
inappropriate above roughly 10% velocity contrast and distinguish
infinitesimal rays from finite-frequency sensitivity volumes. The simulated
contrast is far beyond that regime.

### 4. Matched-envelope delay

A normalized envelope-matching estimator now aligns each reference pulse with
the object trace over nonnegative lags and uses parabolic interpolation for a
sub-sample result. On seed 1 it changed the inconsistent centerline delays to
a coherent range of approximately 49-59 steps, with a median match score of
0.9978.

Using every matched-delay centerline still overestimated permittivity because
short and grazing straight chords do not represent the actual high-contrast
wave path. Thick-strip matched-delay estimates were also poor: their best
common width on seeds 1-5 had 71.37% MAPE.

Conclusion: improved temporal alignment fixes peak association, but geometric
ray selection remains essential.

### 5. Maximum-chord matched-delay estimator

The final method selects only directed rays whose straight centerline target
intersection is at least 99.9% of the maximum available chord. It solves one
global least-squares equation in refractive-index contrast:

\[
\widehat{\Delta n}=
\frac{\sum_i L_{t,i}(c_0\Delta t_i)}{\sum_i L_{t,i}^2},
\qquad
\widehat\epsilon_{r,t}=(n_b+\widehat{\Delta n})^2.
\]

This rule is geometry-only: it does not choose rays or thresholds using the
known target permittivity. Long chords reduce division by a small uncertain
intersection and favor paths that transmit through the target interior rather
than graze an interface.

## Final results

| Seed | True epsr | Original epsr | Original error | Maximum-chord epsr | New error | Selected directed rays |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3.8887 | 2.2765 | 41.46% | 3.7866 | 2.62% | 2 |
| 2 | 4.0236 | 2.5722 | 36.07% | 4.1379 | 2.84% | 4 |
| 3 | 4.2879 | 2.2868 | 46.67% | 4.6287 | 7.95% | 2 |
| 4 | 5.9419 | 3.4763 | 41.49% | 5.8982 | 0.74% | 2 |
| 5 | 3.5388 | 2.0669 | 41.59% | 3.7316 | 5.45% | 2 |
| 6 (holdout) | 2.6459 | 2.6670 | 0.80% | 2.7105 | 2.44% | 4 |
| 7 (holdout) | 5.0321 | 2.1106 | 58.06% | 5.0420 | 0.20% | 4 |
| 8 (holdout) | 4.3580 | 3.1826 | 26.97% | 4.3605 | 0.06% | 2 |

- Development seeds 1-5: 3.92% MAPE, 7.95% worst error.
- Frozen holdout seeds 6-8: 0.90% MAPE, 2.44% worst error.
- All eight cases: 2.79% MAPE, 7.95% worst error.
- Original all-case baseline: 36.64% MAPE, 58.06% worst error.

## Current conclusion and limitations

For the present single homogeneous rectangle, the quantitative information is
recoverable accurately when delay is extracted from the transmitted pulse,
the inversion is performed in refractive index, and only maximum-intersection
paths are used. Reducing strip background pixels alone was not sufficient.

The method currently relies on only one or two independent reciprocal chords,
so it may become fragile under measurement noise, mask error, irregular or
concave shapes, and targets for which no antenna chord crosses a substantial
interior length. It also remains a straight-ray approximation despite the
large contrast. The next defensible extension would trace refracted paths
through the known boundary using Snell's law, or fit the single homogeneous
permittivity directly with a full-wave forward model. Those extensions were
not needed to reach low error in the eight present cases.

## Reproducing the study

Run one improved case:

```matlab
verificationSeed = 1;
forceRecompute = false;
run('inverse_problem/algorithms/td_verification/tds_2d_improved.m')
```

Run or reload all eight cases:

```matlab
run('inverse_problem/algorithms/td_verification/run_experiment_suite.m')
```

Measurements are cached in `results/measurement_seed_XXXX.mat`. Set
`forceRecompute = true` or call `runMeasurementCase(seed, true)` to rerun FDTD.
