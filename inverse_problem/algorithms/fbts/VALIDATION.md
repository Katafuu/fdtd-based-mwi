# Acceleration validation, 2026-09-28

Current evidence is saved under `results/acceleration_validation/` (ignored by
Git so numerical artifacts do not enter source control):

- `tests.mat` and `final_tests.log`: FBTS regressions, ROI equality, geometry,
  threading, checkpoint state capture and batch recovery.
- `forward_tests.mat` and `forward_tests_fixed.log`: original forward-solver
  tests using a native Linux standalone executable. The pre-existing Windows
  path in the example builder was corrected; the test harness accepts
  `FDTD_STANDALONE_EXECUTABLE` without overwriting the supplied Windows binary.
- `serial_parity_bg*_active*.mat` and `serial_parity.log`: serial-build versus
  OpenMP-build reference/actual outputs for six three-iteration controls,
  exactly equal in RX, gradients, sensitivity, image and optimizer histories.
- `controlled.log`: exact one/two-thread high-background controls across all
  four shapes, including full, alternating and singleton sets.
- `window/`: persistent full-grid raw/scattered recording-window scans.
  The 12 ns window failed (worst tail fraction approximately 5.35e-4 against
  1e-4). Longer windows require a completed report before study preparation.

The current machine has Intel integrated graphics and no installed Parallel
Computing Toolbox. CUDA compilation/device testing and process-pool testing
have **not** passed here. Small CPU controls are not substitutes for full-size
physical calibration or intended-PC throughput measurements. No validated
GTX-PC runtime estimate or completed 348-case study is claimed.

The following is the historical pre-acceleration batch validation. Its timings
and grid/medium assumptions are not those of the new high-background study.

---

# FBTS batch validation

Validated locally on 2026-09-25 using MATLAB and the repository's compiled
FDTD MEX. The C solver was not changed. The regression suite command and public
usage are in [README.md](README.md).

Final regression coverage: **18 passed, 0 failed, 1 skipped**. The full-suite
run passed all six standalone/coarse-sensitivity tests; after fixing the
missing-status recovery edge case, all twelve batch tests passed in the affected
suite rerun. The skipped test requires Parallel Computing Toolbox.

## Numerical and data checks

- All 15 nonempty subsets of a four-antenna scene: sliced full-array measurements
  equal independently acquired measurements exactly.
- Cached and independently acquired reconstructions: exact equality of
  non-timing results over three iterations, sensitivity factors 1 and 4,
  full arrays, reduced arrays, and a single active antenna.
- Dihedral catalogs: 29 classes for eight antennas and 223 for twelve,
  excluding the empty active array. Distinct orbit sizes sum to all nonempty
  subsets. Cyclic shifts and reflection merge; arbitrary gap permutations do
  not. Seeded representatives do not alter the caller's RNG state.
- Interrupted four-iteration inversion: resuming after iteration two reproduces
  the uninterrupted conjugate-gradient trajectory exactly. Modified committed
  history is rejected. Failure diagnostics retain the interrupted phase.
- Shared scene data: exact truth/mask/image alignment, complete iteration
  histories, original/local antenna mapping, headless execution, no saved image
  quality metrics, and unchanged caller configuration.
- Moving a whole dataset preserves loading. Analysis loading works without the
  FDTD MEX or buildLib. Saved scenes share identical system files while keeping
  separate acquisitions. Resume verifies hashes and skips completed work.
- Single-antenna/single-iteration MAT histories, changed-file rejection,
  stationary-inversion failure isolation, blocked output directories, and
  coordinator lock exclusion are covered.
- Legacy standalone demo and coarse sensitivity regression tests are retained.

Parallel execution has an optional serial-versus-process-worker comparison.
The local installation lacks Parallel Computing Toolbox, so that test is
skipped. Multiworker numerical behavior and memory usage therefore remain
unvalidated on this machine; serial execution is the default.

## Bounded full-resolution pilot

The pilot used `buildFbtsConfig(struct('seed',42,'numAntennas',8))`, two FBTS
iterations, sensitivity downsample factor 4, and batch options
`struct('disabledCounts',[0 7],'seed',123)`. Both selected cases completed:
the full eight-antenna array and the singleton whose original index was six.

| Measurement | Observed value |
| --- | --- |
| Fine grid / time samples | 400 × 400 / 600 |
| Initial acquisition solves | 8, shared by both cases |
| Iterative fine solves | 36: forward plus adjoint |
| Iterative coarse solves | 36: baseline plus perturbation |
| Batch elapsed time | 121.79 s |
| MATLAB process elapsed time, including startup | 132.66 s |
| Peak process resident memory | 2,553,864 KiB, approximately 2.44 GiB |
| Saved dataset size, including provenance | 3,156,926 bytes, approximately 3.16 MB |

The pilot checked successful status, finite reconstructions, equality between
the final compact history column and the final image, and exact saved truth.
Peak memory was measured using GNU time. These are measurements of this local
two-case run, not a budget or runtime prediction for an entire study. Storage
grows with cases, iterations, DOI pixels and retained RX channels. Workers need
their own field histories and temporary arrays.

The pilot preceded a final resume-only refinement that makes saved catalogs and
manifest statuses authoritative when recovering, including recreating a missing
status file after verifying its result. The full regression suite and then the
affected batch suite were rerun during this refinement; the numerical inversion
code was unchanged.

The synthetic acquisition and inversion use the same point-source forward
model. This validation establishes data-flow equivalence and execution
correctness; it does not establish physical antenna-removal equivalence for a
future solver with antenna loading, or equal reconstruction quality across
orientations of an asymmetric scene.
