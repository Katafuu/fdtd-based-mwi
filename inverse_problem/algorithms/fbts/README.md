# FBTS reconstructions and reproducible antenna batches

For the new 15 cm DOI / permittivity-45, 12-scene study, see
[STUDY_GUIDE.md](STUDY_GUIDE.md). It documents the gated preparation, CPU/OpenMP
and optional CUDA builds, benchmark workflow, and complete raw-data layout.
GPU and process-pool validation require the corresponding hardware/toolbox;
the old pilot timings below do not estimate this new study.

Build the FDTD MEX gateway first. From this directory, the original standalone
workflow remains available:

```matlab
build_cfg;
fbts_demo;
```

The batch default now keeps one reproducibly random representative per circular
layout class under **rotations and reflections**. These are binary bracelets:
29 configurations for 8 antennas, 223 for 12, including the full array and
excluding the no-active-antenna configuration. Arbitrary gap permutations are
not merged. For example, `[1 2 3]` and `[3 2 1]` are equivalent, whereas
`[1 2 3 4]` and `[1 3 2 4]` are not.

Each scene gets one full initial RX scan, consisting of N individual transmitter
solves. Every selected configuration uses `rx_full(A,A,:)`, then starts its own
FBTS inversion from the background. All active antennas transmit and receive,
including self-receiver channels. Transmission starts with the lowest original
active antenna index. Forward, adjoint, and both sensitivity solves remain
independent for every configuration and iteration.

## Run a batch

```matlab
fbtsBuildOptions = struct('seed', 42, 'numAntennas', 12);
build_cfg;
fbtsOptions = struct('numIterations', 15, ...
    'sensitivityDownsampleFactor', 4, ...
    'outputDirectory', "results/antenna_batch");
fbtsBatchOptions = struct('seed', 123, 'symmetry', 'dihedral');
fbts_batch;
```

The script returns `batchSummary`, `batchOutputDirectory`, and `batchManifest`.
The input `cfg` is unchanged. A new output directory must be empty. If omitted,
the output directory is the next `figs/batch_XXXX` folder. The default batch is
headless; no figures are created. The standalone demo still saves its legacy
`cfg` / `results` MAT file and four figures, and now writes its MAT data before
plotting so a plotting error does not discard numerical output.

The function API supports one cfg, a struct array, or a cell array of scenes:

```matlab
scene1 = buildFbtsConfig(struct('seed', 10));
scene2 = buildFbtsConfig(struct('seed', 11));
[summary, directory, manifest] = runFbtsBatch( ...
    {scene1, scene2}, fbtsOptions, fbtsBatchOptions);
```

Alternatively set `fbtsScenes = {scene1, scene2}` before running the script.
A scene can contain multiple target objects. Generation records its exact seed,
RNG state, options, object geometry, physical coordinates, and actual sampled
masks. Explicit builder seeds leave the caller's RNG state unchanged. If the
seed is omitted, the builder draws and records a seed from the caller's RNG.
Imported configurations must supply target definitions and, for nonstandard
rasterization, exact per-object simulation-grid masks. The saved masks and
material properties must reproduce the actual permittivity map.

`buildFbtsConfig` accepts `seed`, `targetOptions`, `Nx`, `Ny`, `dx`, `dy`,
`deltaF`, `numAntennas`, `pmlThicknessRatio`, `pmlPadding`, and `focusPadding`.
`targetOptions` supports the existing generator fields: `numTargetsRange`,
`allowedShapes`, `radiusRange`, `sideRange`, `epsrRange`, `condRange`,
`maxAttempts`, and `allowOverlap`. Batch inversion currently validates lossless,
nonmagnetic materials (zero conductivity and relative permeability one), a
nonempty DOI, zero initial fields, and full-resolution temporal snapshots.
It estimates permittivity only; it does not reconstruct conductivity.

## Batch controls

| Field | Default | Meaning |
| --- | --- | --- |
| `symmetry` | `'dihedral'` | Rotations plus reflections. `'rotation'` retains reflections separately; `'none'` enumerates all subsets. |
| `seed` | `0` | Representative-selection master seed, separate from target generation. |
| `disabledCounts` | `[]` | All counts `0:N-1`; specify e.g. `[0 1 2]` to restrict the study. |
| `maxCandidates` | `2e6` | Fail before enumerating more candidates; raise explicitly or restrict disabled counts. |
| `dryRun` | `false` | Return the catalog and resource estimates without FDTD or output files. |
| `resume` | `false` | Resume the explicit output directory after validating its manifest and hashes. |
| `checkpointEvery` | `1` | Commit a checkpoint every this many complete iterations, and at the final iteration. |
| `saveFigures` | `false` | Export truth, reconstruction, and objective plots after numerical data are saved. |
| `maxWorkers` | `1` | Optional process workers; requires Parallel Computing Toolbox above one. |
| `maxCases` | `Inf` | Limit attempted reconstruction cases in this invocation; remaining cases stay queued. Useful for pilots. |

Representatives are selected uniformly from distinct transforms in each orbit,
using a local stream derived from seed, stable scene ID, and class key. The
choices are independent of execution order and preserved in the manifest.
Cases run with the full array first, then increasing disabled count. Symmetry
pruning validates the nominal uniform circular ordering; use `'none'` for
other arrays. The current index-space circle requires `dx == dy` for this
physical circular-array assumption.

A layout class does **not** imply identical reconstruction quality for its
members around a fixed off-center/asymmetric target. Grid rounding and PML also
limit physical symmetry. Store and use the selected orientation and orbit size
when interpreting results. An unweighted class average differs from an average
over all original subsets, and a single representative does not measure
within-class variability.

## Dry run, pilot, and resume

```matlab
fbtsBatchOptions.dryRun = true;
fbts_batch;                         % prints fine/coarse call and storage estimates
fbtsBatchOptions.dryRun = false;
fbtsBatchOptions.maxCases = 1;
fbts_batch;                         % one completed case, rest queued
fbtsBatchOptions.resume = true;
fbtsBatchOptions.maxCases = Inf;
fbts_batch;                         % verifies and skips completed cases
```

The last call must use the same original scenes, selection settings, numerical
options, and source/MEX identity. A changed numerical iteration limit starts a
new study rather than extending an old checkpoint silently. Operational controls
such as worker count, plot output, or checkpoint interval can change on resume.
For an automatically allocated output directory, set
`fbtsOptions.outputDirectory = batchOutputDirectory` before resuming.

A completed iteration checkpoint stores the current estimate, previous projected
gradient and direction, projection/restart state, scalar histories, solver
runtimes, and checksums of committed image/RX history. A partial iteration is
redone. Completed results and shared files are checked before reuse. Corrupt or
incompatible files are rejected. A coordinator file lock prevents simultaneous
batch coordinators from writing the same directory and is released when the
process exits. Run and checkpoint files are written to temporary files in the
same directory and then renamed; actual durability depends on the filesystem.

Each case writes its own status. The coordinator refreshes `manifest.mat` and
`cases.csv` after each case; resume recovers statuses even if the coordinator was
interrupted before updating the CSV. A reconstruction failure is recorded and
later cases continue. An acquisition failure blocks that scene's pending cases
and leaves unrelated scenes runnable. A plotting failure is reported separately
and does not turn a saved numerical result into a failed reconstruction.

For multiple workers, measure actual peak memory first. One `400 x 400 x 600`
double Ez history alone uses about 768 MB, and fields/temporary arrays coexist.
Dry-run memory estimates are individual-array sizes, not total process RAM.
The optional worker mode uses processes, caps dispatched cases at `maxWorkers`,
and avoids shared MAT writers and nested transmitter parallelism.

## Dataset schema (version 1)

```text
batch_XXXX/
  manifest.mat                     # catalogs, options, identities, resource estimates, statuses
  cases.csv                        # portable relative paths, status and timing summary
  provenance/metadata.mat          # code/build/platform identity
  provenance/source/...            # relevant MATLAB/C sources, including untracked files
  systems/system_<hash>.mat        # one shared system definition per unique system
  scenes/scene_000001/
    truth.mat                      # target definitions, masks and actual material maps
    acquisition.mat                # one full initial RX tensor for this scene
    acquisition_failure.mat        # only if acquiring/validating this scan fails
    cases/n..._dihedral_.../
      result.mat                   # reconstructed maps, optimizer and RX histories
      status.mat                   # authoritative result hash, timing and completion status
      checkpoint.mat               # last committed optimizer state during unfinished runs
      history_work.mat             # streamed image/RX arrays during unfinished runs
      failure.mat                  # diagnostic for a failed case
      iteration_failure.mat        # phase/TX, partial timings and optimizer diagnostics
      figures/                     # optional, generated after result commit
```

Common arrays are not copied into every result. Shared file references are
relative to the batch root and protected by hashes. Move/copy the whole batch
directory to keep a portable dataset. Completed cases remove working history
and checkpoint files after committing the result and status.

| File | Variables and saved contents |
| --- | --- |
| `manifest.mat` | `manifest`: schema and request identity, resolved options, full class/representative catalogs including RNG state, scene identities/references, source identity, resource estimates, status rows and invocation wall times. |
| System MAT | `model`: executable cfg template, dimensions/spacings/time settings/constants, grid metadata, PML profiles, source definition, original antenna geometry/order, exact physical positions, coarse sensitivity geometry/time/PML/interpolation metadata, estimated/prescribed-property declarations. Independent `x_m`, `y_m`, `time_s`, `sourcePulse`, `timeWeight`, all background material maps, `epsr_initial`, `doi_mask`, `valid_domain_mask`, `doi_linear_indices`, and content checksum. |
| `truth.mat` | `scene` and `targets` structs: IDs, generation state/settings, shape/index/physical geometry, material properties, layering and sampling provenance. Independent `epsr_true`, `cond_e_true`, `cond_m_true`, `murx_true`, `mury_true`, `target_masks`, `target_union_mask`, `label_image`, and content checksum. |
| `acquisition.mat` | Independent `rx_full` double `[original TX, original RX, time]`. `measurement` metadata: axes, coordinates, source/time and solver identity, total-Ez convention, noiseless setting, completeness, fingerprints. `measurementTiming`: raw per-TX durations and acquisition wall time. |
| `result.mat` | `caseInfo`: original active/disabled indices, local mapping, first TX, order, class/gaps/orbit/selected transform/RNG, shared references and hashes. `parameters`: resolved FBTS settings and explicit hard-coded numerical rules. Independent `epsr_est`, `epsr_history_doi`, `rx_model_history`, `rx_sensitivity_coarse_history`, last DOI gradient/direction arrays. `history`: cost/per-TX costs, alpha/unconstrained alpha, beta, restarts/reasons, projections/cell counts, derivative, step_a/step_q, finite-difference h, sensitivity dt/Nt, gradient/update norms, evaluated/updated state IDs, channel count, measurement energy and stopping status. `timing`: raw solve/iteration/setup/reconstruction durations. |
| `status.mat` | Completion/failure state, timestamps, completed updates, case identifiers, result hash, solver statistics, setup/reconstruction/save/plot/total durations, plot-only errors and cache hit. Saving duration includes result commit/hash and excludes writing status itself. |
| `checkpoint.mat` | Numerical optimizer state, completed iteration, scalar histories, raw timings, code/input identity and checkpoint/history checksums. |
| `failure.mat` | Error identifier/message/stack, failing stage, available status/timings, completed-iteration count and checkpoint reference. |
| `iteration_failure.mat` | When an iteration fails: phase, iteration/local TX, last completed update, current estimate, scalar histories, partial raw solve timings, and exception/stack. Unperformed solves have NaN timing. Successful recovery removes this diagnostic. |
| Provenance | Source files and hashes, source/MEX identity, Git revision and dirty/untracked status, MATLAB/platform and available hardware facts. Unknown original compiler/build flags are explicitly marked unknown. |

Matrices retain MATLAB `[x,y]` orientation; figures transpose only for display.
Values remain double precision and unnormalized. Material conductivity units
are S/m, physical coordinates are meters, and time is seconds. RX amplitudes
retain solver units; no experimental amplitude calibration is implied.

If `I` iterations complete, `epsr_history_doi` has `I+1` columns, starting with
the initial image. `rx_model_history(:,:,:,j)`, `history.cost(j)`, and the gradient
at iteration j correspond to **estimate state j-1**, before its update. The
final image is state I. There is no extra final forward evaluation hidden in
the saved last cost. Outside-DOI pixels are prescribed background; use
`doi_linear_indices` to expand the compact history onto the full grid.

No Q, DSC, SSIM, PCC, NRMSE, or relative-image-error metrics are computed/saved
by the batch. The raw maps and iteration image history allow later evaluation;
object masks support DSC after choosing a segmentation rule. Choose the ROI,
SSIM window/dynamic range, normalization, and the definition of Q in the later
analysis. `history.step_q` is an optimizer coefficient, not image-quality Q.
Full space-time field movies are not saved by default.

Use the loader to assemble a case without manually resolving references:

```matlab
caseFile = fullfile(batchOutputDirectory, batchSummary.result_file(1));
data = loadFbtsBatchCase(caseFile);   % verifies hashes by default
truth = data.truth.epsr_true;
reconstructed = data.result.epsr_est;
doi = data.system.doi_mask;
rx = data.rx_selected;
```

`data.cfg` is the reconstructed reduced configuration and `data.measurements`
is the shared scan, suitable for `runFbts(data.cfg, options, data.measurements)`.
Loading for analysis requires only the FBTS helper files and MATLAB; the FDTD
binary and `buildLib` are not required. Rerunning the inversion requires the
solver and validates its identity against the saved acquisition.
The default batch uses the same forward model for synthetic truth and inversion;
that modeling assumption should be considered when interpreting image metrics.

## Validation

From the repository root:

```matlab
addpath('buildLib', 'forward_solver/mex', 'inverse_problem/algorithms/fbts');
r = runtests('inverse_problem/algorithms/fbts/tests');
assert(~any([r.Failed]));            % inspect Incomplete for skipped tests
```

The tests exercise real MEX scans, fresh versus cached reconstructions, symmetry
counts/orbits, checksum rejection, interrupted conjugate-gradient recovery,
portable saved data, multiple scenes, deterministic generation, and case failure
isolation. The process-worker comparison is skipped when Parallel Computing
Toolbox is unavailable. [VALIDATION.md](VALIDATION.md) records the regression
results and full-resolution pilot. [BATCH_IMPLEMENTATION_PLAN.md](BATCH_IMPLEMENTATION_PLAN.md)
records the staged design and its reflection amendment.
