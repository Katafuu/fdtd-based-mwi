# Eight-antenna, 12-scene FBTS study

## Status and physical gates

The CPU solver remains the default. The optional CUDA implementation must pass
`validateFbtsCuda` on the NVIDIA PC before it may be used for study results.
A successful build alone is not evidence of numerical correctness or acceleration.
Process workers require an installed and licensed Parallel Computing Toolbox.
CPU execution with one worker does not require that toolbox.

The study has three independently sampled grid locations for each of circle,
triangle, square and hexagon, with placement seed 42. Target permittivity is 2,
background permittivity is 45, and both are lossless and nonmagnetic. The DOI
radius is explicitly 75 mm. The square has a 20 mm edge; the other shapes have
20 mm circumdiameter. Polygon orientations are fixed. Analytical containment
and raster-mask containment are checked. Accepted geometry and RNG states are
saved; resume uses saved scene definitions.

Each scene has 29 antenna subsets, one uniformly sampled member of each
rotation/reflection orbit (seed 123). Cyclic shifts and reversal of the ordered
gap sequence are equivalent; arbitrary gap permutations are not generally
 equivalent. Pruning describes array layouts: rotating only the array relative
to a fixed off-center target need not preserve reconstruction quality.
All active transmitters are retained in ascending original index order.
The full measured RX tensor is acquired once per scene; iterative forward,
adjoint and sensitivity simulations remain independent for each reconstruction.

Preparation first tests 12, 24 and 48 ns recording windows on centered and
near-edge circles, including raw and background-subtracted signals. A doubled
window must confirm the arrivals; 96 ns is only a diagnostic extension of
48 ns. No accepted window exceeds 48 ns. Failure of this gate stops preparation.
Sensitivity factors 1, 2 and 4 are then compared on all four shapes, selecting
the largest with trace and step errors at most 5% against factor 1.

## Local MATLAB setup

Open MATLAB on the intended PC, set the current folder to the repository, then:

```matlab
addpath('buildLib','forward_solver/mex','inverse_problem/algorithms/fbts');
mex -setup C
fdtdBuildOptions = struct('mode','openmp');
run('forward_solver/mex/build_fdtd_mex.m');
fbtsHardwareInfo()
```

For a serial build, replace `openmp` with `serial`. Do not replace a build
under an existing study: checkpoints deliberately reject changed source or
solver identities. Windows needs a MATLAB-supported C/C++ compiler. The build
script supports GCC OpenMP on Linux and supported Microsoft/MinGW compiler
setups on Windows; other platforms need a configured OpenMP toolchain or the
serial build.

On the GTX 1660 Ti PC only, install/enable Parallel Computing Toolbox and a
compatible NVIDIA driver, configure a supported compiler, and also build CUDA:

```matlab
mex -setup C++
gpuDevice()
fdtdBuildOptions = struct('mode','cuda');
run('forward_solver/mex/build_fdtd_mex.m');
fbtsHardwareInfo()
```

The CUDA build currently targets Turing (`sm_75`, including the GTX 1660 Ti),
uses double precision and disables fused multiply-add contraction. It is not
a universal GPU binary. Intel graphics use the CPU backend. Other NVIDIA
architectures require a matching build and fresh validation. CUDA is never
silently replaced by CPU during a run. One worker owns the single GPU;
coarse sensitivity initially stays on the CPU. A VRAM allocation check reserves
20% of currently free device memory. Available memory and actual performance,
not the GPU model alone, determine whether CUDA is useful.

Build all desired backends before preparation, then use a new directory:

```matlab
work = fullfile(pwd,'study_8ant_12scene');
prepared = prepareFbtsStudy(work);
report = benchmarkFbtsStudy(prepared,fullfile(work,'benchmarks'),false);
% On the NVIDIA PC, use true above to validate and benchmark CUDA as well.
```

Calibration and full-length benchmarking can take many hours. Intermediate
calibration scans are saved immediately. Keep the code and binaries unchanged
for a resumed run. Do not interpret timings from another computer as the GTX
PC's timings. Review the validation and benchmark reports before the full run.

## Raw storage

The dataset uses MAT v7.3 files with large arrays as independent variables and
metadata as structs. Physical images retain `[x,y]` orientation. Plotting usually
requires a transpose; data storage does not.

| File | Contents |
| --- | --- |
| `manifest.mat` | Settings, catalog representatives, scene identities, relative references, progress and invocation timing |
| `scenes.csv` | Scene ID, shape, accepted physical center and shared-file paths |
| `cases.csv` | Antenna mappings, success/failure/pending state, timing and result hashes |
| `systems/system_<hash>.mat` | Grid, spacings, timestep, source, weights, PML, antennas, background materials, DOI mask, initial estimate |
| `scenes/scene_XXXXXX/truth.mat` | Geometry, vertices, material properties, target masks, full true material arrays, generation/RNG metadata |
| `scenes/scene_XXXXXX/acquisition.mat` | Full `[TX,RX,time]` measurements, acquisition identity and transmitter timings |
| `scenes/scene_XXXXXX/cases/<class>/result.mat` | Final estimate, DOI image trajectory, modeled RX history, coarse sensitivity history, optimizer history, last gradients/direction, mappings, shared references and raw phase durations |
| `status.mat` | Case state, UTC start/end, final saving/hashing duration and result hash |
| `checkpoint.mat`, `history_work.mat` | Restartable optimizer state and committed history while unfinished |
| `failure.mat`, `iteration_failure.mat` | Failure stage, exception and available interrupted-state diagnostics |
| `provenance/` | Source snapshot, solver identities, build/compiler and execution information |
| `benchmarks/` | Calibration, comparisons, raw samples, hardware, memory observations, execution settings and runtime estimates |

Per-transmitter forward, adjoint and sensitivity solve durations are separate
from MATLAB residual/gradient/update processing. Per-iteration history writing,
checkpoint writing, hashing and wall time are retained. A checkpoint cannot
contain the duration of its own completed write; the latest such entry may be
NaN until a later checkpoint or final result. Final result saving and hashing
are recorded in `status.mat`/`cases.csv`, outside the file being timed/hashed.
Worker time is not study elapsed time: concurrent workers overlap. UTC timestamps
provide chronology; elapsed timers provide durations. Memory snapshots identify
process high-water RSS where available; they are not a measured GPU peak.

Shared truth and acquisition arrays are stored once per scene. No Q, DSC, SSIM,
PCC or NRMSE is saved. Use `loadFbtsBatchCase` to resolve shared arrays for later
analysis. Do not edit shared files or checkpoint histories in place.

After preparation and benchmarking, the complete gated entry point is:

```matlab
[summary, dataset, manifest] = runFbtsStudy(work,false); % CPU only
% GTX PC, CPU and CUDA considered: runFbtsStudy(work,true)
```

It checks saved source/build and machine identities, runs controlled comparisons
and regression tests, selects the measured execution settings, and runs/resumes
the 348 cases. An exception during a physical or numerical gate prevents the
study from launching. Independent reconstruction failures are recorded and the
remaining cases continue. Running the same command again resumes the dataset.
A CPU-only request does not silently select a previously benchmarked CUDA run.

For Linux standalone regression tests, build a native executable from the
existing Makefile source list and set `FDTD_STANDALONE_EXECUTABLE` to its absolute
path. The supplied Windows `.exe` cannot execute natively on Linux. The test
override leaves that executable intact. The default Windows path is unchanged.

Use documentation for your installed MATLAB release when choosing a compiler
and driver: [GPU requirements](https://www.mathworks.com/help/parallel-computing/gpu-computing-requirements.html)
and [mexcuda supported compiler/build options](https://www.mathworks.com/help/parallel-computing/mexcuda.html).
The OpenMP regions contain only C array arithmetic; MATLAB API calls remain on
the calling thread, consistent with [MathWorks thread-safety requirements](https://www.mathworks.com/help/matlab/matlab_external/mex-api-is-not-thread-safe.html).
