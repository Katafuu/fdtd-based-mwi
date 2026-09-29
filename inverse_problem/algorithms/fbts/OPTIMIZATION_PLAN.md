# FBTS CPU/GPU optimization proposal

Status: source and existing pilot timings reviewed on 2026-09-25. This is a
proposal; no additional solver acceleration has been implemented or benchmarked.
Target server specifications are not yet known. The existing process-worker
mode still needs validation on an installation with Parallel Computing Toolbox.

## Evidence and priorities

For the eight-active-antenna, two-iteration 400×400 pilot recorded in
VALIDATION.md, the inversion loop took 86.56 s. Fine forward calls took 38.60 s,
fine adjoint calls 38.32 s, and both coarse sensitivity phases together 0.91 s.
Thus fine forward/adjoint calls account for about 89% of the loop; coarse
sensitivity accounts for about 1%. MEX call timings include setup, field output
allocation/copies and computation; these measurements do not isolate kernel
arithmetic. Instrument those components before attributing the bottleneck to
CPU arithmetic or memory bandwidth.

Current relevant code:

- `runFbtsBatch.m`: process `parfeval` workers already exist, but execute in
  bounded groups with a barrier between groups. Results are collected in
  submission order with `fetchOutputs`.
- `runFbts.m`: transmitters are sequential. Each forward/adjoint call returns
  a full Ez movie; MATLAB then slices a DOI bounding rectangle and constructs
  additional arrays for time reversal, differences, averaging and products.
- `forward_solver/mex/fdtd_mex.c`: allocates full requested histories and copies
  every frame from the C layout into MATLAB's layout.
- `solver.c`, `solver_tr.c`, field-update and PML files: serial time loop and
  spatial loops; no OpenMP or GPU kernels. PML corrections traverse the full
  domain. Each solve rebuilds grid/PML state and prints diagnostics once.
- `build_fdtd_mex.m`: ordinary CPU MEX compilation, no explicit OpenMP/CUDA
  build mode. The standalone Makefile has Windows-specific cleanup commands;
  it is not the portable MATLAB build entry point.

## Stage 1: improve batch throughput

Keep independent reconstruction cases as the primary unit of parallel work.
There are 223 dihedral representatives for twelve antennas, so this workload
usually supplies enough jobs without parallelizing inside an inversion.

Replace group barriers with a rolling queue: submit at most W cases, use
`fetchNext` to receive a finished case, commit its status, then immediately
submit the next. Preserve canonical catalog/CSV order, saved representatives,
the full-array-first submission rule, and ascending TX order within each case.
Completion order can vary. A slow job must not leave other workers idle.

Use `parallel.pool.Constant` or a validated worker-local scene cache to transfer
immutable scene/acquisition data once per worker per scene. This does not create
shared physical memory between worker processes. Maintain one writer per case
and one coordinator for shared files. Document filesystem availability if
extending beyond a single server.

Start with one C solver thread per MATLAB worker. Benchmark W=1,2,4,... within
the machine's core/RAM budget, monitoring throughput, memory bandwidth, resident
memory and storage contention. The measured 2.44 GiB peak for one entire pilot
MATLAB process is a starting observation, not a guaranteed per-worker bound.

If one or very few cases must finish quickly, an alternative is transmitter
parallelism: independent forward→residual→adjoint tasks produce per-TX gradients,
followed by a deterministic ordered reduction; after the direction is known,
independent sensitivity pairs produce step coefficients. Do not nest this inside
case parallelism by default. FBTS iterations and FDTD time steps remain ordered
because later states depend on earlier states. Initial acquisition transmitters
can also run independently, but that once-per-scene work has lower priority.

MathWorks documents asynchronous jobs in [parfeval](https://www.mathworks.com/help/parallel-computing/parallel.pool.parfeval.html)
and worker-local reusable values in [parallel.pool.Constant](https://www.mathworks.com/help/parallel-computing/parallel.pool.constant.html).

## Stage 2: reduce allocations and field-history traffic

First add an optional observation region to the MEX output interface. Continue
propagating fields on the **entire simulation domain**, including PML, but return
only the DOI bounding rectangle or exact DOI samples needed by the gradient.
Retain full-field mode for existing callers and regression comparisons.

Then fuse gradient accumulation into an adjoint-output callback in C: retain
the needed forward DOI history, pair adjoint samples with their exact forward
time indices, and accumulate the trapezoidal product used by the current MATLAB
gradient. This avoids storing a second full adjoint movie and the MATLAB
`flip`, `diff`, averaged-field and product temporaries. Preserve staggering,
signs, endpoints, dt and c0 exactly. Moving/reordering reductions may change
roundoff, so compare gradients and complete trajectories, not just final images.

One 400×400×600 double history is 768,000,000 bytes. An observation-region
history scales with the fraction of spatial cells retained. This change can
increase the number of cases that fit in RAM as well as shorten copies.

After profiling, consider reusable per-worker solver contexts for allocation,
antenna indexing and invariant PML coefficients. Reset all fields and CPML
memory on every independent solve. Recompute coefficients depending on the
current material map. Make diagnostic logging configurable. Any sparse-PML
optimization must prove that omitted coefficients/memory cannot contribute;
do not assume every imported PML profile is a conventional rectangular strip.

## Stage 3: accelerate the CPU solver

Provide separate serial and OpenMP build modes. Parallelize spatial update and
copy loops, preserving synchronization between H updates, H-PML corrections,
E updates, E-PML corrections, source injection and sampling. Keep time steps
sequential. For small coarse grids, avoid threading overhead with a measured
size threshold. Consider a persistent thread team after basic correctness.

Keep MATLAB allocation, error reporting, printing and other MEX/Matrix API
calls on the calling thread. Worker threads operate only on prepared C buffers.
The existing output callbacks call mx APIs, so they cannot simply be executed
inside arbitrary OpenMP worker threads. See [MEX API thread safety](https://www.mathworks.com/help/matlab/matlab_external/mex-api-is-not-thread-safe.html).

Inspect compiler vectorization reports and actual MEX optimization flags.
Preserve the contiguous inner loop, evaluate blocked layout conversion, and
test safe alias/alignment declarations. Precomputing reciprocals or changing
floating-point contraction can alter rounding; retain the original equations
as a reference. Do not enable fast-math or single precision by default.

Tune W worker processes and T solver threads together, initially keeping W×T
within allocated physical cores and accounting for MATLAB's own threads.
Benchmark many single-thread cases against fewer multithread cases. Linear
scaling is not guaranteed: stencil traffic can saturate memory bandwidth.
OpenMP itself does not require Parallel Computing Toolbox; MATLAB process pools
do. Compiler/runtime availability and compatibility still need verification.
[OpenMP](https://www.openmp.org/about/) supports shared-memory C/C++ parallelism.

## Stage 4: optional NVIDIA GPU backend

Add an explicit CUDA backend behind the same solver interface; passing gpuArray
to the existing CPU-only MEX will not port its C loops. MATLAB supports custom
CUDA MEX via [mexcuda](https://www.mathworks.com/help/parallel-computing/run-mex-functions-containing-cuda-code.html).

Keep fields, coefficients, CPML memory, source samples and necessary forward
history on the device across time steps. Port both ordinary and time-reversed
updates, source injection, sampling and gradient accumulation. Transfer small
RX/gradient/result arrays back when needed, avoiding per-step full-field
transfers. This follows NVIDIA's [data-transfer guidance](https://developer.nvidia.com/blog/how-optimize-data-transfers-cuda-cc/).

Start with double precision: finite-difference sensitivity and iterative
gradients need validation before changing precision. Compare GPU double
throughput, memory bandwidth and VRAM requirements for the actual device;
gaming specifications alone do not predict this solver's performance. Small
coarse solves may remain faster on CPU. Measure the full inversion, including
transfers and launch overhead, rather than isolated kernel speed.

Start with one process worker per GPU, assigning a device explicitly. Multiple
workers fighting for one GPU can duplicate memory and reduce throughput.
Batching independent TXs on a GPU is a later utilization experiment, bounded by
VRAM. MathWorks provides [multi-GPU worker patterns](https://www.mathworks.com/help/parallel-computing/run-matlab-functions-on-multiple-gpus.html).

## Portability contract

| Machine | Proposed supported path | Conditions |
| --- | --- | --- |
| Supported Intel/AMD x86-64 CPU, no compatible GPU | Serial or OpenMP CPU solver; optional MATLAB process pool | Matching MATLAB OS/release, compiled MEX, compatible compiler/OpenMP runtime; toolbox for pool |
| Supported Intel/AMD CPU plus compatible NVIDIA GPU | CPU fallback or CUDA backend | MATLAB-supported GPU compute capability, driver, VRAM and CUDA/host compiler combination |
| Intel CPU with Intel integrated/Arc/data-center GPU | CPU backend; Intel GPU unused by CUDA | Intel GPU acceleration requires another backend |
| Intel GPU acceleration required | Custom SYCL/oneAPI or OpenCL library called through a host MEX | Separate development, runtime/toolchain/device support and regression testing |

CPU vendor and GPU vendor are independent choices. An Intel CPU with an NVIDIA
GPU is a normal CUDA host configuration. No design should promise every NVIDIA
GPU, every Intel system, or one universal compiled binary.

MATLAB's native GPU path requires supported NVIDIA hardware; capability ranges
depend on the MATLAB release. Check the target release's [GPU requirements](https://www.mathworks.com/help/parallel-computing/gpu-computing-requirements.html)
and [supported compilers](https://www.mathworks.com/support/requirements/supported-compilers.html).
Windows and Linux MEX binaries are different; rebuild or distribute tested
platform-specific binaries as described in [MEX platform compatibility](https://www.mathworks.com/help/matlab/matlab_external/platform-compatibility.html).
Avoid distributing host-specific AVX/native builds to less capable CPUs.

SYCL can provide source portability across supported vendor backends, but does
not make Intel devices into MATLAB gpuArray devices or guarantee universal
binary/performance portability. Intel describes [cross-vendor SYCL tooling](https://software.seek.intel.com/oneapi-ws-SYCL-portability).
Choose that larger maintenance commitment only if Intel GPU support is an
actual requirement. Otherwise use a broadly compatible CPU fallback and an
optional NVIDIA CUDA backend.

## Acceptance and reproducibility

For each stage measure end-to-end cases/hour, individual-case latency, peak
RAM/VRAM and output/checkpoint overhead. Reuse identical saved scenes and
representatives across variants. Test full, reduced and singleton arrays;
forward/adjoint traces and gradients; sensitivity factors 1 and 4; trajectory,
projection/restart behavior; and interrupted recovery. Require numerical
tolerances justified against double-precision reference results, rather than
assuming bitwise equality across compilers/GPUs/reduction orders.

Record resolved backend, device, precision, compiler/flags, runtime/driver,
thread/worker counts and build hash. Preserve CPU fallback and reproducible
explicit backend selection. Select a backend before acquisition, report the
selection, and never silently switch numerical backends during resume.
Existing acquisition/checkpoint identity checks reject changed solver binaries;
start a new study or implement an explicit validated migration, not a bypass.

Recommended order: rolling case scheduling → observation-region output and
fused gradient → OpenMP/SIMD tuning → CUDA if supported target hardware warrants
it. The CPU/memory work remains useful whether or not a GPU is available.
