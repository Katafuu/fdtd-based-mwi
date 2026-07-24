# Forward solver

This directory contains the independently buildable TMz FDTD C solver, its
MATLAB MEX interface, and file-backed examples.

## Multi-antenna source contract

Both MEX and the standalone runner use sampled sources:

- cfg.antennas.numAntennas: positive integer N.
- cfg.antennas.pos: real full numeric N x 2 matrix of one-based [x y]
  coordinates. Every coordinate must be finite, integer-valued, and inside
  the grid.
- cfg.source.samples: real full numeric N x Nt matrix. Every value must be finite.

Row a of cfg.source.samples belongs to row a of cfg.antennas.pos. At each
time step all antenna samples are added after the electric-field and PML
updates and before snapshots are recorded. Repeated antenna positions are
allowed and their samples accumulate.

MATLAB configurations may retain metadata such as cfg.source.func,
cfg.source.frequency, and cfg.source.location, but the C/MEX solver only
requires cfg.source.samples. The former sourceX, sourceY, and sourceFreq
simulation fields are not accepted.

## Canonical MATLAB array contract

The MEX gateway accepts full, real `double`, `single`, and integer arrays
and converts them internally. Grid and CPML maps use normal MATLAB `[x,y]`
orientation: `epsr`, `cond_e`, `cond_m`, and PML maps are `Nx` by `Ny`;
`murx` is `Nx` by `Ny-1`; and `mury` is `Nx-1` by `Ny`. Sparse,
logical, complex, and incorrectly shaped arrays are rejected. Coordinates,
source samples, and initial fields must also contain only finite values.

No legacy transposed-layout flag or automatic layout detection is provided.

## Standalone build and file contract

On the current Windows toolchain:

    mingw32-make -C forward_solver
    mingw32-make -C forward_solver clean
    mingw32-make -C forward_solver rebuild

The default target is tmzdemo2.exe. The CLI contract is:

    tmzdemo2.exe [input_dir] [output_csv]

With no arguments, paths resolve from the process working directory as
./results/fdtd_input and ./results/ez_field.csv.

grid_config.txt must contain num_antennas N along with the existing grid,
time-step, PML, and snapshot settings. The input directory must also contain:

- antennas.csv: N x 2, one antenna per row.
- source.csv: N x Nt, one sampled waveform per row.
- The existing material CSV files (epsr.csv, murx.csv, mury.csv, cond_e.csv,
  and cond_m.csv).
- PML profile CSV files when file-backed PML profiles are used.

The loader converts these row-wise CSV files to the same time-major internal
source layout produced by MATLAB column-major storage.

## MEX build and output

From the repository root in MATLAB:

    run('forward_solver/mex/build_fdtd_mex.m')

Call the solver with one scalar configuration struct and an optional scalar
logical time-reversal flag. `fdtd_mex(cfg)` and `fdtd_mex(cfg, false)` use
the existing forward solver; `fdtd_mex(cfg, true)` uses the loss-uncompensated
time-reversal update equations. The binary is generated under
`forward_solver/mex/` and is not committed.

The optional output flags are `cfg.returnEz`, `cfg.returnHx`, `cfg.returnHy`,
and `cfg.returnRxSignals`. Receiver output defaults to false. With the scalar
struct call

    result = fdtd_mex(cfg)

`result.rx_signals` is either empty or a real double
`numAntennas`-by-`Nt` matrix. Each column is sampled at all antenna positions
after the electric-field update and all source injections for that time step.
Receiver sampling always covers every time step and is independent of
`snapshotStart` and `snapshotStride`. The positional call remains
`[Ez,Hx,Hy] = fdtd_mex(cfg, tr)`; it does not expose receiver traces.

When `tr` is true, source samples are injected exactly as supplied; callers
must reverse receiver traces in time before passing them as
`cfg.source.samples`. The solver does not automatically negate source samples
or initial magnetic fields. Material conductivity and CPML attenuation remain
positive: this mode reverses the Yee curl signs and CPML field-coupling signs
without introducing gain. The standalone CLI remains forward-only.

An optional `cfg.init` scalar struct seeds the electromagnetic fields before
the first update. It must contain finite real numeric matrices in normal MATLAB
`[x,y]` orientation:

- `cfg.init.Ez`: `Nx` by `Ny`.
- `cfg.init.Hx`: `Nx` by `Ny-1`.
- `cfg.init.Hy`: `Nx-1` by `Ny`.

When `cfg.init` is absent, all fields start at zero. This is a field seed, not
an exact checkpoint: the time index restarts at zero and CPML auxiliary memory
starts at zero. Returned fields use canonical MATLAB orientation:

- `Ez`: `Nx` by `Ny` by `Nframes`.
- `Hx`: `Nx` by `Ny-1` by `Nframes`.
- `Hy`: `Nx-1` by `Ny` by `Nframes`.

Single-frame results remain two-dimensional.

## Examples

- examples/build_fdtd_case.m builds the shared 23-element circular-array
  scenario and activates the sampled 1.5 GHz Ricker pulse at every antenna.
- examples/run_fdtd_mex_case.m runs that scenario through MEX.
- examples/run_fdtd_file_case.m writes the file contract and runs the CLI
  solver.
- examples/build_large_rectangular_pml_case.m builds the larger
  rectangular-PML case with the same source model.
- examples/visualize_fdtd_results.m visualizes workspace or CSV output.
