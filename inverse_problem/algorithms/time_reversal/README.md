# Time reversal

The three main workflows (`tr.mlx`, `tr_rand.mlx`, and `itr_simple.m`) run the
TMz solver directly through `forward_solver/mex/fdtd_mex`. The native MEX path
is the supported solver implementation. `lib/run_fdtd_muliantenna.m` remains as
a focused MATLAB reference used by forward-solver source-matrix regression tests.

## Prerequisite

Build the MEX binary once from the repository root in MATLAB:

```matlab
run('forward_solver/mex/build_fdtd_mex.m')
```

The workflows add `forward_solver/mex` to the MATLAB path and stop with a clear
error if `fdtd_mex` is unavailable.

## Solver contract

Each workflow builds one canonical MATLAB `cfg` and passes it directly to
`fdtd_mex`. Incident and TR copies exist only where their physics or outputs differ:

- `cfg.antennas.pos` as one-based `Nant`-by-2 `[x y]` coordinates.
- `cfg.antennas.txAntennas` as the initial iterative-TR transmitter selection.
- `cfg.source.samples` as an `Nant`-by-`Nt` matrix.
- Canonical `[x,y]` material and CPML maps from the shared builders.
- `returnRxSignals = true` and `returnHx = returnHy = false`.

Forward runs request `Ez` snapshots only when forward histories are enabled.
Time-reversal runs call `fdtd_mex(tr_cfg, true)` with already reversed
receiver traces in `tr_cfg.source.samples`. They request every canonical `Ez`
frame so MATLAB can compute the maximum-power image and the maximum-magnitude
and minimum-entropy focus frames without layout conversion.

The native MEX time-reversal equations and positive material/CPML losses are
authoritative. The legacy MATLAB `lossMode` behavior and conductivity
pre-negation are not reproduced. Magnetic final states and histories are not
requested because none of the main workflows consumes them.

## Workflow behavior

- `build_test.m` creates the deterministic workspace configuration used by
  `tr.mlx`.
- `tr_rand.mlx` builds and solves a random-target practical-TR case.
- `itr_simple.m` feeds the returned TR receiver signals into the next
  simultaneous-source iteration. Run `plot_tr.m` afterward to create its three
  comparison figures when needed.
- Receiver trace processing, scattered-field subtraction, focus summaries,
  post-run figures, and optional `Ez` animations remain in MATLAB.
- Live plotting during a MEX solve is not supported; visualization occurs after
  the solver returns.
