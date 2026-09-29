# FDTD-Based MWI

This repository contains a two-dimensional TMz FDTD forward solver and MATLAB inverse-problem workflows.

## Layout

- `buildLib/` contains shared MATLAB scenario-building packages and antenna-index builders.
- `forward_solver/` contains the standalone C solver, MEX gateway, build system, and examples.
- `inverse_problem/algorithms/time_reversal/` contains the private MATLAB forward/time-reversal solver and runnable workflows.
- `inverse_problem/algorithms/transmission_mwi/` contains the planar and circular transmission-based MWI workflows.

The dependency direction is:

```text
algorithm script -> algorithm lib -> buildLib -> forward_solver MEX/standalone solver
```

The primary time-reversal, transmission-MWI, and hybrid workflows pass canonical
MATLAB configurations directly to `fdtd_mex`; legacy MATLAB solvers remain as reference implementations.

## Quick start

Build the standalone solver on Windows:

```powershell
mingw32-make -C forward_solver
```

Build the MEX gateway from MATLAB:

```matlab
run('forward_solver/mex/build_fdtd_mex.m')
```

Each runnable script resolves and adds only the shared library, its private library, and the MEX folder when needed. Do not add the entire repository with `genpath`.

Generated solver results, executables, MEX binaries, MATLAB backups, and videos are intentionally ignored.
