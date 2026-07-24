# Inverse problems

Active inverse algorithms live under `algorithms/`:

- `time_reversal/` uses a private MATLAB forward/time-reversal solver.
- `transmission_mwi/` uses the forward solver MEX interface and retains separate planar and circular formulations.

General scenario construction belongs in the repository-level `buildLib/`. A shared `inverse_problem/lib/` is intentionally absent because no helper is currently proven to be shared by multiple inverse algorithms with the same semantics. Create it only when such a dependency exists.

Algorithm scripts should add only root `buildLib`, their own `lib`, and `forward_solver/mex` when they actually use MEX.