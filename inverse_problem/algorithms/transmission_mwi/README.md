# Transmission-based MWI

This directory contains the actively retained planar and circular transmission-based microwave-imaging workflows.

## Runnables

```matlab
run('inverse_problem/algorithms/transmission_mwi/tds_2d_planar.m')
run('inverse_problem/algorithms/transmission_mwi/tds_2d_circular.m')
```

The planar workflow uses separate transmitter and receiver lists, midpoint footprint masks, pixel-wise weighted denominators, and `Delta_t = delaysteps * dt`.

The circular workflow keeps one shared antenna list, `Nant x Nant` pair matrices, self-pair exclusion, inward-boresight-relative angles, circular pair weights, and circular-specific reconstruction helpers. It remains a distinct formulation rather than being forced through the planar implementation.

Both scripts preserve their existing numerical constants, source handling, measurement ordering, reconstruction formulas, and plotting calculations. Only repository/path resolution changed.

Shared planar/circular helpers and circular-only helpers are colocated in `lib/` because they belong to this algorithm family; they are not general inverse-problem utilities.