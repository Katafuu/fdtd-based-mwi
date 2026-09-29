#ifndef FDTD_CUDA_SOLVER_H
#define FDTD_CUDA_SOLVER_H
#include "../utility/fdtd-config.h"
typedef struct {
  double *ez, *hx, *hy, *rx;
  int x, y, nx, ny, frames;
} CudaOutputs;
#ifdef __cplusplus
extern "C" {
#endif
int fdtdRunCuda(const Grid *g, const SimConfig *cfg, int reverse,
                const CudaOutputs *outputs, char *error, size_t errorSize);
#ifdef __cplusplus
}
#endif
#endif
