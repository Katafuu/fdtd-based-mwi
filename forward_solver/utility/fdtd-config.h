#ifndef _FDTD_CONFIG_H
#define _FDTD_CONFIG_H

#include <stddef.h>
#include "fdtd-grid1.h"

typedef struct {
  PmlType type;
  int enabled;
  double ax, ay, az;

  const double *kx;
  const double *ky;
  const double *condx;
  const double *condy;
} PmlConfig;

typedef struct {
  int numAntennas;
  const double *pos; /* numAntennas-by-2, MATLAB column-major layout */
} AntennaConfig;

typedef struct {
  const double *samples; /* numAntennas-by-Nt, MATLAB column-major layout */
} SourceConfig;

typedef struct {
  int Nx;
  int Ny;
  int sizeZ;
  int Nt;

  double dx;
  double dy;
  double dt;

  int snapshotStart;
  int snapshotStride;
  int returnEz;
  int returnHx;
  int returnHy;
  int returnRxSignals;

  const double *epsr;
  const double *murx;
  const double *mury;
  const double *cond_e;
  const double *cond_m;

  PmlConfig pml;
  AntennaConfig antennas;
  SourceConfig source;

  int useFileInput;
  char inputDir[512];

  /* Non-NULL only when the standalone file loader owns these buffers. */
  double *ownedAntennaPos;
  double *ownedSourceSamples;
} SimConfig;

void simConfigDefaults(SimConfig *cfg);
void simConfigRelease(SimConfig *cfg);
int simConfigValidate(const SimConfig *cfg, char *err, size_t errSize);
int gridLoadStandaloneConfig(SimConfig *cfg, const char *inputDir);

#endif
