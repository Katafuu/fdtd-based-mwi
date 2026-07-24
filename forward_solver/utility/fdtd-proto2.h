#ifndef _FDTD_PROTO2_H
#define _FDTD_PROTO2_H

#include "fdtd-config.h"
#include "fdtd-grid1.h"
#include "fdtd-output.h"

/* Function prototypes */
void abcInit(Grid *g);
void abc(Grid *g);

void gridInit1d(Grid *g);
void gridInit(Grid *g, const char *inputDir);
Grid *gridCreate(const SimConfig *cfg);
void gridDestroy(Grid *g);
int gridInitFromConfig(Grid *g, const SimConfig *cfg);
void initTMZCavityMode(Grid *g);

void tfsfInit(Grid *g);
void tfsfUpdate(Grid *g);

void updateE2d(Grid *g);
void updateH2d(Grid *g);
void updateE2dTr(Grid *g);
void updateH2dTr(Grid *g);
void pml_update_magnetic(Grid *g);
void pml_update_electric(Grid *g);

#endif
