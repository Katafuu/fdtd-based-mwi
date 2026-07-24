#ifndef _PML_H
#define _PML_H

#include "../utility/fdtd-config.h"
#include "../utility/fdtd-grid1.h"

int pmlInit(Grid *g, const SimConfig *cfg);
void pmlDestroy(Grid *g);
void pmlUpdateH(Grid *g);
void pmlUpdateE(Grid *g);
void pmlUpdateHTr(Grid *g);
void pmlUpdateETr(Grid *g);

#endif
