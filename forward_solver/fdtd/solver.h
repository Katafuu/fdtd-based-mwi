#ifndef _FDTD_SOLVER_H
#define _FDTD_SOLVER_H

#include "../utility/fdtd-config.h"
#include "../utility/fdtd-grid1.h"
#include "../utility/fdtd-output.h"

void fdtdRun(Grid *g, const SimConfig *cfg, OutputHandler *out);
void fdtdRunTr(Grid *g, const SimConfig *cfg, OutputHandler *out);

#endif
