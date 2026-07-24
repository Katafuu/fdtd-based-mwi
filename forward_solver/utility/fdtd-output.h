#ifndef _FDTD_OUTPUT_H
#define _FDTD_OUTPUT_H

#include "fdtd-config.h"
#include "fdtd-grid1.h"

typedef struct OutputHandler OutputHandler;

struct OutputHandler {
  void *ctx;
  int (*begin)(OutputHandler *out, const Grid *g, const SimConfig *cfg);
  int (*record)(OutputHandler *out, const Grid *g, int timeStep);
  int (*end)(OutputHandler *out, const Grid *g);
};

int fileOutputInit(OutputHandler *out, const char *filename);
void fileOutputDestroy(OutputHandler *out);

#endif
