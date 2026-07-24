#include "fdtd-output.h"
#include "fdtd-macro-tmz.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  FILE *fp;
  char filename[512];
  int snapshotStart;
  int snapshotStride;
} FileOutputContext;

static int shouldRecord(const FileOutputContext *ctx, int timeStep) {
  if (ctx->snapshotStride <= 0)
    return timeStep == ctx->snapshotStart;

  return timeStep >= ctx->snapshotStart &&
    (timeStep - ctx->snapshotStart) % ctx->snapshotStride == 0;
}

static int fileOutputBegin(OutputHandler *out, const Grid *grid,
                           const SimConfig *cfg) {
  FileOutputContext *ctx = (FileOutputContext *)out->ctx;
  Grid *g = (Grid *)grid;

  ctx->snapshotStart = cfg->snapshotStart;
  ctx->snapshotStride = cfg->snapshotStride;

  ctx->fp = fopen(ctx->filename, "w");
  if (ctx->fp == NULL) {
    fprintf(stderr, "fileOutputBegin: could not open %s for writing.\n",
            ctx->filename);
    return 0;
  }

  fprintf(ctx->fp, "%d,%d\n", NxG(g), NyG(g));
  return 1;
}

static int fileOutputRecord(OutputHandler *out, const Grid *grid,
                            int timeStep) {
  FileOutputContext *ctx = (FileOutputContext *)out->ctx;
  Grid *g = (Grid *)grid;
  int mm, nn;

  if (ctx->fp == NULL)
    return 0;

  if (!shouldRecord(ctx, timeStep))
    return 1;

  for (mm = 0; mm < NxG(g); mm++) {
    for (nn = 0; nn < NyG(g); nn++) {
      fprintf(ctx->fp, "%.17g", Ez(mm, nn));

      if (!(mm == NxG(g) - 1 && nn == NyG(g) - 1))
        fprintf(ctx->fp, ",");
    }
  }

  fprintf(ctx->fp, "\n");
  return 1;
}

static int fileOutputEnd(OutputHandler *out, const Grid *grid) {
  FileOutputContext *ctx = (FileOutputContext *)out->ctx;
  (void)grid;

  if (ctx->fp != NULL) {
    fclose(ctx->fp);
    ctx->fp = NULL;
  }

  return 1;
}

int fileOutputInit(OutputHandler *out, const char *filename) {
  FileOutputContext *ctx;

  memset(out, 0, sizeof(*out));
  ctx = (FileOutputContext *)calloc(1, sizeof(*ctx));
  if (ctx == NULL) {
    perror("fileOutputInit");
    return 0;
  }

  snprintf(ctx->filename, sizeof(ctx->filename), "%s", filename);
  out->ctx = ctx;
  out->begin = fileOutputBegin;
  out->record = fileOutputRecord;
  out->end = fileOutputEnd;

  return 1;
}

void fileOutputDestroy(OutputHandler *out) {
  FileOutputContext *ctx;

  if (out == NULL || out->ctx == NULL)
    return;

  ctx = (FileOutputContext *)out->ctx;
  if (ctx->fp != NULL)
    fclose(ctx->fp);
  free(ctx);
  memset(out, 0, sizeof(*out));
}
