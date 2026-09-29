#include "mex.h"
#include "../utility/fdtd-config.h"
#include "../utility/fdtd-grid1.h"
#include "../utility/fdtd-macro-tmz.h"
#include "../utility/fdtd-output.h"
#include "../utility/fdtd-proto2.h"
#include "../fdtd/solver.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

typedef struct {
  mxArray *resultStruct;
  mxArray *ez;
  mxArray *hx;
  mxArray *hy;
  mxArray *rxSignals;
  int frame;
  int numFrames;
  int asStruct;
  int nlhs;
  mxArray **plhs;
  int returnEz;
  int returnHx;
  int returnHy;
  int returnRxSignals;
  const SimConfig *cfg;
} MexOutputBuffer;

typedef struct {
  int provided;
  const double *ez;
  const double *hx;
  const double *hy;
} MexInitialState;

typedef struct {
  double *epsr;
  double *murx;
  double *mury;
  double *condE;
  double *condM;
  double *condx;
  double *condy;
  double *kx;
  double *ky;
  double *antennaPos;
  double *sourceSamples;
  double *initEz;
  double *initHx;
  double *initHy;
} MexOwnedBuffers;

static const mxArray *getField(const mxArray *s, const char *name) {
  if (s == NULL || !mxIsStruct(s))
    return NULL;
  return mxGetField(s, 0, name);
}

static int getIntField(const mxArray *s, const char *name, int *out) {
  const mxArray *v = getField(s, name);
  if (v == NULL)
    return 0;
  if (!mxIsNumeric(v) || mxIsComplex(v) || mxGetNumberOfElements(v) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "%s must be a real scalar.", name);
  *out = (int)mxGetScalar(v);
  return 1;
}

static int getRequiredPositiveIntField(const mxArray *s, const char *name) {
  const mxArray *v = getField(s, name);
  double value;

  if (v == NULL)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "%s is required.", name);
  if (!mxIsNumeric(v) || mxIsComplex(v) || mxGetNumberOfElements(v) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "%s must be a positive integer scalar.", name);

  value = mxGetScalar(v);
  if (!isfinite(value) || value < 1.0 || value > (double)INT_MAX ||
      value != floor(value))
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "%s must be a positive integer scalar.", name);

  return (int)value;
}

static int getDoubleField(const mxArray *s, const char *name, double *out) {
  const mxArray *v = getField(s, name);
  if (v == NULL)
    return 0;
  if (!mxIsNumeric(v) || mxIsComplex(v) || mxGetNumberOfElements(v) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "%s must be a real scalar.", name);
  *out = mxGetScalar(v);
  return 1;
}

static int getBoolField(const mxArray *s, const char *name, int *out) {
  const mxArray *v = getField(s, name);
  if (v == NULL)
    return 0;
  if (mxIsLogical(v)) {
    *out = mxIsLogicalScalarTrue(v);
    return 1;
  }
  if (!mxIsNumeric(v) || mxIsComplex(v) || mxGetNumberOfElements(v) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "%s must be logical or scalar numeric.", name);
  *out = mxGetScalar(v) != 0.0;
  return 1;
}

static int parseTrFlag(const mxArray *value) {
  if (!mxIsLogical(value) || mxGetNumberOfElements(value) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidInput",
                      "tr must be a scalar MATLAB logical.");

  return mxIsLogicalScalarTrue(value);
}

static double numericValueAt(const mxArray *value, size_t index) {
  const void *data = mxGetData(value);

  switch (mxGetClassID(value)) {
  case mxDOUBLE_CLASS:
    return ((const double *)data)[index];
  case mxSINGLE_CLASS:
    return (double)((const float *)data)[index];
  case mxINT8_CLASS:
    return (double)((const int8_t *)data)[index];
  case mxUINT8_CLASS:
    return (double)((const uint8_t *)data)[index];
  case mxINT16_CLASS:
    return (double)((const int16_t *)data)[index];
  case mxUINT16_CLASS:
    return (double)((const uint16_t *)data)[index];
  case mxINT32_CLASS:
    return (double)((const int32_t *)data)[index];
  case mxUINT32_CLASS:
    return (double)((const uint32_t *)data)[index];
  case mxINT64_CLASS:
    return (double)((const int64_t *)data)[index];
  case mxUINT64_CLASS:
    return (double)((const uint64_t *)data)[index];
  default:
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "Unsupported numeric array class.");
    return 0.0;
  }
}

static void validateNumericMatrix(const mxArray *value, const char *name,
                                  mwSize expectedRows,
                                  mwSize expectedCols) {
  if (!mxIsNumeric(value) || mxIsLogical(value) || mxIsSparse(value) ||
      mxIsComplex(value))
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "%s must be a full, real numeric matrix.", name);
  if (mxGetNumberOfDimensions(value) > 2 ||
      mxGetM(value) != expectedRows || mxGetN(value) != expectedCols)
    mexErrMsgIdAndTxt(
        "fdtd_mex:InvalidConfig",
        "%s must have size %llu-by-%llu in MATLAB [x,y] orientation.",
        name, (unsigned long long)expectedRows,
        (unsigned long long)expectedCols);
}

static double *copyNumericMatrix(const mxArray *value, const char *name,
                                 mwSize rows, mwSize cols,
                                 int spatialLayout) {
  double *copy;
  size_t count = (size_t)rows * (size_t)cols;
  size_t row, col;

  validateNumericMatrix(value, name, rows, cols);
  copy = (double *)mxMalloc(count * sizeof(double));
  if (copy == NULL)
    mexErrMsgIdAndTxt("fdtd_mex:Allocation", "Could not allocate %s.", name);

  if (spatialLayout) {
    for (row = 0; row < (size_t)rows; row++)
      for (col = 0; col < (size_t)cols; col++)
        copy[row * (size_t)cols + col] =
            numericValueAt(value, row + col * (size_t)rows);
  } else {
    for (row = 0; row < count; row++)
      copy[row] = numericValueAt(value, row);
  }

  return copy;
}

static double *getNumericMatrixField(const mxArray *s, const char *fieldName,
                                     const char *qualifiedName,
                                     mwSize expectedRows,
                                     mwSize expectedCols,
                                     int required, int spatialLayout) {
  const mxArray *value = getField(s, fieldName);

  if (value == NULL) {
    if (required)
      mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "%s is required.",
                        qualifiedName);
    return NULL;
  }

  return copyNumericMatrix(value, qualifiedName, expectedRows, expectedCols,
                           spatialLayout);
}

static void releaseOwnedBuffers(MexOwnedBuffers *buffers) {
  if (buffers == NULL)
    return;

  mxFree(buffers->epsr);
  mxFree(buffers->murx);
  mxFree(buffers->mury);
  mxFree(buffers->condE);
  mxFree(buffers->condM);
  mxFree(buffers->condx);
  mxFree(buffers->condy);
  mxFree(buffers->kx);
  mxFree(buffers->ky);
  mxFree(buffers->antennaPos);
  mxFree(buffers->sourceSamples);
  mxFree(buffers->initEz);
  mxFree(buffers->initHx);
  mxFree(buffers->initHy);
  memset(buffers, 0, sizeof(*buffers));
}

static void validateFiniteArray(const double *values, size_t count,
                                const char *name) {
  size_t idx;

  for (idx = 0; idx < count; idx++) {
    if (!isfinite(values[idx]))
      mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                        "%s must contain only finite values.", name);
  }
}

static void parseMexInitialState(MexInitialState *state,
                                 MexOwnedBuffers *buffers,
                                 const mxArray *mxCfg, int nx, int ny) {
  const mxArray *init = getField(mxCfg, "init");
  size_t ezCount, hxCount, hyCount;

  memset(state, 0, sizeof(*state));
  if (init == NULL)
    return;
  if (!mxIsStruct(init) || mxGetNumberOfElements(init) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "cfg.init must be a scalar struct.");

  buffers->initEz = getNumericMatrixField(
      init, "Ez", "cfg.init.Ez", (mwSize)nx, (mwSize)ny, 1, 0);
  buffers->initHx = getNumericMatrixField(
      init, "Hx", "cfg.init.Hx", (mwSize)nx, (mwSize)(ny - 1), 1, 0);
  buffers->initHy = getNumericMatrixField(
      init, "Hy", "cfg.init.Hy", (mwSize)(nx - 1), (mwSize)ny, 1, 0);
  state->ez = buffers->initEz;
  state->hx = buffers->initHx;
  state->hy = buffers->initHy;

  ezCount = (size_t)nx * (size_t)ny;
  hxCount = (size_t)nx * (size_t)(ny - 1);
  hyCount = (size_t)(nx - 1) * (size_t)ny;
  validateFiniteArray(state->ez, ezCount, "cfg.init.Ez");
  validateFiniteArray(state->hx, hxCount, "cfg.init.Hx");
  validateFiniteArray(state->hy, hyCount, "cfg.init.Hy");
  state->provided = 1;
}

static void applyMexInitialState(Grid *g, const MexInitialState *state) {
  int x, y;

  if (!state->provided)
    return;

  for (x = 0; x < NxG(g); x++)
    for (y = 0; y < NyG(g); y++)
      EzG(g, x, y) = state->ez[(size_t)y * (size_t)NxG(g) +
                                (size_t)x];

  for (x = 0; x < NxG(g); x++)
    for (y = 0; y < NyG(g) - 1; y++)
      HxG(g, x, y) = state->hx[(size_t)y * (size_t)NxG(g) +
                                (size_t)x];

  for (x = 0; x < NxG(g) - 1; x++)
    for (y = 0; y < NyG(g); y++)
      HyG(g, x, y) = state->hy[(size_t)y * (size_t)(NxG(g) - 1) +
                                (size_t)x];
}

static void parsePmlType(PmlConfig *pml, const mxArray *typeValue) {
  char *type;

  if (typeValue == NULL)
    return;
  if (!mxIsChar(typeValue))
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "cfg.pml.type must be a string.");

  type = mxArrayToString(typeValue);
  if (type == NULL)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "Could not read cfg.pml.type.");

  if (strcmp(type, "none") == 0 || strcmp(type, "off") == 0) {
    pml->type = PML_TYPE_NONE;
    pml->enabled = 0;
  } else if (strcmp(type, "cpml") == 0 || strcmp(type, "pml") == 0) {
    pml->type = PML_TYPE_CPML;
    pml->enabled = 1;
  } else {
    mxFree(type);
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "Unknown cfg.pml.type.");
  }

  mxFree(type);
}

static void parseMexConfig(SimConfig *cfg, MexOwnedBuffers *buffers,
                           const mxArray *mxCfg) {
  const mxArray *grid;
  const mxArray *pml;
  const mxArray *antennas;
  const mxArray *source;
  char err[256];

  if (mxCfg == NULL)
    return;
  if (!mxIsStruct(mxCfg) || mxGetNumberOfElements(mxCfg) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "Input must be a scalar struct.");

  getIntField(mxCfg, "Nx", &cfg->Nx);
  getIntField(mxCfg, "Ny", &cfg->Ny);
  getIntField(mxCfg, "sizeZ", &cfg->sizeZ);
  getIntField(mxCfg, "Nt", &cfg->Nt);
  getDoubleField(mxCfg, "dx", &cfg->dx);
  getDoubleField(mxCfg, "dy", &cfg->dy);
  getDoubleField(mxCfg, "dt", &cfg->dt);
  getIntField(mxCfg, "snapshotStart", &cfg->snapshotStart);
  getIntField(mxCfg, "snapshotStride", &cfg->snapshotStride);
  getBoolField(mxCfg, "returnEz", &cfg->returnEz);
  getBoolField(mxCfg, "returnHx", &cfg->returnHx);
  getBoolField(mxCfg, "returnHy", &cfg->returnHy);
  getBoolField(mxCfg, "returnRxSignals", &cfg->returnRxSignals);

  if (cfg->Nx < 2 || cfg->Ny < 2 || cfg->Nt < 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "Nx and Ny must be at least 2 and Nt must be positive.");

  grid = getField(mxCfg, "grid");
  if (grid != NULL) {
    if (!mxIsStruct(grid) || mxGetNumberOfElements(grid) != 1)
      mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "cfg.grid must be a scalar struct.");
    buffers->epsr = getNumericMatrixField(
        grid, "epsr", "cfg.grid.epsr", (mwSize)cfg->Nx, (mwSize)cfg->Ny,
        0, 1);
    buffers->murx = getNumericMatrixField(
        grid, "murx", "cfg.grid.murx", (mwSize)cfg->Nx,
        (mwSize)(cfg->Ny - 1), 0, 1);
    buffers->mury = getNumericMatrixField(
        grid, "mury", "cfg.grid.mury", (mwSize)(cfg->Nx - 1),
        (mwSize)cfg->Ny, 0, 1);
    buffers->condE = getNumericMatrixField(
        grid, "cond_e", "cfg.grid.cond_e", (mwSize)cfg->Nx,
        (mwSize)cfg->Ny, 0, 1);
    buffers->condM = getNumericMatrixField(
        grid, "cond_m", "cfg.grid.cond_m", (mwSize)cfg->Nx,
        (mwSize)cfg->Ny, 0, 1);
    cfg->epsr = buffers->epsr;
    cfg->murx = buffers->murx;
    cfg->mury = buffers->mury;
    cfg->cond_e = buffers->condE;
    cfg->cond_m = buffers->condM;
  }

  pml = getField(mxCfg, "pml");
  if (pml != NULL) {
    if (!mxIsStruct(pml) || mxGetNumberOfElements(pml) != 1)
      mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "cfg.pml must be a scalar struct.");
    parsePmlType(&cfg->pml, getField(pml, "type"));
    getBoolField(pml, "enabled", &cfg->pml.enabled);
    getDoubleField(pml, "ax", &cfg->pml.ax);
    getDoubleField(pml, "ay", &cfg->pml.ay);
    getDoubleField(pml, "az", &cfg->pml.az);
    buffers->condx = getNumericMatrixField(
        pml, "condx", "cfg.pml.condx", (mwSize)cfg->Nx,
        (mwSize)cfg->Ny, 0, 1);
    buffers->condy = getNumericMatrixField(
        pml, "condy", "cfg.pml.condy", (mwSize)cfg->Nx,
        (mwSize)cfg->Ny, 0, 1);
    buffers->kx = getNumericMatrixField(
        pml, "kx", "cfg.pml.kx", (mwSize)cfg->Nx,
        (mwSize)cfg->Ny, 0, 1);
    buffers->ky = getNumericMatrixField(
        pml, "ky", "cfg.pml.ky", (mwSize)cfg->Nx,
        (mwSize)cfg->Ny, 0, 1);
    cfg->pml.condx = buffers->condx;
    cfg->pml.condy = buffers->condy;
    cfg->pml.kx = buffers->kx;
    cfg->pml.ky = buffers->ky;
  }

  antennas = getField(mxCfg, "antennas");
  if (antennas == NULL || !mxIsStruct(antennas) ||
      mxGetNumberOfElements(antennas) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "cfg.antennas must be a scalar struct.");
  cfg->antennas.numAntennas = getRequiredPositiveIntField(
      antennas, "numAntennas");
  buffers->antennaPos = getNumericMatrixField(
      antennas, "pos", "cfg.antennas.pos",
      (mwSize)cfg->antennas.numAntennas, 2, 1, 0);
  cfg->antennas.pos = buffers->antennaPos;

  source = getField(mxCfg, "source");
  if (source == NULL || !mxIsStruct(source) ||
      mxGetNumberOfElements(source) != 1)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig",
                      "cfg.source must be a scalar struct.");
  buffers->sourceSamples = getNumericMatrixField(
      source, "samples", "cfg.source.samples",
      (mwSize)cfg->antennas.numAntennas, (mwSize)cfg->Nt, 1, 0);
  cfg->source.samples = buffers->sourceSamples;

  if (!simConfigValidate(cfg, err, sizeof(err)))
    mexErrMsgIdAndTxt("fdtd_mex:InvalidConfig", "%s", err);
}

static int frameCount(const SimConfig *cfg) {
  int first, count;

  if (cfg->snapshotStride <= 0)
    return 1;
  if (cfg->snapshotStart >= cfg->Nt)
    return 0;

  first = cfg->snapshotStart < 0 ? 0 : cfg->snapshotStart;
  count = ((cfg->Nt - 1) - first) / cfg->snapshotStride + 1;
  return count < 0 ? 0 : count;
}

static int shouldRecordMex(const MexOutputBuffer *buf, int timeStep) {
  const SimConfig *cfg = buf->cfg;
  int start;

  if (cfg->snapshotStride <= 0)
    return timeStep == cfg->Nt - 1;

  start = cfg->snapshotStart < 0 ? 0 : cfg->snapshotStart;
  return timeStep >= start && (timeStep - start) % cfg->snapshotStride == 0;
}

static mxArray *createFieldArray(int rows, int cols, int frames) {
  if (frames <= 1)
    return mxCreateDoubleMatrix((mwSize)rows, (mwSize)cols, mxREAL);
  else {
    mwSize dims[3];
    dims[0] = (mwSize)rows;
    dims[1] = (mwSize)cols;
    dims[2] = (mwSize)frames;
    return mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  }
}

static void copyCanonicalField(mxArray *dst, const double *src,
                               int rows, int cols, int frame) {
  double *out = mxGetPr(dst);
  size_t frameOffset = (size_t)frame * (size_t)rows * (size_t)cols;
  int row, col;

  for (row = 0; row < rows; row++)
    for (col = 0; col < cols; col++)
      out[frameOffset + (size_t)row + (size_t)col * (size_t)rows] =
          src[(size_t)row * (size_t)cols + (size_t)col];
}

static int mexOutputBegin(OutputHandler *out, const Grid *grid,
                          const SimConfig *cfg) {
  MexOutputBuffer *buf = (MexOutputBuffer *)out->ctx;
  Grid *g = (Grid *)grid;
  const char *fields[] = {"Ez", "Hx", "Hy", "rx_signals", "Nx",
                          "Ny", "Nt", "snapshotStride", "dt", "dx",
                          "dy"};

  buf->cfg = cfg;
  buf->frame = 0;
  buf->numFrames = frameCount(cfg);

  if (buf->returnEz)
    buf->ez = createFieldArray(NxG(g), NyG(g), buf->numFrames);
  if (buf->returnHx)
    buf->hx = createFieldArray(NxG(g), NyG(g) - 1, buf->numFrames);
  if (buf->returnHy)
    buf->hy = createFieldArray(NxG(g) - 1, NyG(g), buf->numFrames);
  if (buf->returnRxSignals)
    buf->rxSignals = mxCreateDoubleMatrix(
        (mwSize)cfg->antennas.numAntennas, (mwSize)cfg->Nt, mxREAL);

  if (buf->asStruct && buf->nlhs > 0) {
    buf->resultStruct = mxCreateStructMatrix(1, 1, 11, fields);
    mxSetField(buf->resultStruct, 0, "Ez", buf->ez != NULL ? buf->ez : mxCreateDoubleMatrix(0, 0, mxREAL));
    mxSetField(buf->resultStruct, 0, "Hx", buf->hx != NULL ? buf->hx : mxCreateDoubleMatrix(0, 0, mxREAL));
    mxSetField(buf->resultStruct, 0, "Hy", buf->hy != NULL ? buf->hy : mxCreateDoubleMatrix(0, 0, mxREAL));
    mxSetField(buf->resultStruct, 0, "rx_signals",
               buf->rxSignals != NULL ? buf->rxSignals :
                                        mxCreateDoubleMatrix(0, 0, mxREAL));
    mxSetField(buf->resultStruct, 0, "Nx", mxCreateDoubleScalar((double)NxG(g)));
    mxSetField(buf->resultStruct, 0, "Ny", mxCreateDoubleScalar((double)NyG(g)));
    mxSetField(buf->resultStruct, 0, "Nt", mxCreateDoubleScalar((double)NtG(g)));
    mxSetField(buf->resultStruct, 0, "snapshotStride", mxCreateDoubleScalar((double)cfg->snapshotStride));
    mxSetField(buf->resultStruct, 0, "dt", mxCreateDoubleScalar(Dt));
    mxSetField(buf->resultStruct, 0, "dx", mxCreateDoubleScalar(Dx));
    mxSetField(buf->resultStruct, 0, "dy", mxCreateDoubleScalar(Dy));
    buf->plhs[0] = buf->resultStruct;
  } else {
    if (buf->nlhs >= 1)
      buf->plhs[0] = buf->ez != NULL ? buf->ez : mxCreateDoubleMatrix(0, 0, mxREAL);
    if (buf->nlhs >= 2)
      buf->plhs[1] = buf->hx != NULL ? buf->hx : mxCreateDoubleMatrix(0, 0, mxREAL);
    if (buf->nlhs >= 3)
      buf->plhs[2] = buf->hy != NULL ? buf->hy : mxCreateDoubleMatrix(0, 0, mxREAL);
  }

  return 1;
}

static int mexOutputRecord(OutputHandler *out, const Grid *grid, int timeStep) {
  MexOutputBuffer *buf = (MexOutputBuffer *)out->ctx;
  Grid *g = (Grid *)grid;

  if (buf->rxSignals != NULL) {
    double *rx = mxGetPr(buf->rxSignals);
    const int numAntennas = buf->cfg->antennas.numAntennas;
    int antenna;

    for (antenna = 0; antenna < numAntennas; antenna++) {
      const int x = (int)buf->cfg->antennas.pos[antenna] - 1;
      const int y =
          (int)buf->cfg->antennas.pos[numAntennas + antenna] - 1;
      rx[(size_t)timeStep * (size_t)numAntennas + (size_t)antenna] =
          EzG(g, x, y);
    }
  }

  if (!shouldRecordMex(buf, timeStep))
    return 1;
  if (buf->frame >= buf->numFrames)
    return 1;

  if (buf->ez != NULL)
    copyCanonicalField(buf->ez, g->ez, NxG(g), NyG(g), buf->frame);
  if (buf->hx != NULL)
    copyCanonicalField(buf->hx, g->hx, NxG(g), NyG(g) - 1, buf->frame);
  if (buf->hy != NULL)
    copyCanonicalField(buf->hy, g->hy, NxG(g) - 1, NyG(g), buf->frame);

  buf->frame++;
  return 1;
}

static int mexOutputEnd(OutputHandler *out, const Grid *grid) {
  (void)out;
  (void)grid;
  return 1;
}

int mexOutputInit(OutputHandler *out, MexOutputBuffer *buf,
                  const Grid *g, const SimConfig *cfg) {
  (void)g;
  (void)cfg;
  memset(out, 0, sizeof(*out));
  out->ctx = buf;
  out->begin = mexOutputBegin;
  out->record = mexOutputRecord;
  out->end = mexOutputEnd;
  return 1;
}

void snapshotStoreForMex(Grid *g, int timeStep, MexOutputBuffer *out) {
  OutputHandler handler;
  handler.ctx = out;
  handler.record = mexOutputRecord;
  handler.begin = NULL;
  handler.end = NULL;
  mexOutputRecord(&handler, g, timeStep);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  SimConfig cfg;
  Grid *g;
  OutputHandler out;
  MexOutputBuffer mexOut;
  MexInitialState initialState;
  MexOwnedBuffers ownedBuffers;
  int tr = 0;

  if (nrhs < 1 || nrhs > 2)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidInput",
                      "Expected a config struct and optional logical tr flag.");
  if (nlhs > 3)
    mexErrMsgIdAndTxt("fdtd_mex:InvalidOutput", "Use result = fdtd_mex(cfg, tr) or [Ez,Hx,Hy] = fdtd_mex(cfg, tr).");
  if (nrhs == 2)
    tr = parseTrFlag(prhs[1]);

  simConfigDefaults(&cfg);
  memset(&ownedBuffers, 0, sizeof(ownedBuffers));
  parseMexConfig(&cfg, &ownedBuffers, prhs[0]);
  parseMexInitialState(&initialState, &ownedBuffers, prhs[0], cfg.Nx, cfg.Ny);

  if (nlhs == 0) {
    cfg.returnEz = 0;
    cfg.returnHx = 0;
    cfg.returnHy = 0;
    cfg.returnRxSignals = 0;
  } else if (nlhs > 1) {
    cfg.returnEz = 1;
    cfg.returnHx = nlhs >= 2;
    cfg.returnHy = nlhs >= 3;
    cfg.returnRxSignals = 0;
  }

  g = gridCreate(&cfg);
  if (g == NULL) {
    releaseOwnedBuffers(&ownedBuffers);
    mexErrMsgIdAndTxt("fdtd_mex:GridInit", "Could not initialize FDTD grid.");
  }
  applyMexInitialState(g, &initialState);

  memset(&mexOut, 0, sizeof(mexOut));
  mexOut.asStruct = nlhs <= 1;
  mexOut.nlhs = nlhs;
  mexOut.plhs = plhs;
  mexOut.returnEz = cfg.returnEz;
  mexOut.returnHx = cfg.returnHx;
  mexOut.returnHy = cfg.returnHy;
  mexOut.returnRxSignals = cfg.returnRxSignals;

  mexOutputInit(&out, &mexOut, g, &cfg);
  if (tr)
    fdtdRunTr(g, &cfg, &out);
  else
    fdtdRun(g, &cfg, &out);

  gridDestroy(g);
  releaseOwnedBuffers(&ownedBuffers);
}
