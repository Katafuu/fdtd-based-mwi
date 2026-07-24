#include "fdtd-config.h"
#include "fdtd-io.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void setError(char *err, size_t errSize, const char *msg) {
  if (err != NULL && errSize > 0)
    snprintf(err, errSize, "%s", msg);
}

void simConfigDefaults(SimConfig *cfg) {
  memset(cfg, 0, sizeof(*cfg));

  cfg->Nx = 250;
  cfg->Ny = 250;
  cfg->sizeZ = 1;
  cfg->Nt = 400;
  cfg->dx = 0.01;
  cfg->dy = cfg->dx;
  cfg->dt = 0.0;

  cfg->snapshotStart = 0;
  cfg->snapshotStride = 0; /* generic default: final field only */
  cfg->returnEz = 1;
  cfg->returnHx = 0;
  cfg->returnHy = 0;
  cfg->returnRxSignals = 0;

  cfg->pml.type = PML_TYPE_NONE;
  cfg->pml.enabled = 0;
  cfg->pml.ax = 1.0;
  cfg->pml.ay = 1.0;
  cfg->pml.az = 1.0;
}

void simConfigRelease(SimConfig *cfg) {
  if (cfg == NULL)
    return;

  free(cfg->ownedAntennaPos);
  free(cfg->ownedSourceSamples);
  cfg->ownedAntennaPos = NULL;
  cfg->ownedSourceSamples = NULL;
  cfg->antennas.pos = NULL;
  cfg->source.samples = NULL;
}

int simConfigValidate(const SimConfig *cfg, char *err, size_t errSize) {
  if (cfg->Nx < 2 || cfg->Ny < 2) {
    setError(err, errSize, "Nx and Ny must both be at least 2.");
    return 0;
  }

  if (cfg->Nt < 1) {
    setError(err, errSize, "Nt must be positive.");
    return 0;
  }

  if (cfg->dx <= 0.0 || cfg->dy <= 0.0) {
    setError(err, errSize, "dx and dy must be positive.");
    return 0;
  }

  if (cfg->dt < 0.0) {
    setError(err, errSize, "dt must be nonnegative; use 0 to auto-compute it.");
    return 0;
  }

  if (cfg->antennas.numAntennas < 1) {
    setError(err, errSize, "cfg.antennas.numAntennas must be positive.");
    return 0;
  }

  if (cfg->antennas.pos == NULL) {
    setError(err, errSize, "cfg.antennas.pos is required.");
    return 0;
  }

  if (cfg->source.samples == NULL) {
    setError(err, errSize, "cfg.source.samples is required.");
    return 0;
  }

  {
    int antenna;
    const int numAntennas = cfg->antennas.numAntennas;

    for (antenna = 0; antenna < numAntennas; antenna++) {
      const double x = cfg->antennas.pos[antenna];
      const double y = cfg->antennas.pos[numAntennas + antenna];

      if (!isfinite(x) || !isfinite(y) || x != floor(x) || y != floor(y)) {
        setError(err, errSize,
                 "cfg.antennas.pos must contain finite integer coordinates.");
        return 0;
      }
      if (x < 1.0 || x > (double)cfg->Nx ||
          y < 1.0 || y > (double)cfg->Ny) {
        setError(err, errSize,
                 "cfg.antennas.pos contains a coordinate outside the grid.");
        return 0;
      }
    }
  }

  {
    size_t sample;
    const size_t sampleCount = (size_t)cfg->antennas.numAntennas *
                               (size_t)cfg->Nt;

    for (sample = 0; sample < sampleCount; sample++) {
      if (!isfinite(cfg->source.samples[sample])) {
        setError(err, errSize,
                 "cfg.source.samples must contain only finite values.");
        return 0;
      }
    }
  }

  if (cfg->pml.type != PML_TYPE_NONE && cfg->pml.type != PML_TYPE_CPML) {
    setError(err, errSize, "unknown PML type.");
    return 0;
  }

  if (cfg->snapshotStride < 0) {
    setError(err, errSize, "snapshotStride must be nonnegative.");
    return 0;
  }

  setError(err, errSize, "");
  return 1;
}

static void parsePmlType(PmlConfig *pml, const char *value) {
  if (strcmp(value, "none") == 0 || strcmp(value, "off") == 0 ||
      strcmp(value, "0") == 0) {
    pml->type = PML_TYPE_NONE;
    pml->enabled = 0;
  } else if (strcmp(value, "cpml") == 0 || strcmp(value, "pml") == 0 ||
             strcmp(value, "1") == 0) {
    pml->type = PML_TYPE_CPML;
    pml->enabled = 1;
  }
}

int gridLoadStandaloneConfig(SimConfig *cfg, const char *inputDir) {
  FILE *fp;
  char path[512];
  char line[256];
  double *rawAntennaPos = NULL;
  double *rawSourceSamples = NULL;
  int antenna, timeStep;
  int numAntennas;

  /* Preserve the historical standalone CSV behavior: write every step. */
  cfg->snapshotStart = 0;
  cfg->snapshotStride = 1;

  if (inputDir == NULL || inputDir[0] == '\0')
    return 1;

  snprintf(cfg->inputDir, sizeof(cfg->inputDir), "%s", inputDir);
  cfg->useFileInput = 1;
  cfg->pml.type = PML_TYPE_CPML;
  cfg->pml.enabled = 0;

  fdtdJoinPath(path, sizeof(path), inputDir, "grid_config.txt");
  fp = fopen(path, "r");
  if (fp == NULL) {
    fprintf(stderr, "gridLoadStandaloneConfig: could not open %s.\n", path);
    return 0;
  }

  while (fgets(line, sizeof(line), fp) != NULL) {
    char key[64];
    char value[128];
    char *comment = strchr(line, '#');

    if (comment != NULL)
      *comment = '\0';

    if (sscanf(line, "%63s %127s", key, value) != 2)
      continue;

    if (strcmp(key, "Nx") == 0)
      cfg->Nx = atoi(value);
    else if (strcmp(key, "Ny") == 0)
      cfg->Ny = atoi(value);
    else if (strcmp(key, "Nt") == 0)
      cfg->Nt = atoi(value);
    else if (strcmp(key, "dx") == 0)
      cfg->dx = atof(value);
    else if (strcmp(key, "dy") == 0)
      cfg->dy = atof(value);
    else if (strcmp(key, "dt") == 0)
      cfg->dt = atof(value);
    else if (strcmp(key, "num_antennas") == 0)
      cfg->antennas.numAntennas = atoi(value);
    else if (strcmp(key, "ax") == 0)
      cfg->pml.ax = atof(value);
    else if (strcmp(key, "ay") == 0)
      cfg->pml.ay = atof(value);
    else if (strcmp(key, "az") == 0)
      cfg->pml.az = atof(value);
    else if (strcmp(key, "pml_type") == 0 || strcmp(key, "pmlType") == 0)
      parsePmlType(&cfg->pml, value);
    else if (strcmp(key, "pml_enabled") == 0 || strcmp(key, "pmlEnabled") == 0)
      cfg->pml.enabled = atoi(value) != 0;
    else if (strcmp(key, "snapshot_start") == 0 || strcmp(key, "snapshotStart") == 0)
      cfg->snapshotStart = atoi(value);
    else if (strcmp(key, "snapshot_stride") == 0 || strcmp(key, "snapshotStride") == 0)
      cfg->snapshotStride = atoi(value);
  }

  fclose(fp);

  numAntennas = cfg->antennas.numAntennas;
  if (cfg->Nx < 2 || cfg->Ny < 2) {
    fprintf(stderr,
            "gridLoadStandaloneConfig: Nx and Ny must be at least 2.\n");
    return 0;
  }
  if (cfg->Nt < 1) {
    fprintf(stderr, "gridLoadStandaloneConfig: Nt must be positive.\n");
    return 0;
  }
  if (numAntennas < 1) {
    fprintf(stderr,
            "gridLoadStandaloneConfig: num_antennas must be positive.\n");
    return 0;
  }

  rawAntennaPos = (double *)malloc((size_t)numAntennas * 2U * sizeof(double));
  rawSourceSamples = (double *)malloc((size_t)numAntennas *
                                     (size_t)cfg->Nt * sizeof(double));
  cfg->ownedAntennaPos = (double *)malloc((size_t)numAntennas *
                                         2U * sizeof(double));
  cfg->ownedSourceSamples = (double *)malloc((size_t)numAntennas *
                                             (size_t)cfg->Nt * sizeof(double));
  if (rawAntennaPos == NULL || rawSourceSamples == NULL ||
      cfg->ownedAntennaPos == NULL || cfg->ownedSourceSamples == NULL) {
    fprintf(stderr,
            "gridLoadStandaloneConfig: could not allocate source inputs.\n");
    free(rawAntennaPos);
    free(rawSourceSamples);
    simConfigRelease(cfg);
    return 0;
  }

  fdtdJoinPath(path, sizeof(path), inputDir, "antennas.csv");
  if (!fdtdReadMatrixCsv(path, rawAntennaPos, numAntennas, 2)) {
    free(rawAntennaPos);
    free(rawSourceSamples);
    simConfigRelease(cfg);
    return 0;
  }
  for (antenna = 0; antenna < numAntennas; antenna++) {
    cfg->ownedAntennaPos[antenna] = rawAntennaPos[2 * antenna];
    cfg->ownedAntennaPos[numAntennas + antenna] =
        rawAntennaPos[2 * antenna + 1];
  }

  fdtdJoinPath(path, sizeof(path), inputDir, "source.csv");
  if (!fdtdReadMatrixCsv(path, rawSourceSamples, numAntennas, cfg->Nt)) {
    free(rawAntennaPos);
    free(rawSourceSamples);
    simConfigRelease(cfg);
    return 0;
  }
  for (timeStep = 0; timeStep < cfg->Nt; timeStep++)
    for (antenna = 0; antenna < numAntennas; antenna++)
      cfg->ownedSourceSamples[timeStep * numAntennas + antenna] =
          rawSourceSamples[antenna * cfg->Nt + timeStep];

  free(rawAntennaPos);
  free(rawSourceSamples);
  cfg->antennas.pos = cfg->ownedAntennaPos;
  cfg->source.samples = cfg->ownedSourceSamples;
  return 1;
}
