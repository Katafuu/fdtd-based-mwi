#include "../utility/fdtd-macro-tmz.h"
#include "../utility/fdtd-alloc1.h"
#include "../utility/fdtd-config.h"
#include "../utility/fdtd-io.h"
#include "../utility/fdtd-proto2.h"
#include "pml.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void gridSetConfig(Grid *g, const SimConfig *cfg) {
  double c0;

  Type = tmZGrid;
  NxG(g) = cfg->Nx;
  NyG(g) = cfg->Ny;
  SizeZ = cfg->sizeZ;
  Time = 0;
  NtG(g) = cfg->Nt;
  Dx = cfg->dx;
  Dy = cfg->dy;
  Dt = cfg->dt;
  Cdtds = 0.0;

  Eps0 = 8.854187817e-12;
  Mu0 = 4.0e-7 * 3.14159265358979323846;

  c0 = 1.0 / sqrt(Mu0 * Eps0);
  if (Dt <= 0.0)
    Dt = 0.99 / (c0 * sqrt(1.0 / (Dx * Dx) + 1.0 / (Dy * Dy)));
  SourceFreq = c0 / (20.0 * Dx);
  Cdtds = c0 * Dt / Dx;

}

static void gridAllocFields(Grid *g) {
  ALLOC_2D(g->hx,   NxG(g), NyG(g) - 1, double);
  ALLOC_2D(g->chxh, NxG(g), NyG(g) - 1, double);
  ALLOC_2D(g->chxe, NxG(g), NyG(g) - 1, double);
  ALLOC_2D(g->hy,   NxG(g) - 1, NyG(g), double);
  ALLOC_2D(g->chyh, NxG(g) - 1, NyG(g), double);
  ALLOC_2D(g->chye, NxG(g) - 1, NyG(g), double);
  ALLOC_2D(g->ez,   NxG(g), NyG(g), double);
  ALLOC_2D(g->ceze, NxG(g), NyG(g), double);
  ALLOC_2D(g->cezh, NxG(g), NyG(g), double);

  ALLOC_2D(g->epsr,   NxG(g), NyG(g), double);
  ALLOC_2D(g->murx,   NxG(g), NyG(g) - 1, double);
  ALLOC_2D(g->mury,   NxG(g) - 1, NyG(g), double);
  ALLOC_2D(g->cond_e, NxG(g), NyG(g), double);
  ALLOC_2D(g->cond_m, NxG(g), NyG(g), double);
}

static void gridSetDefaultGeometry(Grid *g) {
  int mm, nn;

  for (mm = 0; mm < NxG(g); mm++)
    for (nn = 0; nn < NyG(g); nn++) {
      Epsr(mm, nn) = 1.0;
      CondE(mm, nn) = 0.0;
      CondM(mm, nn) = 0.0;
    }

  for (mm = 0; mm < NxG(g); mm++)
    for (nn = 0; nn < NyG(g) - 1; nn++)
      Murx(mm, nn) = 1.0;

  for (mm = 0; mm < NxG(g) - 1; mm++)
    for (nn = 0; nn < NyG(g); nn++)
      Mury(mm, nn) = 1.0;
}

static void copyIfProvided(double *dst, const double *src, int count) {
  if (src != NULL)
    memcpy(dst, src, (size_t)count * sizeof(double));
}

static int gridLoadGeometryFiles(Grid *g, const char *inputDir) {
  char path[512];

  fdtdJoinPath(path, sizeof(path), inputDir, "epsr.csv");
  if (!fdtdReadMatrixCsv(path, g->epsr, NxG(g), NyG(g)))
    return 0;

  fdtdJoinPath(path, sizeof(path), inputDir, "murx.csv");
  if (!fdtdReadMatrixCsv(path, g->murx, NxG(g), NyG(g) - 1))
    return 0;

  fdtdJoinPath(path, sizeof(path), inputDir, "mury.csv");
  if (!fdtdReadMatrixCsv(path, g->mury, NxG(g) - 1, NyG(g)))
    return 0;

  fdtdJoinPath(path, sizeof(path), inputDir, "cond_e.csv");
  if (!fdtdReadMatrixCsv(path, g->cond_e, NxG(g), NyG(g)))
    return 0;

  fdtdJoinPath(path, sizeof(path), inputDir, "cond_m.csv");
  return fdtdReadMatrixCsv(path, g->cond_m, NxG(g), NyG(g));
}

static int gridApplyGeometryConfig(Grid *g, const SimConfig *cfg) {
  gridSetDefaultGeometry(g);

  if (cfg->useFileInput && cfg->inputDir[0] != '\0' &&
      !gridLoadGeometryFiles(g, cfg->inputDir))
    return 0;

  copyIfProvided(g->epsr, cfg->epsr, NxG(g) * NyG(g));
  copyIfProvided(g->murx, cfg->murx, NxG(g) * (NyG(g) - 1));
  copyIfProvided(g->mury, cfg->mury, (NxG(g) - 1) * NyG(g));
  copyIfProvided(g->cond_e, cfg->cond_e, NxG(g) * NyG(g));
  copyIfProvided(g->cond_m, cfg->cond_m, NxG(g) * NyG(g));
  return 1;
}

static void gridBuildUpdateCoefficients(Grid *g) {
  double eps, mu, sigma, sigma_m;
  double denom_e, denom_m;
  int mm, nn;

  for (mm = 0; mm < NxG(g); mm++)
    for (nn = 0; nn < NyG(g); nn++) {
      eps = Eps0 * Epsr(mm, nn);
      sigma = CondE(mm, nn);
      denom_e = 1.0 + sigma * Dt / (2.0 * eps);

      Ceze(mm, nn) = (1.0 - sigma * Dt / (2.0 * eps)) / denom_e;
      Cezh(mm, nn) = (Dt / (eps * Dx)) / denom_e;
    }

  for (mm = 0; mm < NxG(g); mm++)
    for (nn = 0; nn < NyG(g) - 1; nn++) {
      mu = Mu0 * Murx(mm, nn);
      sigma_m = CondM(mm, nn);
      denom_m = 1.0 + sigma_m * Dt / (2.0 * mu);

      Chxh(mm, nn) = (1.0 - sigma_m * Dt / (2.0 * mu)) / denom_m;
      Chxe(mm, nn) = (Dt / (mu * Dy)) / denom_m;
    }

  for (mm = 0; mm < NxG(g) - 1; mm++)
    for (nn = 0; nn < NyG(g); nn++) {
      mu = Mu0 * Mury(mm, nn);
      sigma_m = CondM(mm, nn);
      denom_m = 1.0 + sigma_m * Dt / (2.0 * mu);

      Chyh(mm, nn) = (1.0 - sigma_m * Dt / (2.0 * mu)) / denom_m;
      Chye(mm, nn) = (Dt / (mu * Dx)) / denom_m;
    }
}

Grid *gridCreate(const SimConfig *cfg) {
  Grid *g;

  ALLOC_1D(g, 1, Grid);
  if (!gridInitFromConfig(g, cfg)) {
    gridDestroy(g);
    return NULL;
  }

  return g;
}

int gridInitFromConfig(Grid *g, const SimConfig *cfg) {
  char err[256];

  if (!simConfigValidate(cfg, err, sizeof(err))) {
    fprintf(stderr, "gridInitFromConfig: %s\n", err);
    return 0;
  }

  memset(g, 0, sizeof(*g));
  gridSetConfig(g, cfg);


  gridAllocFields(g);
  if (!gridApplyGeometryConfig(g, cfg))
    return 0;
  if (!pmlInit(g, cfg))
    return 0;
  gridBuildUpdateCoefficients(g);

  return 1;
}

void gridDestroy(Grid *g) {
  if (g == NULL)
    return;

  free(g->hx);
  free(g->chxh);
  free(g->chxe);
  free(g->hy);
  free(g->chyh);
  free(g->chye);
  free(g->hz);
  free(g->chzh);
  free(g->chze);
  free(g->ex);
  free(g->cexe);
  free(g->cexh);
  free(g->ey);
  free(g->ceye);
  free(g->ceyh);
  free(g->ez);
  free(g->ceze);
  free(g->cezh);
  free(g->epsr);
  free(g->murx);
  free(g->mury);
  free(g->cond_e);
  free(g->cond_m);
  pmlDestroy(g);
  free(g);
}

void gridInit(Grid *g, const char *inputDir) {
  SimConfig cfg;

  simConfigDefaults(&cfg);
  if (!gridLoadStandaloneConfig(&cfg, inputDir))
    exit(1);
  if (!gridInitFromConfig(g, &cfg)) {
    simConfigRelease(&cfg);
    exit(1);
  }
  simConfigRelease(&cfg);
}

void initTMZCavityMode(Grid *g) {
  const double pi = 3.14159265358979323846;
  const int modeX = 1;
  const int modeY = 1;
  const double e0 = 1.0;
  const double lx = (NxG(g) - 1) * Dx;
  const double ly = (NyG(g) - 1) * Dy;
  const double kx = modeX * pi / lx;
  const double ky = modeY * pi / ly;
  const double eps = Eps0 * Epsr(0, 0);
  const double mu = Mu0 * Murx(0, 0);
  const double omega = sqrt(kx * kx + ky * ky) / sqrt(mu * eps);
  const double hTimeFactor = sin(-0.5 * omega * Dt);
  int mm, nn;

  if (fabs(Murx(0, 0) - Mury(0, 0)) > 1e-12) {
    fprintf(stderr,
      "initTMZCavityMode: murx and mury must match for isotropic cavity mode.\n");
    exit(1);
  }

  for (mm = 0; mm < NxG(g); mm++) {
    const double x = mm * Dx;

    for (nn = 0; nn < NyG(g); nn++) {
      const double y = nn * Dy;
      Ez(mm, nn) = e0 * sin(kx * x) * sin(ky * y);
    }
  }

  for (mm = 0; mm < NxG(g); mm++) {
    const double x = mm * Dx;

    for (nn = 0; nn < NyG(g) - 1; nn++) {
      const double y = (nn + 0.5) * Dy;
      Hx(mm, nn) = -e0 * ky / (mu * omega) *
        sin(kx * x) * cos(ky * y) * hTimeFactor;
    }
  }

  for (mm = 0; mm < NxG(g) - 1; mm++) {
    const double x = (mm + 0.5) * Dx;

    for (nn = 0; nn < NyG(g); nn++) {
      const double y = nn * Dy;
      Hy(mm, nn) = e0 * kx / (mu * omega) *
        cos(kx * x) * sin(ky * y) * hTimeFactor;
    }
  }
}
