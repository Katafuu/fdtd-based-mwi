#include "pml.h"
#include "pml-private.h"
#include "../utility/fdtd-alloc1.h"
#include "../utility/fdtd-io.h"
#include "../utility/fdtd-macro-tmz.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int pmlTryReadMatrixCsv(const char *filename, double *data,
                               int rows, int cols) {
  FILE *fp = fopen(filename, "r");

  if (fp == NULL)
    return 0;

  fclose(fp);
  if (!fdtdReadMatrixCsv(filename, data, rows, cols))
    return -1;
  return 1;
}

static void pmlAllocProfiles(Grid *g) {
  ALLOC_2D(g->pml.kx,    NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.ky,    NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.condx, NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.condy, NxG(g), NyG(g), double);
}

static void pmlAllocRuntime(Grid *g) {
  ALLOC_2D(g->pml.bx, NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.by, NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.cx, NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.cy, NxG(g), NyG(g), double);

  ALLOC_2D(g->pml.psiEzX, NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.psiEzY, NxG(g), NyG(g), double);
  ALLOC_2D(g->pml.psiHxY, NxG(g), NyG(g) - 1, double);
  ALLOC_2D(g->pml.psiHyX, NxG(g) - 1, NyG(g), double);
}

static void pmlSetProfileDefaults(Grid *g) {
  int mm, nn;

  for (mm = 0; mm < NxG(g); mm++)
    for (nn = 0; nn < NyG(g); nn++) {
      Kx(mm, nn) = 1.0;
      Ky(mm, nn) = 1.0;
      Condx(mm, nn) = 0.0;
      Condy(mm, nn) = 0.0;
    }
}

static void pmlCopyIfProvided(double *dst, const double *src, int count,
                              int *enabled) {
  if (src != NULL) {
    memcpy(dst, src, (size_t)count * sizeof(double));
    *enabled = 1;
  }
}

static int pmlLoadProfileFiles(Grid *g, const SimConfig *cfg, int *enabled) {
  char path[512];
  int status;

  if (!cfg->useFileInput || cfg->inputDir[0] == '\0')
    return 1;

  fdtdJoinPath(path, sizeof(path), cfg->inputDir, "condx.csv");
  status = pmlTryReadMatrixCsv(path, g->pml.condx, NxG(g), NyG(g));
  if (status < 0)
    return 0;
  if (status > 0)
    *enabled = 1;

  fdtdJoinPath(path, sizeof(path), cfg->inputDir, "condy.csv");
  status = pmlTryReadMatrixCsv(path, g->pml.condy, NxG(g), NyG(g));
  if (status < 0)
    return 0;
  if (status > 0)
    *enabled = 1;

  fdtdJoinPath(path, sizeof(path), cfg->inputDir, "kx.csv");
  if (pmlTryReadMatrixCsv(path, g->pml.kx, NxG(g), NyG(g)) < 0)
    return 0;

  fdtdJoinPath(path, sizeof(path), cfg->inputDir, "ky.csv");
  if (pmlTryReadMatrixCsv(path, g->pml.ky, NxG(g), NyG(g)) < 0)
    return 0;

  return 1;
}

void pmlInitCoefficients(Grid *g) {
  int mm, nn;

  for (mm = 0; mm < NxG(g); mm++)
    for (nn = 0; nn < NyG(g); nn++) {
      BxPml(mm, nn) = pml_b(Condx(mm, nn), Kx(mm, nn), Ax, Eps0, Dt);
      ByPml(mm, nn) = pml_b(Condy(mm, nn), Ky(mm, nn), Ay, Eps0, Dt);

      CxPml(mm, nn) = pml_c(Condx(mm, nn), Kx(mm, nn), Ax, BxPml(mm, nn));
      CyPml(mm, nn) = pml_c(Condy(mm, nn), Ky(mm, nn), Ay, ByPml(mm, nn));
    }

  g->pml.coefficientsInitialized = 1;
}

void pmlPrintDiagnostics(Grid *g) {
  double condxMax = 0.0, condyMax = 0.0;
  double kxMin = 1.0, kxMax = 0.0, kyMin = 1.0, kyMax = 0.0;
  double cxMin = 0.0, cxMax = 0.0, cyMin = 0.0, cyMax = 0.0;
  double bxMin = 1.0, bxMax = 0.0, byMin = 1.0, byMax = 0.0;
  int condxNonzero = 0, condyNonzero = 0;
  int mm, nn;

  if (g->pml.diagnosticsPrinted || !g->pml.enabled)
    return;

  for (mm = 0; mm < NxG(g); mm++) {
    for (nn = 0; nn < NyG(g); nn++) {
      if (Condx(mm, nn) != 0.0)
        condxNonzero++;
      if (Condy(mm, nn) != 0.0)
        condyNonzero++;

      if (Condx(mm, nn) > condxMax)
        condxMax = Condx(mm, nn);
      if (Condy(mm, nn) > condyMax)
        condyMax = Condy(mm, nn);

      if (Kx(mm, nn) < kxMin)
        kxMin = Kx(mm, nn);
      if (Kx(mm, nn) > kxMax)
        kxMax = Kx(mm, nn);
      if (Ky(mm, nn) < kyMin)
        kyMin = Ky(mm, nn);
      if (Ky(mm, nn) > kyMax)
        kyMax = Ky(mm, nn);

      if (CxPml(mm, nn) < cxMin)
        cxMin = CxPml(mm, nn);
      if (CxPml(mm, nn) > cxMax)
        cxMax = CxPml(mm, nn);
      if (CyPml(mm, nn) < cyMin)
        cyMin = CyPml(mm, nn);
      if (CyPml(mm, nn) > cyMax)
        cyMax = CyPml(mm, nn);

      if (BxPml(mm, nn) < bxMin)
        bxMin = BxPml(mm, nn);
      if (BxPml(mm, nn) > bxMax)
        bxMax = BxPml(mm, nn);
      if (ByPml(mm, nn) < byMin)
        byMin = ByPml(mm, nn);
      if (ByPml(mm, nn) > byMax)
        byMax = ByPml(mm, nn);
    }
  }

  printf("PML diagnostics:\n");
  printf("  grid: Nx=%d Ny=%d Dx=%.17g Dy=%.17g Dt=%.17g\n",
         NxG(g), NyG(g), Dx, Dy, Dt);
  printf("  Condx: nonzero=%d max=%.17g\n", condxNonzero, condxMax);
  printf("  Condy: nonzero=%d max=%.17g\n", condyNonzero, condyMax);
  printf("  Kx: min=%.17g max=%.17g\n", kxMin, kxMax);
  printf("  Ky: min=%.17g max=%.17g\n", kyMin, kyMax);
  printf("  CxPml: min=%.17g max=%.17g\n", cxMin, cxMax);
  printf("  CyPml: min=%.17g max=%.17g\n", cyMin, cyMax);
  printf("  BxPml: min=%.17g max=%.17g\n", bxMin, bxMax);
  printf("  ByPml: min=%.17g max=%.17g\n", byMin, byMax);

  g->pml.diagnosticsPrinted = 1;
}

int pmlInit(Grid *g, const SimConfig *cfg) {
  int enabled = cfg->pml.enabled;

  g->pml.type = cfg->pml.type;
  g->pml.enabled = cfg->pml.enabled;
  g->pml.coefficientsInitialized = 0;
  g->pml.diagnosticsPrinted = 0;
  g->pml.ax = cfg->pml.ax;
  g->pml.ay = cfg->pml.ay;
  g->pml.az = cfg->pml.az;

  pmlAllocProfiles(g);
  pmlSetProfileDefaults(g);

  if (!pmlLoadProfileFiles(g, cfg, &enabled))
    return 0;
  pmlCopyIfProvided(g->pml.condx, cfg->pml.condx, NxG(g) * NyG(g), &enabled);
  pmlCopyIfProvided(g->pml.condy, cfg->pml.condy, NxG(g) * NyG(g), &enabled);
  pmlCopyIfProvided(g->pml.kx, cfg->pml.kx, NxG(g) * NyG(g), &enabled);
  pmlCopyIfProvided(g->pml.ky, cfg->pml.ky, NxG(g) * NyG(g), &enabled);

  if (enabled && g->pml.type == PML_TYPE_NONE)
    g->pml.type = PML_TYPE_CPML;
  g->pml.enabled = enabled && g->pml.type != PML_TYPE_NONE;

  if (g->pml.enabled) {
    pmlAllocRuntime(g);
    pmlInitCoefficients(g);
  }

  return 1;
}

void pmlDestroy(Grid *g) {
  free(g->pml.kx);
  free(g->pml.ky);
  free(g->pml.condx);
  free(g->pml.condy);
  free(g->pml.bx);
  free(g->pml.by);
  free(g->pml.cx);
  free(g->pml.cy);
  free(g->pml.psiEzX);
  free(g->pml.psiEzY);
  free(g->pml.psiHxY);
  free(g->pml.psiHyX);
  memset(&g->pml, 0, sizeof(g->pml));
}

void pmlUpdateH(Grid *g) {
  int mm, nn;

  if (!g->pml.enabled)
    return;

  if (!g->pml.coefficientsInitialized)
    pmlInitCoefficients(g);
  pmlPrintDiagnostics(g);

  #ifdef _OPENMP

  #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)

  #endif

  for (mm = 0; mm < NxG(g); mm++) {
    for (nn = 0; nn < NyG(g) - 1; nn++) {
      upd_psi_hx_y(g, mm, nn);
      Hx(mm, nn) -= (Dt / (Mu0 * Murx(mm, nn))) * PsiHxY(mm, nn);
    }
  }

  #ifdef _OPENMP

  #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)

  #endif

  for (mm = 0; mm < NxG(g) - 1; mm++) {
    for (nn = 0; nn < NyG(g); nn++) {
      upd_psi_hy_x(g, mm, nn);
      Hy(mm, nn) += (Dt / (Mu0 * Mury(mm, nn))) * PsiHyX(mm, nn);
    }
  }
}

void pmlUpdateE(Grid *g) {
  int mm, nn;

  if (!g->pml.enabled)
    return;

  if (!g->pml.coefficientsInitialized)
    pmlInitCoefficients(g);

  #ifdef _OPENMP

  #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)

  #endif

  for (mm = 1; mm < NxG(g) - 1; mm++) {
    for (nn = 1; nn < NyG(g) - 1; nn++) {
      upd_psi_ez_x(g, mm, nn);
      upd_psi_ez_y(g, mm, nn);
      Ez(mm, nn) += (Dt / (Eps0 * Epsr(mm, nn))) *
        (PsiEzX(mm, nn) - PsiEzY(mm, nn));
    }
  }
}

void pml_update_magnetic(Grid *g) {
  pmlUpdateH(g);
}

void pml_update_electric(Grid *g) {
  pmlUpdateE(g);
}
