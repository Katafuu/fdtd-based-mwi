#ifndef _PML_PRIVATE_H
#define _PML_PRIVATE_H

#include "../utility/fdtd-macro-tmz.h"
#include <math.h>

void pmlInitCoefficients(Grid *g);
void pmlPrintDiagnostics(Grid *g);

#define PML_UPDATE_PSI(PSI, C, B, DERIV) \
  ((PSI) = (C) * (DERIV) + (B) * (PSI))

#define UPDATE_PSI_X(M, N, PSI, DERIV) \
  PML_UPDATE_PSI(PSI, CxPml(M, N), BxPml(M, N), DERIV)
#define UPDATE_PSI_Y(M, N, PSI, DERIV) \
  PML_UPDATE_PSI(PSI, CyPml(M, N), ByPml(M, N), DERIV)

static inline double pml_b(double cond, double kappa, double a,
                           double eps0, double dt) {
  return exp(-((a / eps0) + (cond / (kappa * eps0))) * dt);
}

static inline double pml_c(double cond, double kappa, double a, double b) {
  double denom = cond * kappa + kappa * kappa * a;

  if (denom == 0.0)
    return 0.0;

  return (cond / denom) * (b - 1.0);
}

static inline void upd_psi_ez_x(Grid *g, int mm, int nn) {
  UPDATE_PSI_X(
      mm, nn, PsiEzX(mm, nn),
      (Hy(mm, nn) - Hy(mm - 1, nn)) / Dx);
}

static inline void upd_psi_ez_y(Grid *g, int mm, int nn) {
  UPDATE_PSI_Y(
      mm, nn, PsiEzY(mm, nn),
      (Hx(mm, nn) - Hx(mm, nn - 1)) / Dy);
}

static inline void upd_psi_hx_y(Grid *g, int mm, int nn) {
  UPDATE_PSI_Y(
      mm, nn, PsiHxY(mm, nn),
      (Ez(mm, nn + 1) - Ez(mm, nn)) / Dy);
}

static inline void upd_psi_hy_x(Grid *g, int mm, int nn) {
  UPDATE_PSI_X(
      mm, nn, PsiHyX(mm, nn),
      (Ez(mm + 1, nn) - Ez(mm, nn)) / Dx);
}

#endif
