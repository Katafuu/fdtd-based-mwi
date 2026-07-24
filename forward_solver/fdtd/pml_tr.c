#include "pml.h"
#include "pml-private.h"
#include "../utility/fdtd-macro-tmz.h"

void pmlUpdateHTr(Grid *g) {
  int mm, nn;

  if (!g->pml.enabled)
    return;

  if (!g->pml.coefficientsInitialized)
    pmlInitCoefficients(g);
  pmlPrintDiagnostics(g);

  for (mm = 0; mm < NxG(g); mm++) {
    for (nn = 0; nn < NyG(g) - 1; nn++) {
      upd_psi_hx_y(g, mm, nn);
      Hx(mm, nn) += (Dt / (Mu0 * Murx(mm, nn))) * PsiHxY(mm, nn);
    }
  }

  for (mm = 0; mm < NxG(g) - 1; mm++) {
    for (nn = 0; nn < NyG(g); nn++) {
      upd_psi_hy_x(g, mm, nn);
      Hy(mm, nn) -= (Dt / (Mu0 * Mury(mm, nn))) * PsiHyX(mm, nn);
    }
  }
}

void pmlUpdateETr(Grid *g) {
  int mm, nn;

  if (!g->pml.enabled)
    return;

  if (!g->pml.coefficientsInitialized)
    pmlInitCoefficients(g);

  for (mm = 1; mm < NxG(g) - 1; mm++) {
    for (nn = 1; nn < NyG(g) - 1; nn++) {
      upd_psi_ez_x(g, mm, nn);
      upd_psi_ez_y(g, mm, nn);
      Ez(mm, nn) -= (Dt / (Eps0 * Epsr(mm, nn))) *
        (PsiEzX(mm, nn) - PsiEzY(mm, nn));
    }
  }
}
