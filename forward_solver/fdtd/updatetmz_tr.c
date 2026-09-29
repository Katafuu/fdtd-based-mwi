#include "../utility/fdtd-macro-tmz.h"

/* update magnetic field using the time-reversed curl signs */
void updateH2dTr(Grid *g) {
  int mm, nn;

  if (Type == oneDGrid) {
    #ifdef _OPENMP
    #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)
    #endif
    for (mm = 0; mm < NxG(g) - 1; mm++) {
      Hy1(mm) = Chyh1(mm) * Hy1(mm)
        - Chye1(mm) * (Ez1(mm + 1) - Ez1(mm));
    }
  } else {
    #ifdef _OPENMP
    #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)
    #endif
    for (mm = 0; mm < NxG(g); mm++) {
      for (nn = 0; nn < NyG(g) - 1; nn++) {
        Hx(mm, nn) = Chxh(mm, nn) * Hx(mm, nn)
          + Chxe(mm, nn) * (Ez(mm, nn + 1) - Ez(mm, nn))
            * 1.0 / Ky(mm, nn);
      }
    }

    #ifdef _OPENMP

    #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)

    #endif

    for (mm = 0; mm < NxG(g) - 1; mm++) {
      for (nn = 0; nn < NyG(g); nn++) {
        Hy(mm, nn) = Chyh(mm, nn) * Hy(mm, nn)
          - Chye(mm, nn) * (Ez(mm + 1, nn) - Ez(mm, nn))
            * 1.0 / Kx(mm, nn);
      }
    }
  }
}

/* update electric field using the time-reversed curl sign */
void updateE2dTr(Grid *g) {
  int mm, nn;

  if (Type == oneDGrid) {
    #ifdef _OPENMP
    #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)
    #endif
    for (mm = 1; mm < NxG(g) - 1; mm++) {
      Ez1(mm) = Ceze1(mm) * Ez1(mm)
        - Cezh1(mm) * (Hy1(mm) - Hy1(mm - 1));
    }
  } else {
    #ifdef _OPENMP
    #pragma omp parallel for private(nn) schedule(static) num_threads(g->solverThreads) if(g->solverThreads > 1)
    #endif
    for (mm = 1; mm < NxG(g) - 1; mm++) {
      for (nn = 1; nn < NyG(g) - 1; nn++) {
        Ez(mm, nn) = Ceze(mm, nn) * Ez(mm, nn) -
          Cezh(mm, nn) *
            ((Hy(mm, nn) - Hy(mm - 1, nn)) * 1.0 / Kx(mm, nn) -
             (Hx(mm, nn) - Hx(mm, nn - 1)) * 1.0 / Ky(mm, nn));
      }
    }
  }
}
