/* Second-order ABC for TMz grid. */
#include <math.h>
#include "../utility/fdtd-alloc1.h"
#include "../utility/fdtd-macro-tmz.h"

/* Define macros for arrays that store the previous values of the
 * fields. For each one of these arrays the three arguments are as
 * follows:
 *
 *   first argument: spatial displacement from the boundary
 *   second argument: displacement back in time
 *   third argument: distance from either the bottom (if EzLeft or
 *                     EzRight) or left (if EzTop or EzBottom) side
 *                     of grid
 *
 */
#define EzLeft(M, Q, N)    ezLeft[(N) * 6 + (Q) * 3 + (M)]
#define EzRight(M, Q, N) ezRight[(N) * 6 + (Q) * 3 + (M)]
#define EzTop(N, Q, M)        ezTop[(M) * 6 + (Q) * 3 + (N)]
#define EzBottom(N, Q, M) ezBottom[(M) * 6 + (Q) * 3 + (N)]

static int initDone = 0;
static double coef0, coef1, coef2;
static double *ezLeft, *ezRight, *ezTop, *ezBottom;

void abcInit(Grid *g) {
  double temp1, temp2;

    initDone = 1;

    /* allocate memory for ABC arrays */
    ALLOC_1D(ezLeft, NyG(g) * 6, double);
    ALLOC_1D(ezRight, NyG(g) * 6, double);
    ALLOC_1D(ezTop, NxG(g) * 6, double);
    ALLOC_1D(ezBottom, NxG(g) * 6, double);

    /* calculate ABC coefficients */
    temp1 = sqrt(Cezh(0, 0) * Chye(0, 0));
    temp2 = 1.0 / temp1 + 2.0 + temp1;
    coef0 = -(1.0 / temp1 - 2.0 + temp1) / temp2;
    coef1 = -2.0 * (temp1 - 1.0 / temp1) / temp2;
    coef2 = 4.0 * (temp1 + 1.0 / temp1) / temp2;

    return;
}

void abc(Grid *g)
{
    int mm, nn;

    /* ABC at left side of grid */
    for (nn = 0; nn < NyG(g); nn++) {
      Ez(0, nn) = coef0 * (Ez(2, nn) + EzLeft(0, 1, nn))
       + coef1 * (EzLeft(0, 0, nn) + EzLeft(2, 0, nn)
                  - Ez(1, nn) - EzLeft(1, 1, nn))
       + coef2 * EzLeft(1, 0, nn) - EzLeft(2, 1, nn);

        /* memorize old fields */
        for (mm = 0; mm < 3; mm++) {
          EzLeft(mm, 1, nn) = EzLeft(mm, 0, nn);
          EzLeft(mm, 0, nn) = Ez(mm, nn);
        }
    }

    /* ABC at right side of grid */
    for (nn = 0; nn < NyG(g); nn++) {
      Ez(NxG(g) - 1, nn) = coef0 * (Ez(NxG(g) - 3, nn) + EzRight(0, 1, nn))
       + coef1 * (EzRight(0, 0, nn) + EzRight(2, 0, nn)
                  - Ez(NxG(g) - 2, nn) - EzRight(1, 1, nn))
       + coef2 * EzRight(1, 0, nn) - EzRight(2, 1, nn);

        /* memorize old fields */
        for (mm = 0; mm < 3; mm++) {
          EzRight(mm, 1, nn) = EzRight(mm, 0, nn);
          EzRight(mm, 0, nn) = Ez(NxG(g) - 1 - mm, nn);
        }
    }

    /* ABC at bottom of grid */
    for (mm = 0; mm < NxG(g); mm++) {
      Ez(mm, 0) = coef0 * (Ez(mm, 2) + EzBottom(0, 1, mm))
       + coef1 * (EzBottom(0, 0, mm) + EzBottom(2, 0, mm)
                  - Ez(mm, 1) - EzBottom(1, 1, mm))
       + coef2 * EzBottom(1, 0, mm) - EzBottom(2, 1, mm);

        /* memorize old fields */
        for (nn = 0; nn < 3; nn++) {
          EzBottom(nn, 1, mm) = EzBottom(nn, 0, mm);
          EzBottom(nn, 0, mm) = Ez(mm, nn);
        }
    }

    /* ABC at top of grid */
    for (mm = 0; mm < NxG(g); mm++) {
        Ez(mm, NyG(g) - 1) = coef0 * (Ez(mm, NyG(g) - 3) + EzTop(0, 1, mm))
         + coef1 * (EzTop(0, 0, mm) + EzTop(2, 0, mm)
                    - Ez(mm, NyG(g) - 2) - EzTop(1, 1, mm))
         + coef2 * EzTop(1, 0, mm) - EzTop(2, 1, mm);

        /* memorize old fields */
        for (nn = 0; nn < 3; nn++) {
          EzTop(nn, 1, mm) = EzTop(nn, 0, mm);
          EzTop(nn, 0, mm) = Ez(mm, NyG(g) - 1 - nn);
        }
    }

    return;
}
