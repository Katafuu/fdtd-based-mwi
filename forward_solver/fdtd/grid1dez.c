#include <math.h>
#include "../utility/fdtd-macro-tmz.h"
#include "../utility/fdtd-alloc1.h"

#define NLOSS            20   // number of lossy cells at end of 1D grid
#define MAX_LOSS         0.35 // maximum loss factor in lossy layer

void gridInit1d(Grid *g) {
  double imp0 = 377.0, depthInLayer, lossFactor;
  int mm;

  NxG(g) += NLOSS;            // size of domain
  Type = oneDGrid;           // set grid type

  ALLOC_1D(g->hy,   NxG(g) - 1, double);
  ALLOC_1D(g->chyh, NxG(g) - 1, double);
  ALLOC_1D(g->chye, NxG(g) - 1, double);
  ALLOC_1D(g->ez,   NxG(g), double);
  ALLOC_1D(g->ceze, NxG(g), double);
  ALLOC_1D(g->cezh, NxG(g), double);

  /* set the electric- and magnetic-field update coefficients */
  for (mm = 0; mm < NxG(g) - 1; mm++) {
    if (mm < NxG(g) - 1 - NLOSS) {
      Ceze1(mm) = 1.0;
      Cezh1(mm) = Cdtds * imp0;
      Chyh1(mm) = 1.0;
      Chye1(mm) = Cdtds / imp0;
    } else {
      depthInLayer = mm - (NxG(g) - 1 - NLOSS) + 0.5;
      lossFactor = MAX_LOSS * pow(depthInLayer / NLOSS, 2);
      Ceze1(mm) = (1.0 - lossFactor) / (1.0 + lossFactor);
      Cezh1(mm) = Cdtds * imp0 / (1.0 + lossFactor);
      depthInLayer += 0.5;
            lossFactor = MAX_LOSS * pow(depthInLayer / NLOSS, 2);
            Chyh1(mm) = (1.0 - lossFactor) / (1.0 + lossFactor);
            Chye1(mm) = Cdtds / imp0 / (1.0 + lossFactor);
        }
    }

    return;
}
