#include "ezinc.h"
#define M_PI 3.14159265358979323846

static double sourceDt = 0.0, sourceDx = 0.0, sourceC0 = 0.0;

/* initialize source-function variables */
void ezIncInit(Grid *g) {
  sourceDt = Dt;
  sourceDx = Dx;
  sourceC0 = 1.0 / sqrt(Mu0 * Eps0);
  return;
}

/* calculate source function at given time and location */
double ezInc(double time, double location, double sourceFreq) {
  double arg;

  if (sourceDt <= 0.0 || sourceDx <= 0.0 || sourceC0 <= 0.0) {
    fprintf(stderr,
        "ezInc: ezIncInit() must be called before ezInc.\n"
        "       Source time/space constants must be positive.\n");
    exit(-1);
  }

  if (sourceFreq <= 0.0) {
    fprintf(stderr,
        "ezInc: source frequency must be positive.\n");
    exit(-1);
  }

  arg = M_PI * (sourceFreq * (time * sourceDt - location * sourceDx / sourceC0) - 1.0);
  arg = arg * arg;

  return (1.0 - 2.0 * arg) * exp(-arg);
}