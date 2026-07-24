#ifndef _EZINC_H
#define _EZINC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../utility/fdtd-macro-tmz.h"

void ezIncInit(Grid *g);
double ezInc(double time, double location, double sourceFreq);

#endif
