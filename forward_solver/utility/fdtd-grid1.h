#ifndef _FDTD_GRID1_H
#define _FDTD_GRID1_H

enum GRIDTYPE {oneDGrid, teZGrid, tmZGrid, threeDGrid};

typedef enum {
   PML_TYPE_NONE = 0,
   PML_TYPE_CPML = 1
} PmlType;

typedef struct {
   int enabled;
   PmlType type;
   int coefficientsInitialized;
   int diagnosticsPrinted;

   double ax, ay, az;

   double *kx, *ky;
   double *condx, *condy;

   double *bx, *by;
   double *cx, *cy;

   double *psiEzX, *psiEzY;
   double *psiHxY, *psiHyX;
} PmlState;

struct Grid {
   double *hx, *chxh, *chxe;
   double *hy, *chyh, *chye;
   double *hz, *chzh, *chze; // not needed in TM
   double *ex, *cexe, *cexh; // not needed in TM
   double *ey, *ceye, *ceyh; // not needed in TM
   double *ez, *ceze, *cezh;
   int Nx, Ny, sizeZ;
   int time, Nt;
   int type;
   double sourceFreq; // retained for the separate TFSF implementation
   double cdtds; // Courant number
   double dx, dy, dt;
   double eps0, mu0; // standard constants
   double *epsr, *murx, *mury;
   double *cond_e, *cond_m;

   PmlState pml;
};

typedef struct Grid Grid;

#endif
