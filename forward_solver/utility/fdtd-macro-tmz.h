#ifndef _FDTD_MACRO_TMZ_H
#define _FDTD_MACRO_TMZ_H

#include "fdtd-grid1.h"

/* macros that permit the "Grid" to be specified */
/* one-dimensional grid */
#define Hy1G(G, M)     G->hy[M]
#define Chyh1G(G, M)   G->chyh[M]
#define Chye1G(G, M)   G->chye[M]

#define Ez1G(G, M)     G->ez[M]
#define Ceze1G(G, M)   G->ceze[M]
#define Cezh1G(G, M)   G->cezh[M]

/* TMz grid */
#define HxG(G, M, N)   G->hx[(M) * (NyG(G)-1) + (N)]
#define ChxhG(G, M, N) G->chxh[(M) * (NyG(G)-1) + (N)]
#define ChxeG(G, M, N) G->chxe[(M) * (NyG(G)-1) + (N)]

#define HyG(G, M, N)   G->hy[(M) * NyG(G) + (N)]
#define ChyhG(G, M, N) G->chyh[(M) * NyG(G) + (N)]
#define ChyeG(G, M, N) G->chye[(M) * NyG(G) + (N)]

#define EzG(G, M, N)   G->ez[(M) * NyG(G) + (N)]
#define CezeG(G, M, N) G->ceze[(M) * NyG(G) + (N)]
#define CezhG(G, M, N) G->cezh[(M) * NyG(G) + (N)]

#define EpsrG(G, M, N)  G->epsr[(M) * NyG(G) + (N)]
#define MurxG(G, M, N)  G->murx[(M) * (NyG(G)-1) + (N)]
#define MuryG(G, M, N)  G->mury[(M) * NyG(G) + (N)]
#define CondEG(G, M, N) G->cond_e[(M) * NyG(G) + (N)]
#define CondMG(G, M, N) G->cond_m[(M) * NyG(G) + (N)]

#define KxG(G, M, N)    G->pml.kx[(M) * NyG(G) + (N)]
#define KyG(G, M, N)    G->pml.ky[(M) * NyG(G) + (N)]
#define CondxG(G, M, N) G->pml.condx[(M) * NyG(G) + (N)]
#define CondyG(G, M, N) G->pml.condy[(M) * NyG(G) + (N)]

#define BxPmlG(G, M, N) G->pml.bx[(M) * NyG(G) + (N)]
#define ByPmlG(G, M, N) G->pml.by[(M) * NyG(G) + (N)]
#define CxPmlG(G, M, N) G->pml.cx[(M) * NyG(G) + (N)]
#define CyPmlG(G, M, N) G->pml.cy[(M) * NyG(G) + (N)]

#define PsiEzXG(G, M, N) G->pml.psiEzX[(M) * NyG(G) + (N)]
#define PsiEzYG(G, M, N) G->pml.psiEzY[(M) * NyG(G) + (N)]
#define PsiHxYG(G, M, N) G->pml.psiHxY[(M) * (NyG(G)-1) + (N)]
#define PsiHyXG(G, M, N) G->pml.psiHyX[(M) * NyG(G) + (N)]

#define NxG(G)      G->Nx
#define NyG(G)      G->Ny
#define SizeZG(G)      G->sizeZ
#define TimeG(G)       G->time
#define NtG(G)    G->Nt
#define SourceFreqG(G) G->sourceFreq
#define CdtdsG(G)      G->cdtds
#define DxG(G)         G->dx
#define DyG(G)         G->dy
#define DtG(G)         G->dt
#define TypeG(G)       G->type
#define AxG(G)         G->pml.ax
#define AyG(G)         G->pml.ay
#define AzG(G)         G->pml.az
#define Eps0G(G)       G->eps0
#define Mu0G(G)        G->mu0

/* macros that assume the "Grid" is "g" */
/* one-dimensional grid */
#define Hy1(M)         Hy1G(g, M)
#define Chyh1(M)       Chyh1G(g, M)
#define Chye1(M)       Chye1G(g, M)

#define Ez1(M)         Ez1G(g, M)
#define Ceze1(M)       Ceze1G(g, M)
#define Cezh1(M)       Cezh1G(g, M)

/* TMz grid */
#define Hx(M, N)       HxG(g, M, N)
#define Chxh(M, N)     ChxhG(g, M, N)
#define Chxe(M, N)     ChxeG(g, M, N)

#define Hy(M, N)       HyG(g, M, N)
#define Chyh(M, N)     ChyhG(g, M, N)
#define Chye(M, N)     ChyeG(g, M, N)

#define Ez(M, N)       EzG(g, M, N)
#define Ceze(M, N)     CezeG(g, M, N)
#define Cezh(M, N)     CezhG(g, M, N)

#define Epsr(M, N)     EpsrG(g, M, N)
#define Murx(M, N)     MurxG(g, M, N)
#define Mury(M, N)     MuryG(g, M, N)
#define CondE(M, N)    CondEG(g, M, N)
#define CondM(M, N)    CondMG(g, M, N)

#define Kx(M, N)       KxG(g, M, N)
#define Ky(M, N)       KyG(g, M, N)
#define Condx(M, N)    CondxG(g, M, N)
#define Condy(M, N)    CondyG(g, M, N)

#define BxPml(M, N)    BxPmlG(g, M, N)
#define ByPml(M, N)    ByPmlG(g, M, N)
#define CxPml(M, N)    CxPmlG(g, M, N)
#define CyPml(M, N)    CyPmlG(g, M, N)
#define PsiEzX(M, N)   PsiEzXG(g, M, N)
#define PsiEzY(M, N)   PsiEzYG(g, M, N)
#define PsiHxY(M, N)   PsiHxYG(g, M, N)
#define PsiHyX(M, N)   PsiHyXG(g, M, N)

#define SizeZ           SizeZG(g)
#define Time            TimeG(g)
#define SourceFreq      SourceFreqG(g)
#define Type            TypeG(g)

#define Cdtds           CdtdsG(g)
#define Dx              DxG(g)
#define Dy              DyG(g)
#define Dt              DtG(g)

#define Ax              AxG(g)
#define Ay              AyG(g)
#define Az              AzG(g)

#define Eps0            Eps0G(g)
#define Mu0             Mu0G(g)

#endif      /* matches #ifndef _FDTD_MACRO_TMZ_H */
