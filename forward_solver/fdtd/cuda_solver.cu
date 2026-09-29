/* Optional double-precision CUDA backend. CPU equations and update order are
 * retained. No MATLAB API calls occur in kernels or device-buffer helpers. */
#include "cuda_solver.h"
#include <cuda_runtime.h>
#include <vector>
#include <stdexcept>
#include <cstdio>

static void checked(cudaError_t status) {
  if (status != cudaSuccess) throw std::runtime_error(cudaGetErrorString(status));
}
struct Buffers {
  std::vector<void *> owned;
  ~Buffers() { for (void *p : owned) cudaFree(p); }
  double *allocate(size_t n, const double *source = nullptr) {
    if (!n) return nullptr;
    double *p = nullptr; checked(cudaMalloc((void **)&p,n*sizeof(double)));
    owned.push_back(p);
    if (source) checked(cudaMemcpy(p,source,n*sizeof(double),cudaMemcpyHostToDevice));
    else checked(cudaMemset(p,0,n*sizeof(double)));
    return p;
  }
};

__global__ static void magnetic(Grid g, int reverse, int pml) {
  int i = blockIdx.x*blockDim.x+threadIdx.x;
  if (i >= g.Nx*g.Ny) return;
  int x = i/g.Ny, y = i%g.Ny;
  if (y < g.Ny-1) {
    int h = x*(g.Ny-1)+y;
    if (!pml) {
      double term = g.chxe[h]*(g.ez[i+1]-g.ez[i])*1.0/g.pml.ky[i];
      g.hx[h] = reverse ? g.chxh[h]*g.hx[h]+term : g.chxh[h]*g.hx[h]-term;
    } else {
      g.pml.psiHxY[h] = g.pml.cy[i]*((g.ez[i+1]-g.ez[i])/g.dy)+g.pml.by[i]*g.pml.psiHxY[h];
      double term = (g.dt/(g.mu0*g.murx[h]))*g.pml.psiHxY[h];
      if (reverse) g.hx[h] += term; else g.hx[h] -= term;
    }
  }
  if (x < g.Nx-1) {
    if (!pml) {
      double term = g.chye[i]*(g.ez[i+g.Ny]-g.ez[i])*1.0/g.pml.kx[i];
      g.hy[i] = reverse ? g.chyh[i]*g.hy[i]-term : g.chyh[i]*g.hy[i]+term;
    } else {
      g.pml.psiHyX[i] = g.pml.cx[i]*((g.ez[i+g.Ny]-g.ez[i])/g.dx)+g.pml.bx[i]*g.pml.psiHyX[i];
      double term = (g.dt/(g.mu0*g.mury[i]))*g.pml.psiHyX[i];
      if (reverse) g.hy[i] -= term; else g.hy[i] += term;
    }
  }
}
__global__ static void electric(Grid g, int reverse, int pml) {
  int i = blockIdx.x*blockDim.x+threadIdx.x;
  if (i >= g.Nx*g.Ny) return;
  int x = i/g.Ny, y = i%g.Ny, h = x*(g.Ny-1)+y;
  if (x == 0 || y == 0 || x == g.Nx-1 || y == g.Ny-1) return;
  if (!pml) {
    double term = g.cezh[i]*((g.hy[i]-g.hy[i-g.Ny])*1.0/g.pml.kx[i]-
                             (g.hx[h]-g.hx[h-1])*1.0/g.pml.ky[i]);
    g.ez[i] = reverse ? g.ceze[i]*g.ez[i]-term : g.ceze[i]*g.ez[i]+term;
  } else {
    g.pml.psiEzX[i] = g.pml.cx[i]*((g.hy[i]-g.hy[i-g.Ny])/g.dx)+g.pml.bx[i]*g.pml.psiEzX[i];
    g.pml.psiEzY[i] = g.pml.cy[i]*((g.hx[h]-g.hx[h-1])/g.dy)+g.pml.by[i]*g.pml.psiEzY[i];
    double term = (g.dt/(g.eps0*g.epsr[i]))*(g.pml.psiEzX[i]-g.pml.psiEzY[i]);
    if (reverse) g.ez[i] -= term; else g.ez[i] += term;
  }
}
__global__ static void sources(Grid g, const double *positions, const double *samples,
                               int antennas, int time) {
  // Preserve injection order, including coincident source coordinates.
  for (int a=0; a<antennas; ++a) {
    int x=(int)positions[a]-1, y=(int)positions[antennas+a]-1;
    g.ez[x*g.Ny+y] += samples[(size_t)time*antennas+a];
  }
}
__global__ static void record(Grid g, CudaOutputs out, const double *positions,
                              int antennas, int time, int frame) {
  int i = blockIdx.x*blockDim.x+threadIdx.x;
  if (out.rx && i < antennas) {
    int x=(int)positions[i]-1, y=(int)positions[antennas+i]-1;
    out.rx[(size_t)time*antennas+i]=g.ez[x*g.Ny+y];
  }
  if (frame < 0) return;
  if (out.ez && i < out.nx*out.ny) {
    int x=i%out.nx, y=i/out.nx;
    out.ez[(size_t)frame*out.nx*out.ny+i]=g.ez[(x+out.x)*g.Ny+y+out.y];
  }
  if (out.hx && i < g.Nx*(g.Ny-1)) {
    int x=i%g.Nx, y=i/g.Nx;
    out.hx[(size_t)frame*g.Nx*(g.Ny-1)+i]=g.hx[x*(g.Ny-1)+y];
  }
  if (out.hy && i < (g.Nx-1)*g.Ny) {
    int x=i%(g.Nx-1), y=i/(g.Nx-1);
    out.hy[(size_t)frame*(g.Nx-1)*g.Ny+i]=g.hy[x*g.Ny+y];
  }
}

extern "C" int fdtdRunCuda(const Grid *host, const SimConfig *cfg, int reverse,
                            const CudaOutputs *output, char *error, size_t errorSize) {
  try {
    Grid g=*host; Buffers buffers;
    size_t n=(size_t)g.Nx*g.Ny, hx=(size_t)g.Nx*(g.Ny-1), hy=(size_t)(g.Nx-1)*g.Ny;
    size_t ezCount=output->ez ? (size_t)output->nx*output->ny*output->frames : 0;
    size_t hxCount=output->hx ? hx*output->frames : 0;
    size_t hyCount=output->hy ? hy*output->frames : 0;
    size_t rxCount=output->rx ? (size_t)cfg->antennas.numAntennas*cfg->Nt : 0;
    size_t freeBytes,totalBytes; checked(cudaMemGetInfo(&freeBytes,&totalBytes));
    size_t fieldCount=6*n+4*hx+4*hy+(g.pml.enabled ? 8*n+hx+hy : 0);
    size_t inputs=(size_t)cfg->antennas.numAntennas*(cfg->Nt+2);
    if (8.0*(fieldCount+inputs+ezCount+hxCount+hyCount+rxCount) > .8*freeBytes)
      throw std::runtime_error("Insufficient free VRAM for double-precision fields and requested histories (20% reserve).");
#define COPY(member,count) g.member=buffers.allocate(count,host->member)
    COPY(ez,n); COPY(ceze,n); COPY(cezh,n); COPY(epsr,n);
    COPY(hx,hx); COPY(chxh,hx); COPY(chxe,hx); COPY(murx,hx);
    COPY(hy,hy); COPY(chyh,hy); COPY(chye,hy); COPY(mury,hy);
    COPY(pml.kx,n); COPY(pml.ky,n);
    if (g.pml.enabled) {
      COPY(pml.bx,n); COPY(pml.by,n); COPY(pml.cx,n); COPY(pml.cy,n);
      COPY(pml.psiEzX,n); COPY(pml.psiEzY,n); COPY(pml.psiHxY,hx); COPY(pml.psiHyX,hy);
    }
#undef COPY
    double *positions=buffers.allocate(2*cfg->antennas.numAntennas,cfg->antennas.pos);
    double *samples=buffers.allocate((size_t)cfg->Nt*cfg->antennas.numAntennas,cfg->source.samples);
    CudaOutputs device=*output;
    device.ez=buffers.allocate(ezCount); device.hx=buffers.allocate(hxCount);
    device.hy=buffers.allocate(hyCount); device.rx=buffers.allocate(rxCount);
    size_t lanes=n > (size_t)cfg->antennas.numAntennas ? n : (size_t)cfg->antennas.numAntennas;
    int blocks=(int)((lanes+255)/256), frame=0;
    for (int t=0;t<cfg->Nt;++t) {
      magnetic<<<blocks,256>>>(g,reverse,0);
      if (g.pml.enabled) magnetic<<<blocks,256>>>(g,reverse,1);
      electric<<<blocks,256>>>(g,reverse,0);
      if (g.pml.enabled) electric<<<blocks,256>>>(g,reverse,1);
      sources<<<1,1>>>(g,positions,samples,cfg->antennas.numAntennas,t);
      int first=cfg->snapshotStart<0 ? 0 : cfg->snapshotStart;
      bool snapshot=cfg->snapshotStride<=0 ? t==cfg->Nt-1 : t>=first && (t-first)%cfg->snapshotStride==0;
      int f=snapshot && frame<output->frames ? frame++ : -1;
      record<<<blocks,256>>>(g,device,positions,cfg->antennas.numAntennas,t,f);
      checked(cudaGetLastError());
    }
    checked(cudaDeviceSynchronize());
#define DOWNLOAD(member,count) if(count) checked(cudaMemcpy(output->member,device.member,(count)*sizeof(double),cudaMemcpyDeviceToHost))
    DOWNLOAD(ez,ezCount); DOWNLOAD(hx,hxCount); DOWNLOAD(hy,hyCount); DOWNLOAD(rx,rxCount);
#undef DOWNLOAD
    return 1;
  } catch (const std::exception &exception) {
    std::snprintf(error,errorSize,"%s",exception.what()); return 0;
  }
}
