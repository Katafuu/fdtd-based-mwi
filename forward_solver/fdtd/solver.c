#include "solver.h"
#include "pml.h"
#include "../utility/fdtd-macro-tmz.h"
#include "../utility/fdtd-proto2.h"

void fdtdRun(Grid *g, const SimConfig *cfg, OutputHandler *out) {
  int antenna;

  if (out != NULL && out->begin != NULL)
    out->begin(out, g, cfg);


  for (Time = 0; Time < NtG(g); Time++) {
    updateH2d(g);
    pmlUpdateH(g);

    updateE2d(g);
    pmlUpdateE(g);

    for (antenna = 0; antenna < cfg->antennas.numAntennas; antenna++) {
      const int sourceX = (int)cfg->antennas.pos[antenna] - 1;
      const int sourceY =
          (int)cfg->antennas.pos[cfg->antennas.numAntennas + antenna] - 1;
      const size_t sampleIndex =
          (size_t)Time * (size_t)cfg->antennas.numAntennas +
          (size_t)antenna;

      EzG(g, sourceX, sourceY) += cfg->source.samples[sampleIndex];
    }

    if (out != NULL && out->record != NULL)
      out->record(out, g, Time);
  }

  if (out != NULL && out->end != NULL)
    out->end(out, g);
}
