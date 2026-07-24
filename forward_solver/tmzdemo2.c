/* TMz simulation with optional file-backed MATLAB geometry input. */

#include <stdio.h>
#include <time.h>
#include "utility/fdtd-config.h"
#include "utility/fdtd-output.h"
#include "utility/fdtd-proto2.h"
#include "fdtd/solver.h"

int main(int argc, char *argv[])
{
  SimConfig cfg;
  Grid *g;
  OutputHandler out;
  char *inputDir = "./results/fdtd_input";
  char *outputFile = "./results/ez_field.csv";
  clock_t start, end;
  double time_spent;

  if (argc > 3) {
    fprintf(stderr, "Usage: %s [input_dir] [output_csv]\n", argv[0]);
    return 1;
  }

  if (argc >= 2)
    inputDir = argv[1];
  if (argc >= 3)
    outputFile = argv[2];

  simConfigDefaults(&cfg);
  if (!gridLoadStandaloneConfig(&cfg, inputDir)) {
    simConfigRelease(&cfg);
    return 1;
  }

  g = gridCreate(&cfg);
  if (g == NULL) {
    simConfigRelease(&cfg);
    return 1;
  }

  if (!fileOutputInit(&out, outputFile)) {
    gridDestroy(g);
    simConfigRelease(&cfg);
    return 1;
  }

  start = clock();
  fdtdRun(g, &cfg, &out);
  end = clock();

  fileOutputDestroy(&out);
  gridDestroy(g);
  simConfigRelease(&cfg);

  time_spent = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Execution time: %f seconds\n", time_spent);

  return 0;
}
