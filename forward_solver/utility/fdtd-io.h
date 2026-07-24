#ifndef _FDTD_IO_H
#define _FDTD_IO_H

#include <stddef.h>

void fdtdJoinPath(char *out, size_t outSize, const char *dir,
                  const char *filename);
int fdtdReadMatrixCsv(const char *filename, double *data, int rows, int cols);

#endif
