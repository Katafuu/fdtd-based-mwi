#include "fdtd-io.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fdtdJoinPath(char *out, size_t outSize, const char *dir,
                  const char *filename) {
  size_t len;
  char sep = '/';

  if (dir == NULL || dir[0] == '\0') {
    snprintf(out, outSize, "%s", filename);
    return;
  }

  len = strlen(dir);
  if (dir[len - 1] == '/' || dir[len - 1] == '\\')
    sep = '\0';

  if (sep == '\0')
    snprintf(out, outSize, "%s%s", dir, filename);
  else
    snprintf(out, outSize, "%s%c%s", dir, sep, filename);
}

int fdtdReadMatrixCsv(const char *filename, double *data, int rows, int cols) {
  FILE *fp;
  char line[65536];
  int row = 0;
  int col = 0;

  if (rows < 1 || cols < 1) {
    fprintf(stderr, "fdtdReadMatrixCsv: matrix dimensions must be positive.\n");
    return 0;
  }

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "fdtdReadMatrixCsv: could not open %s.\n", filename);
    return 0;
  }

  while (fgets(line, sizeof(line), fp) != NULL) {
    char *ptr = line;

    while (*ptr != '\0') {
      char *end;
      double value;

      while (*ptr == ',' || *ptr == ' ' || *ptr == '\t' || *ptr == '\r')
        ptr++;

      if (*ptr == '\n') {
        if (col != cols) {
          fprintf(stderr,
                  "fdtdReadMatrixCsv: row %d in %s has %d values; expected %d.\n",
                  row + 1, filename, col, cols);
          fclose(fp);
          return 0;
        }
        row++;
        col = 0;
        ptr++;
        continue;
      }

      if (*ptr == '\0')
        break;

      if (row >= rows || col >= cols) {
        fprintf(stderr,
                "fdtdReadMatrixCsv: %s has more than %d-by-%d values.\n",
                filename, rows, cols);
        fclose(fp);
        return 0;
      }

      value = strtod(ptr, &end);
      if (end == ptr) {
        fprintf(stderr, "fdtdReadMatrixCsv: invalid numeric value in %s.\n",
                filename);
        fclose(fp);
        return 0;
      }

      data[(size_t)row * (size_t)cols + (size_t)col] = value;
      col++;
      ptr = end;
    }
  }

  if (ferror(fp)) {
    fprintf(stderr, "fdtdReadMatrixCsv: error while reading %s.\n", filename);
    fclose(fp);
    return 0;
  }
  fclose(fp);

  if (col > 0) {
    if (col != cols) {
      fprintf(stderr,
              "fdtdReadMatrixCsv: row %d in %s has %d values; expected %d.\n",
              row + 1, filename, col, cols);
      return 0;
    }
    row++;
  }

  if (row != rows) {
    fprintf(stderr,
            "fdtdReadMatrixCsv: %s has %d rows; expected %d.\n",
            filename, row, rows);
    return 0;
  }

  return 1;
}
