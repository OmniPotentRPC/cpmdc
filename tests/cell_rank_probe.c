/* Reads cells, nine reals per line in strtod syntax, and prints for each
 * "state rc volume": cpmdc_cell_rank3_state(), then cpmdc_cell_volume3()'s
 * return code and volume in %a, so a test can compare the bits. */
#include "cpmdc_cell_rank.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  char line[1024];
  while (fgets(line, sizeof(line), stdin)) {
    double cell[9];
    char *p = line;
    for (int i = 0; i < 9; ++i) {
      char *end = NULL;
      cell[i] = strtod(p, &end);
      if (end == p) {
        fprintf(stderr, "cell_rank_probe: cannot read entry %d of: %s", i, line);
        return 2;
      }
      p = end;
    }
    double volume = 0.0;
    int rc = cpmdc_cell_volume3(cell, &volume);
    printf("%d %d %a\n", cpmdc_cell_rank3_state(cell), rc, volume);
  }
  return 0;
}
