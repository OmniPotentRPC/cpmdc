/* Reads cells, nine reals per line in strtod syntax, and prints
 * cpmdc_cell_rank3_state() for each, so a test can compare it with
 * cpmdc._rank on the same bits. */
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
    printf("%d\n", cpmdc_cell_rank3_state(cell));
  }
  return 0;
}
