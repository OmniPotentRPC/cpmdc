#include "cpmdc_params.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int periodic_vector_deck(const char *deck) {
  if (!strstr(deck, "CELL VECTORS") || strstr(deck, "\n SYMMETRY\n") ||
      strstr(deck, "HOCKNEY")) {
    fprintf(stderr, "CELL VECTORS must omit SYMMETRY and isolated defaults\n%s", deck);
    return 0;
  }
  return 1;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  FILE *file = fopen(argv[1], "rb");
  if (!file) return 2;
  if (fseek(file, 0, SEEK_END) != 0) { fclose(file); return 2; }
  long count = ftell(file);
  if (count <= 0 || count > (1 << 20)) { fclose(file); return 2; }
  rewind(file);
  unsigned char *bytes = malloc((size_t)count);
  if (!bytes) { fclose(file); return 2; }
  int read_ok = fread(bytes, 1, (size_t)count, file) == (size_t)count;
  fclose(file);
  if (!read_ok) { free(bytes); return 2; }
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(bytes, (size_t)count, &arena, &root) != 0) {
    free(bytes); return 2;
  }
  char deck[CPMDC_BLOCKS];
  int ok = cpmdc_params_render_input_deck(root, deck, sizeof(deck)) == 0 &&
           periodic_vector_deck(deck) && strstr(deck, "1 0 0 0 2 0 0 0 3");
  const double positions[6] = {0, 0, 0, .74, 0, 0};
  const int numbers[2] = {1, 1};
  const double cell[9] = {10, .5, 0, 0, 12, 0, 0, 0, 14};
  if (cpmdc_params_render_deck_with_geometry(root, 2, positions, numbers,
                                            cell, 1, deck, sizeof(deck)) != 0 ||
      !periodic_vector_deck(deck)) ok = 0;
  cpmdc_params_release(&arena);
  free(bytes);
  return ok ? 0 : 1;
}
