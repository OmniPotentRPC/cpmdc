/* A scalar CPMDParams message names no system section. The cell-less deck
 * leaves symmetry unset. A force call with a positive volume chooses cubic
 * or orthorhombic symmetry, or cell vectors when the box is tilted. A zero
 * box stays an isolated Hockney cluster. */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);

static int g_failed = 0;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  long n;
  unsigned char *buf;
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  buf = (unsigned char *)malloc((size_t)n);
  if (!buf) {
    fclose(fp);
    return NULL;
  }
  if (fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static void require_has(const char *deck, const char *needle, const char *what) {
  if (deck && needle && strstr(deck, needle))
    return;
  fprintf(stderr, "%s: missing [%s]\n%s\n", what, needle ? needle : "",
          deck ? deck : "");
  g_failed = 1;
}

static void require_absent(const char *deck, const char *needle,
                           const char *what) {
  if (!deck || !needle || !strstr(deck, needle))
    return;
  fprintf(stderr, "%s: unexpected [%s]\n%s\n", what, needle, deck);
  g_failed = 1;
}

static int render_case(CPMDParams_ptr root, const double *cell, int has_cell,
                       char *deck, size_t cap) {
  double pos[6] = {0.0, 0.0, 0.0, 0.74, 0.0, 0.0};
  int z[2] = {1, 1};
  return cpmdc_params_render_deck_with_geometry(root, 2, pos, z, cell, has_cell,
                                                deck, cap);
}

static int compose_case(const double *cell, int has_cell, char *deck,
                        int cap) {
  double pos[6] = {0.0, 0.0, 0.0, 0.74, 0.0, 0.0};
  int z[2] = {1, 1};
  int deck_len = 0;
  int rc = cpmdc_embed_compose_cold_deck(2, pos, z, cell, has_cell, deck, cap,
                                         &deck_len);
  if (rc != 1 || deck_len < 1) {
    fprintf(stderr, "compose failed rc=%d len=%d has_cell=%d\n%s\n", rc,
            deck_len, has_cell, deck);
    g_failed = 1;
    return -1;
  }
  deck[deck_len] = '\0';
  return 0;
}

int main(int argc, char **argv) {
  size_t n = 0;
  unsigned char *msg;
  struct capn arena;
  CPMDParams_ptr root;
  char deck[CPMDC_BLOCKS];
  double orth[9] = {10.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 15.0};
  double cube[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  double tilt[9] = {10.0, 0.5, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0, 14.0};
  double zero[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if (argc < 2) {
    fprintf(stderr, "usage: %s params_scalar_blyp.bin\n", argv[0]);
    return 2;
  }
  msg = read_file(argv[1], &n);
  if (!msg) {
    fprintf(stderr, "failed to read %s\n", argv[1]);
    return 1;
  }
  if (cpmdc_params_root(msg, n, &arena, &root) != 0) {
    fprintf(stderr, "params parse failed\n");
    free(msg);
    return 1;
  }

  memset(deck, 0, sizeof(deck));
  if (cpmdc_params_render_input_deck(root, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "render_input_deck failed\n");
    g_failed = 1;
  } else {
    require_has(deck, "&SYSTEM\n", "cell-less render");
    require_absent(deck, "SYMMETRY", "cell-less render");
    require_absent(deck, "HOCKNEY", "cell-less render");
  }

  memset(deck, 0, sizeof(deck));
  if (render_case(root, orth, 1, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "orthorhombic render failed\n");
    g_failed = 1;
  } else {
    require_has(deck, " SYMMETRY\n  8\n", "orthorhombic render");
    require_absent(deck, "HOCKNEY", "orthorhombic render");
    require_absent(deck, " SYMMETRY\n  0\n", "orthorhombic render");
    require_absent(deck, " SYMMETRY\n  1\n", "orthorhombic render");
  }

  memset(deck, 0, sizeof(deck));
  if (render_case(root, cube, 1, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "cubic render failed\n");
    g_failed = 1;
  } else {
    require_has(deck, " SYMMETRY\n  1\n", "cubic render");
    require_absent(deck, "HOCKNEY", "cubic render");
    require_absent(deck, " SYMMETRY\n  0\n", "cubic render");
    require_absent(deck, " SYMMETRY\n  8\n", "cubic render");
  }

  memset(deck, 0, sizeof(deck));
  if (render_case(root, tilt, 1, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "tilted render failed\n");
    g_failed = 1;
  } else {
    require_has(deck, "CELL VECTORS", "tilted render");
    require_absent(deck, "SYMMETRY", "tilted render");
    require_absent(deck, "HOCKNEY", "tilted render");
  }

  memset(deck, 0, sizeof(deck));
  if (render_case(root, zero, 1, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "zero-cell render failed\n");
    g_failed = 1;
  } else {
    require_has(deck, " SYMMETRY\n  0\n", "zero-cell render");
    require_has(deck, "POISSON SOLVER HOCKNEY", "zero-cell render");
  }

  cpmdc_params_release(&arena);

  if (cpmdc_set_params(msg, n) != 0) {
    fprintf(stderr, "set_params failed\n");
    free(msg);
    return 1;
  }
  free(msg);

  memset(deck, 0, sizeof(deck));
  if (compose_case(orth, 1, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, " SYMMETRY\n  8\n", "orthorhombic compose");
    require_absent(deck, "HOCKNEY", "orthorhombic compose");
    require_absent(deck, " SYMMETRY\n  0\n", "orthorhombic compose");
    require_absent(deck, " SYMMETRY\n  1\n", "orthorhombic compose");
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(zero, 0, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, " SYMMETRY\n  0\n", "missing-box compose");
    require_has(deck, "POISSON SOLVER HOCKNEY", "missing-box compose");
    require_has(deck, "12.000000", "missing-box compose");
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(tilt, 1, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, "CELL VECTORS", "tilted compose");
    require_absent(deck, "SYMMETRY", "tilted compose");
    require_absent(deck, "HOCKNEY", "tilted compose");
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(cube, 1, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, " SYMMETRY\n  1\n", "cubic compose");
    require_absent(deck, "HOCKNEY", "cubic compose");
    require_absent(deck, " SYMMETRY\n  0\n", "cubic compose");
    require_absent(deck, " SYMMETRY\n  8\n", "cubic compose");
  }

  return g_failed ? 1 : 0;
}
