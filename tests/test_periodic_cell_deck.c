/* A scalar CPMDParams message names no system section. The cell-less deck
 * leaves symmetry unset. A force call with a positive volume chooses cubic
 * or orthorhombic symmetry, or cell vectors when the box is tilted. A zero
 * box stays an isolated Hockney cluster. */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <ctype.h>
#include <float.h>
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

static void require_lines_within_80(const char *deck, const char *what) {
  if (!deck)
    return;
  for (const char *line = deck; *line;) {
    const char *end = strchr(line, '\n');
    size_t n = end ? (size_t)(end - line) : strlen(line);
    if (n > 80) {
      fprintf(stderr, "%s: line length %zu exceeds 80: %.*s\n", what, n,
              (int)n, line);
      g_failed = 1;
    }
    line += n + (end ? 1 : 0);
  }
}

static void require_diagonal_cell_fields(const char *deck, const char *what) {
  const char *cell = deck ? strstr(deck, " CELL\n") : NULL;
  double parsed[6] = {0};
  if (!cell) {
    fprintf(stderr, "%s: missing CELL section\n%s\n", what,
            deck ? deck : "");
    g_failed = 1;
    return;
  }
  cell += strlen(" CELL\n");
  if (sscanf(cell, "%lf %lf %lf %lf %lf %lf", &parsed[0], &parsed[1],
             &parsed[2], &parsed[3], &parsed[4], &parsed[5]) != 6 ||
      parsed[0] <= 0.0 || parsed[1] <= 0.0 || parsed[2] <= 0.0 ||
      parsed[3] != 0.0 || parsed[4] != 0.0 || parsed[5] != 0.0) {
    fprintf(stderr, "%s: unexpected CELL fields\n%s\n", what, deck);
    g_failed = 1;
  }
}

/* CELL VECTORS must read back as the very cell that was accepted. */
static void require_vectors_preserved(const char *deck, const double *cell,
                                      const char *what) {
  const char *vectors = strstr(deck, "CELL VECTORS\n");
  double parsed[9];
  const char *cursor;
  if (!vectors) {
    fprintf(stderr, "%s: could not read CELL VECTORS\n%s\n", what, deck);
    g_failed = 1;
    return;
  }
  cursor = vectors + strlen("CELL VECTORS\n");
  for (int i = 0; i < 9; ++i) {
    char *end = NULL;
    parsed[i] = strtod(cursor, &end);
    if (end == cursor || !isspace((unsigned char)*end)) {
      fprintf(stderr, "%s: cell[%d] has no delimited number\n%s\n", what, i,
              deck);
      g_failed = 1;
      return;
    }
    cursor = end;
  }
  for (int i = 0; i < 9; ++i) {
    if (parsed[i] != cell[i]) {
      fprintf(stderr, "%s: cell[%d] changed from %.17g to %.17g\n%s\n", what,
              i, cell[i], parsed[i], deck);
      g_failed = 1;
    }
  }
  require_absent(deck, "SYMMETRY", what);
  require_absent(deck, "HOCKNEY", what);
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
  double small_tilt[9] = {1.0, 0.0, 0.0, 0.3, 2.58e-12, 0.0,
                          0.2, 0.1, 0.8};
  double zero[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  /* A strained cell and cells whose last vector approaches the rank cutoff. */
  const double accepted[][9] = {
      {10.0, -5.0e-5, 0.0, -5.0e-5, 10.0, 0.0, 0.0, 0.0, 10.0},
      {1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.125, 0.25, 2.58e-12},
      {1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.125, 0.25, 1.0e-13},
      {1.0, 0.0, 0.0, 0.25, 1.0, 0.0, 0.125, 0.25, 6.0 * DBL_EPSILON},
  };
  const size_t n_accepted = sizeof(accepted) / sizeof(accepted[0]);

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
    require_has(deck, "0.5", "tilted render");
    require_absent(deck, "SYMMETRY", "tilted render");
    require_absent(deck, "HOCKNEY", "tilted render");
  }

  memset(deck, 0, sizeof(deck));
  if (render_case(root, small_tilt, 1, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "small tilted render failed\n");
    g_failed = 1;
  } else {
    require_has(deck, "CELL VECTORS", "small tilted render");
    require_absent(deck, "SYMMETRY", "small tilted render");
    require_absent(deck, "HOCKNEY", "small tilted render");
  }

  for (size_t i = 0; i < n_accepted; ++i) {
    memset(deck, 0, sizeof(deck));
    if (render_case(root, accepted[i], 1, deck, sizeof(deck)) != 0) {
      fprintf(stderr, "accepted-cell render %zu failed\n", i);
      g_failed = 1;
    } else {
      require_vectors_preserved(deck, accepted[i], "accepted-cell render");
    }
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
    require_lines_within_80(deck, "orthorhombic compose");
    require_diagonal_cell_fields(deck, "orthorhombic compose");
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(small_tilt, 1, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, "CELL VECTORS", "small tilted compose");
    require_absent(deck, "SYMMETRY", "small tilted compose");
    require_absent(deck, "HOCKNEY", "small tilted compose");
    require_lines_within_80(deck, "small tilted compose");
  }

  for (size_t i = 0; i < n_accepted; ++i) {
    memset(deck, 0, sizeof(deck));
    if (compose_case(accepted[i], 1, deck, (int)sizeof(deck)) == 0) {
      require_vectors_preserved(deck, accepted[i], "accepted-cell compose");
      require_lines_within_80(deck, "accepted-cell compose");
    }
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(zero, 0, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, " SYMMETRY\n  0\n", "missing-box compose");
    require_has(deck, "POISSON SOLVER HOCKNEY", "missing-box compose");
    require_has(deck, "1.2000000000000000E+001", "missing-box compose");
    require_lines_within_80(deck, "missing-box compose");
    require_diagonal_cell_fields(deck, "missing-box compose");
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(tilt, 1, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, "CELL VECTORS", "tilted compose");
    require_has(deck, "5.0000000000000000E-001", "tilted compose");
    require_absent(deck, "SYMMETRY", "tilted compose");
    require_absent(deck, "HOCKNEY", "tilted compose");
    require_lines_within_80(deck, "tilted compose");
  }

  memset(deck, 0, sizeof(deck));
  if (compose_case(cube, 1, deck, (int)sizeof(deck)) == 0) {
    require_has(deck, " SYMMETRY\n  1\n", "cubic compose");
    require_absent(deck, "HOCKNEY", "cubic compose");
    require_absent(deck, " SYMMETRY\n  0\n", "cubic compose");
    require_absent(deck, " SYMMETRY\n  8\n", "cubic compose");
    require_lines_within_80(deck, "cubic compose");
    require_diagonal_cell_fields(deck, "cubic compose");
  }

  return g_failed ? 1 : 0;
}
