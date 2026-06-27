#include "cpmdc_params.h"
#include "cpmdc_features.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned char *slurp(const char *p, size_t *n) {
  FILE *f = fopen(p, "rb");
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);
  long L = ftell(f);
  rewind(f);
  unsigned char *b = malloc((size_t)L);
  *n = fread(b, 1, (size_t)L, f);
  fclose(f);
  return b;
}

static int check_deck(const char *bin, const char **need, int nneed) {
  size_t sz = 0;
  unsigned char *msg = slurp(bin, &sz);
  if (!msg) { fprintf(stderr, "missing %s\n", bin); return -1; }
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(msg, sz, &arena, &root) != 0) { free(msg); return -1; }
  char deck[CPMDC_BLOCKS];
  if (cpmdc_params_render_input_deck(root, deck, sizeof(deck)) != 0) {
    cpmdc_params_release(&arena);
    free(msg);
    return -1;
  }
  for (int i = 0; i < nneed; ++i) {
    if (!strstr(deck, need[i])) {
      fprintf(stderr, "missing token [%s] in deck from %s\n", need[i], bin);
      fprintf(stderr, "--- deck ---\n%s\n", deck);
      cpmdc_params_release(&arena);
      free(msg);
      return -1;
    }
  }
  cpmdc_params_release(&arena);
  free(msg);
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s cp_md.bin dft_func.bin cpmd_geometry.bin\n", argv[0]);
    return 2;
  }
  /* inventory must resolve CP and DFT catalog ids */
  if (!cpmdc_feature_find("catalog.section.PIMD") ||
      !cpmdc_feature_find("catalog.section.VDW") ||
      !cpmdc_feature_find("catalog.cpmd.EMASS") ||
      !cpmdc_feature_find("catalog.dft.FUNCTIONAL_BLYP") ||
      !cpmdc_feature_find("catalog.dft.FUNCTIONAL_PBE0")) {
    fprintf(stderr, "inventory find failed for CP/DFT catalog ids\n");
    return 1;
  }
  const char *cp_need[] = {
      "&CPMD", "MOLECULAR DYNAMICS", "EMASS", "400.0", "NOSE IONS",
      "ANNEALING IONS", "TIMESTEP", "RESTART WAVEFUNCTION", "TRAJECTORY",
      "&SYSTEM", "SCALE", "0.5", "CHARGE", "1", "&DFT", "FUNCTIONAL BLYP",
  };
  const char *dft_need[] = {
      "&DFT", "FUNCTIONAL PBE0", "LSD", "GC-CUTOFF", "&PIMD", "TROTTER",
      "&VDW", "VDW CORRECTION",
  };
  const char *geometry_need[] = {
      "&CPMD", "OPTIMIZE GEOMETRY", "CONVERGENCE ORBITALS", "2e-06",
      "CONVERGENCE GEOMETRY", "0.0001", "MAXSTEP", "50", "MAXITER", "12",
      "TIMESTEP", "3", "EMASS", "450", "&SYSTEM", "CUTOFF", "80", "&DFT",
      "FUNCTIONAL PBE",
  };
  if (check_deck(argv[1], cp_need, (int)(sizeof(cp_need) / sizeof(cp_need[0]))) != 0)
    return 1;
  if (check_deck(argv[2], dft_need, (int)(sizeof(dft_need) / sizeof(dft_need[0]))) != 0)
    return 1;
  if (check_deck(argv[3], geometry_need,
                 (int)(sizeof(geometry_need) / sizeof(geometry_need[0]))) != 0)
    return 1;
  printf("OK cp_md, dft_multi, and cpmd_geometry render + inventory finds\n");
  return 0;
}
