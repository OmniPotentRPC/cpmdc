#include "cpmdc_params.h"
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

static int render_deck(const char *binpath, char *deck, size_t decksz) {
  size_t sz = 0;
  unsigned char *msg = slurp(binpath, &sz);
  if (!msg) return -1;
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(msg, sz, &arena, &root) != 0) {
    free(msg);
    return -1;
  }
  int rc = cpmdc_params_render_input_deck(root, deck, decksz);
  cpmdc_params_release(&arena);
  free(msg);
  return rc;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params50.bin params70.bin\n", argv[0]);
    return 2;
  }
  char d50[CPMDC_BLOCKS], d70[CPMDC_BLOCKS];
  if (render_deck(argv[1], d50, sizeof(d50)) != 0 ||
      render_deck(argv[2], d70, sizeof(d70)) != 0) {
    fprintf(stderr, "render failed\n");
    return 1;
  }
  int has50 = strstr(d50, "50") != NULL;
  int has70 = strstr(d70, "70") != NULL;
  int differ = strcmp(d50, d70) != 0;
  printf("deck50_has_50=%d deck70_has_70=%d decks_differ=%d\n", has50, has70,
         differ);
  if (strstr(d50, "inputBlocks") || strstr(d50, "&") == NULL) {
    /* still require CUTOFF present */
  }
  printf("deck50_snip_CUTOFF: ");
  const char *p = strstr(d50, "CUTOFF");
  if (p)
    fwrite(p, 1, 40 < strlen(p) ? 40 : strlen(p), stdout);
  printf("\n");
  if (!has50 || !has70 || !differ)
    return 1;
  return 0;
}
