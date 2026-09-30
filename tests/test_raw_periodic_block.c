#include "cpmdc_params.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned char *slurp(const char *path, size_t *n) {
  FILE *f = fopen(path, "rb");
  if (!f)
    return NULL;
  if (fseek(f, 0, SEEK_END) != 0) {
    fclose(f);
    return NULL;
  }
  long len = ftell(f);
  if (len <= 0) {
    fclose(f);
    return NULL;
  }
  rewind(f);
  unsigned char *buf = malloc((size_t)len);
  if (!buf) {
    fclose(f);
    return NULL;
  }
  *n = fread(buf, 1, (size_t)len, f);
  fclose(f);
  return buf;
}

static int count_token(const char *deck, const char *token) {
  int n = 0;
  size_t len = strlen(token);
  for (const char *p = deck; (p = strstr(p, token)) != NULL; p += len)
    n++;
  return n;
}

/* A raw block that already names the periodic nitride deck is the deck.
 * The scalar fields still fill an absent section, and must not open a
 * second &SYSTEM with the isolated Hockney solver. */
int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "params_raw_periodic_block.bin";
  size_t sz = 0;
  unsigned char *msg = slurp(path, &sz);
  if (!msg) {
    fprintf(stderr, "missing %s\n", path);
    return 1;
  }
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(msg, sz, &arena, &root) != 0) {
    fprintf(stderr, "decode failed\n");
    free(msg);
    return 1;
  }
  char deck[CPMDC_BLOCKS];
  int rc = cpmdc_params_render_input_deck(root, deck, sizeof(deck));
  cpmdc_params_release(&arena);
  free(msg);
  if (rc != 0) {
    fprintf(stderr, "render failed\n");
    return 1;
  }
  if (count_token(deck, "&SYSTEM") != 1 || count_token(deck, "&CPMD") != 1 ||
      count_token(deck, "&DFT") != 1 || count_token(deck, "&ATOMS") != 1) {
    fprintf(stderr, "a raw section was opened again\n%s\n", deck);
    return 1;
  }
  if (strstr(deck, " SYMMETRY\n  1\n") == NULL ||
      strstr(deck, " SYMMETRY\n  0\n") != NULL ||
      strstr(deck, "POISSON SOLVER HOCKNEY") != NULL) {
    fprintf(stderr, "periodic system was replaced\n%s\n", deck);
    return 1;
  }
  if (strstr(deck, "NEWCODE") == NULL || strstr(deck, "GC-CUTOFF") == NULL ||
      strstr(deck, "*Si_MT_BLYP.psp KLEINMAN-BYLANDER") == NULL ||
      count_token(deck, "*Si_MT_BLYP.psp") != 1) {
    fprintf(stderr, "group keywords did not survive\n%s\n", deck);
    return 1;
  }
  return 0;
}
