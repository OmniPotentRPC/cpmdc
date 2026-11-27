/* The renderer names the scheme on an EXTRAPOLATE WFN line explicitly.
 * A line that names none of ASPC, POLY, and POLYNOMIAL already selects
 * ASPC, matching CPMD's own default. STORE alone renders
 * EXTRAPOLATE WFN STORE ASPC. Order 0 writes that line with no order
 * line, and POLY at order 0 stays EXTRAPOLATE WFN POLY. */
#include "cpmdc_params.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static capn_text text_of(const char *s) {
  capn_text t = {s ? (int)strlen(s) : 0, s ? s : "", NULL};
  return t;
}

static int render(const char *options, int order, char *deck, size_t cap) {
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  CPMDParams_ptr ptr = new_CPMDParams(root.seg);
  struct CPMDParams view;
  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
  struct CPMDCpmdSection cpmd;
  CPMDCpmdSection_ptr cp = new_CPMDCpmdSection(root.seg);
  struct CPMDInputSection sec;
  int rc;

  memset(&view, 0, sizeof(view));
  view.task = text_of("gradient");
  view.cutOffRy = 70.0;
  view.multiplicity = 1;
  memset(&cpmd, 0, sizeof(cpmd));
  cpmd.optimizeWavefunction = 1;
  if (options)
    cpmd.extrapolateWfnOptions = text_of(options);
  cpmd.extrapolateWfnOrder = order;
  write_CPMDCpmdSection(&cpmd, cp);
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_cpmd;
  sec.cpmd = cp;
  set_CPMDInputSection(&sec, secs, 0);
  view.inputSections = secs;
  write_CPMDParams(&view, ptr);
  rc = cpmdc_params_render_input_deck(ptr, deck, cap);
  capn_free(&arena);
  return rc;
}

static int keyword_line_is(const char *deck, const char *want) {
  const char *hit = strstr(deck, " EXTRAPOLATE WFN");
  const char *nl;
  char line[256];
  size_t n;
  if (!hit)
    return 0;
  nl = strchr(hit, '\n');
  if (!nl)
    return 0;
  n = (size_t)(nl - hit);
  if (n >= sizeof(line))
    return 0;
  memcpy(line, hit, n);
  line[n] = '\0';
  return strcmp(line, want) == 0;
}

static int expect_line(const char *options, int order, const char *want_line,
                       const char *want_order) {
  char deck[CPMDC_BLOCKS];
  const char *order_at;
  if (render(options, order, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "render failed for options [%s] order %d\n",
            options ? options : "", order);
    return 1;
  }
  if (!keyword_line_is(deck, want_line)) {
    fprintf(stderr, "keyword line for [%s] order %d\nwant [%s]\ndeck:\n%s\n",
            options ? options : "", order, want_line, deck);
    return 1;
  }
  order_at = strstr(deck, want_line);
  if (!order_at || strncmp(order_at + strlen(want_line), want_order,
                           strlen(want_order)) != 0) {
    fprintf(stderr, "order line for [%s] order %d\nwant [%s]\ndeck:\n%s\n",
            options ? options : "", order, want_order, deck);
    return 1;
  }
  return 0;
}

/* An order payload is a line of two spaces and digits. Order 0 omits it. */
static int order_line_follows(const char *after) {
  const char *p;
  if (!after || after[0] != '\n')
    return 0;
  p = after + 1;
  if (p[0] != ' ' || p[1] != ' ')
    return 0;
  p += 2;
  if (*p < '0' || *p > '9')
    return 0;
  while (*p >= '0' && *p <= '9')
    p++;
  return *p == '\n';
}

static int expect_line_without_order(const char *options, const char *want_line) {
  char deck[CPMDC_BLOCKS];
  const char *at;
  if (render(options, 0, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "render failed for options [%s] order 0\n", options);
    return 1;
  }
  if (!keyword_line_is(deck, want_line)) {
    fprintf(stderr, "keyword line for [%s] order 0\nwant [%s]\ndeck:\n%s\n",
            options, want_line, deck);
    return 1;
  }
  at = strstr(deck, want_line);
  if (!at || order_line_follows(at + strlen(want_line))) {
    fprintf(stderr, "order line for [%s] order 0\ndeck:\n%s\n", options, deck);
    return 1;
  }
  return 0;
}

static int expect_no_extrapolation(const char *options, int order) {
  char deck[CPMDC_BLOCKS];
  if (render(options, order, deck, sizeof(deck)) != 0) {
    fprintf(stderr, "render failed for options [%s] order %d\n",
            options ? options : "", order);
    return 1;
  }
  if (strstr(deck, "EXTRAPOLATE") != NULL) {
    fprintf(stderr, "extrapolation line for [%s] order %d\ndeck:\n%s\n",
            options ? options : "", order, deck);
    return 1;
  }
  return 0;
}

int main(void) {
  if (expect_line(NULL, 4, " EXTRAPOLATE WFN ASPC", "\n  4\n") != 0)
    return 1;
  if (expect_line("STORE", 4, " EXTRAPOLATE WFN STORE ASPC", "\n  4\n") != 0)
    return 1;
  if (expect_line("POLY", 4, " EXTRAPOLATE WFN POLY", "\n  4\n") != 0)
    return 1;
  if (expect_line("POLYNOMIAL", 3, " EXTRAPOLATE WFN POLYNOMIAL", "\n  3\n") !=
      0)
    return 1;
  if (expect_line("ASPC", 4, " EXTRAPOLATE WFN ASPC", "\n  4\n") != 0)
    return 1;
  if (expect_line("STORE ASPC CSTEPS=3", 5,
                  " EXTRAPOLATE WFN STORE ASPC CSTEPS=3", "\n  5\n") != 0)
    return 1;
  if (expect_no_extrapolation(NULL, 0) != 0)
    return 1;
  if (expect_line_without_order("STORE", " EXTRAPOLATE WFN STORE ASPC") != 0)
    return 1;
  if (expect_line_without_order("POLY", " EXTRAPOLATE WFN POLY") != 0)
    return 1;
  return 0;
}
