/* The effective functional, cutoff, charge, and multiplicity of a CPMDParams
 * message with scalar overrides: both deck renderers print exactly the values
 * cpmdc_params_effective_config_ov reports. A value a typed field sets beats
 * the override, the override beats a schema default, and a typed section
 * beats the top-level field. A typed field set to its schema default (BLYP,
 * 70 Ry, charge 0, multiplicity 1) reads as unset. */
#include "cpmdc_params.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static capn_text text_of(const char *s) {
  capn_text t = {s ? (int)strlen(s) : 0, s ? s : "", NULL};
  return t;
}

struct Case {
  const char *name;
  const char *root_functional;
  double root_cutoff_ry;
  int root_charge;
  int root_mult;
  int has_dft;
  const char *dft_functional;
  int has_sys;
  double sys_cutoff_ry;
  int sys_charge;
  int sys_mult;
  const char *ov_functional;
  double ov_cutoff_ry;
  int ov_charge;
  int ov_mult;
  const char *want_functional;
  double want_cutoff_ry;
  int want_charge;
  int want_mult;
};

struct Scalars {
  char functional[64];
  double cutoff_ry;
  int charge;
  int mult;
};

static int count_text(const char *hay, const char *needle) {
  int n = 0;
  size_t len = strlen(needle);
  for (const char *p = hay; (p = strstr(p, needle)) != NULL; p += len)
    n++;
  return n;
}

/* The scalars a deck states. CPMD defaults are 0 for CHARGE and 1 for
 * MULTIPLICITY when the deck omits them. A keyword that appears twice is a
 * failure: the second line would hide behind the first. */
static int read_deck(const char *deck, struct Scalars *out) {
  memset(out, 0, sizeof(*out));
  if (count_text(deck, " FUNCTIONAL ") != 1 ||
      count_text(deck, " CUTOFF\n  ") != 1 ||
      count_text(deck, " CHARGE\n  ") > 1 ||
      count_text(deck, " MULTIPLICITY\n  ") > 1)
    return -1;
  const char *fn = strstr(deck, " FUNCTIONAL ") + strlen(" FUNCTIONAL ");
  size_t n = strcspn(fn, "\n");
  if (n >= sizeof(out->functional))
    return -1;
  memcpy(out->functional, fn, n);
  out->functional[n] = '\0';
  out->cutoff_ry =
      strtod(strstr(deck, " CUTOFF\n  ") + strlen(" CUTOFF\n  "), NULL);
  const char *ch = strstr(deck, " CHARGE\n  ");
  out->charge = ch ? (int)strtol(ch + strlen(" CHARGE\n  "), NULL, 10) : 0;
  const char *mu = strstr(deck, " MULTIPLICITY\n  ");
  out->mult = mu ? (int)strtol(mu + strlen(" MULTIPLICITY\n  "), NULL, 10) : 1;
  return 0;
}

static int same(const struct Scalars *a, const struct Scalars *b) {
  return strcmp(a->functional, b->functional) == 0 &&
         fabs(a->cutoff_ry - b->cutoff_ry) < 1e-9 && a->charge == b->charge &&
         a->mult == b->mult;
}

static void build(struct capn *arena, const struct Case *c,
                  CPMDParams_ptr *out) {
  capn_init_malloc(arena);
  capn_ptr root = capn_root(arena);
  CPMDParams_ptr ptr = new_CPMDParams(root.seg);
  struct CPMDParams view;
  memset(&view, 0, sizeof(view));
  view.task = text_of("gradient");
  if (c->root_functional)
    view.functional = text_of(c->root_functional);
  view.cutOffRy = c->root_cutoff_ry > 0.0 ? c->root_cutoff_ry : 70.0;
  view.charge = c->root_charge;
  view.multiplicity = c->root_mult > 0 ? c->root_mult : 1;
  int nsec = (c->has_dft ? 1 : 0) + (c->has_sys ? 1 : 0);
  if (nsec > 0) {
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, nsec);
    int idx = 0;
    if (c->has_dft) {
      struct CPMDDftSection dft;
      memset(&dft, 0, sizeof(dft));
      if (c->dft_functional)
        dft.functional = text_of(c->dft_functional);
      CPMDDftSection_ptr dp = new_CPMDDftSection(root.seg);
      write_CPMDDftSection(&dft, dp);
      struct CPMDInputSection sec;
      memset(&sec, 0, sizeof(sec));
      sec.which = CPMDInputSection_dft;
      sec.dft = dp;
      set_CPMDInputSection(&sec, secs, idx++);
    }
    if (c->has_sys) {
      struct CPMDSystemSection sys;
      memset(&sys, 0, sizeof(sys));
      sys.angstrom = 1;
      sys.cutOffRy = c->sys_cutoff_ry;
      sys.charge = c->sys_charge;
      sys.multiplicity = c->sys_mult;
      CPMDSystemSection_ptr sp = new_CPMDSystemSection(root.seg);
      write_CPMDSystemSection(&sys, sp);
      struct CPMDInputSection sec;
      memset(&sec, 0, sizeof(sec));
      sec.which = CPMDInputSection_system;
      sec.system = sp;
      set_CPMDInputSection(&sec, secs, idx++);
    }
    view.inputSections = secs;
  }
  write_CPMDParams(&view, ptr);
  *out = ptr;
}

static int check_case(const struct Case *c) {
  struct capn arena;
  CPMDParams_ptr params;
  build(&arena, c, &params);

  CPMDCScalarOverrides ov;
  memset(&ov, 0, sizeof(ov));
  ov.functional = c->ov_functional;
  ov.cutoff_ry = c->ov_cutoff_ry;
  if (c->ov_charge != 0) {
    ov.has_charge = 1;
    ov.charge = c->ov_charge;
  }
  if (c->ov_mult != 0) {
    ov.has_multiplicity = 1;
    ov.multiplicity = c->ov_mult;
  }

  struct Scalars want;
  memset(&want, 0, sizeof(want));
  snprintf(want.functional, sizeof(want.functional), "%s", c->want_functional);
  want.cutoff_ry = c->want_cutoff_ry;
  want.charge = c->want_charge;
  want.mult = c->want_mult;

  struct Scalars config;
  memset(&config, 0, sizeof(config));
  int rc = cpmdc_params_effective_config_ov(
      params, &ov, config.functional, sizeof(config.functional),
      &config.cutoff_ry, &config.charge, &config.mult);

  char input_deck[CPMDC_BLOCKS];
  char geometry_deck[CPMDC_BLOCKS];
  struct Scalars from_input, from_geometry;
  const double positions[6] = {0.0, 0.0, 0.0, 0.757, 0.586, 0.0};
  const int atomic_numbers[2] = {8, 1};
  int rc_input = cpmdc_params_render_input_deck_ov(params, &ov, input_deck,
                                                   sizeof(input_deck));
  int rc_geometry = cpmdc_params_render_deck_with_geometry_ov(
      params, &ov, 2, positions, atomic_numbers, NULL, 0, geometry_deck,
      sizeof(geometry_deck));
  int failures = 0;
  if (rc != 0 || rc_input != 0 || rc_geometry != 0 ||
      read_deck(input_deck, &from_input) != 0 ||
      read_deck(geometry_deck, &from_geometry) != 0 ||
      !same(&config, &want) || !same(&from_input, &want) ||
      !same(&from_geometry, &want)) {
    fprintf(stderr,
            "%s: want %s %.10g %d %d | config %s %.10g %d %d | input deck %s "
            "%.10g %d %d | geometry deck %s %.10g %d %d (rc %d %d %d)\n",
            c->name, want.functional, want.cutoff_ry, want.charge, want.mult,
            config.functional, config.cutoff_ry, config.charge, config.mult,
            from_input.functional, from_input.cutoff_ry, from_input.charge,
            from_input.mult, from_geometry.functional,
            from_geometry.cutoff_ry, from_geometry.charge, from_geometry.mult,
            rc, rc_input, rc_geometry);
    failures = 1;
  }
  capn_free(&arena);
  return failures;
}

static void test_overrides_fill_typed_defaults(void **state) {
  (void)state;
  const struct Case cases[] = {
      {.name = "dft at default text",
       .has_dft = 1, .dft_functional = "BLYP", .ov_functional = "PBE",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "dft unset text",
       .has_dft = 1, .ov_functional = "PBE",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "system at default 70",
       .has_sys = 1, .sys_cutoff_ry = 70.0, .ov_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "system cutoff unset",
       .has_sys = 1, .ov_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "dft and system at defaults",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 70.0, .ov_functional = "PBE", .ov_cutoff_ry = 80.0,
       .want_functional = "PBE", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "all four from the override",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 70.0, .sys_mult = 1, .ov_functional = "PBE",
       .ov_cutoff_ry = 50.0, .ov_charge = -1, .ov_mult = 3,
       .want_functional = "PBE", .want_cutoff_ry = 50.0, .want_charge = -1,
       .want_mult = 3},
      {.name = "no sections, override only",
       .ov_functional = "PBE", .ov_cutoff_ry = 80.0, .ov_charge = 2,
       .ov_mult = 2,
       .want_functional = "PBE", .want_cutoff_ry = 80.0, .want_charge = 2,
       .want_mult = 2},
      {.name = "no override, defaults",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 70.0,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_mult = 1},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += check_case(&cases[i]);
  assert_int_equal(failures, 0);
}

static void test_typed_value_beats_override(void **state) {
  (void)state;
  const struct Case cases[] = {
      {.name = "typed dft functional",
       .has_dft = 1, .dft_functional = "LDA", .ov_functional = "PBE",
       .want_functional = "LDA", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "typed system cutoff",
       .has_sys = 1, .sys_cutoff_ry = 50.0, .ov_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 50.0, .want_mult = 1},
      {.name = "typed system charge and multiplicity",
       .has_sys = 1, .sys_charge = 4, .sys_mult = 2, .ov_charge = -1,
       .ov_mult = 3,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_charge = 4,
       .want_mult = 2},
      {.name = "typed dft, default system",
       .has_dft = 1, .dft_functional = "LDA", .has_sys = 1,
       .sys_cutoff_ry = 70.0, .ov_functional = "PBE", .ov_cutoff_ry = 80.0,
       .want_functional = "LDA", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "default dft, typed system cutoff",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 50.0, .ov_functional = "PBE", .ov_cutoff_ry = 80.0,
       .want_functional = "PBE", .want_cutoff_ry = 50.0, .want_mult = 1},
      {.name = "typed dft and system",
       .has_dft = 1, .dft_functional = "LDA", .has_sys = 1,
       .sys_cutoff_ry = 50.0, .sys_charge = 1, .sys_mult = 2,
       .ov_functional = "PBE", .ov_cutoff_ry = 80.0, .ov_charge = -1,
       .ov_mult = 3,
       .want_functional = "LDA", .want_cutoff_ry = 50.0, .want_charge = 1,
       .want_mult = 2},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += check_case(&cases[i]);
  assert_int_equal(failures, 0);
}

static void test_top_level_value_with_default_sections(void **state) {
  (void)state;
  const struct Case cases[] = {
      {.name = "root functional, default dft",
       .root_functional = "PBE", .has_dft = 1, .dft_functional = "BLYP",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "root functional beats override",
       .root_functional = "PBE", .has_dft = 1, .dft_functional = "BLYP",
       .ov_functional = "LDA",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "root cutoff beats override",
       .root_cutoff_ry = 50.0, .has_sys = 1, .sys_cutoff_ry = 70.0,
       .ov_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 50.0, .want_mult = 1},
      {.name = "root charge and multiplicity beat override",
       .root_charge = 2, .root_mult = 3, .has_sys = 1, .ov_charge = -1,
       .ov_mult = 2,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_charge = 2,
       .want_mult = 3},
      {.name = "root multiplicity, default system",
       .root_mult = 3, .has_sys = 1,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_mult = 3},
      {.name = "typed section beats the root field",
       .root_functional = "PBE", .root_cutoff_ry = 50.0, .has_dft = 1,
       .dft_functional = "LDA", .has_sys = 1, .sys_cutoff_ry = 60.0,
       .want_functional = "LDA", .want_cutoff_ry = 60.0, .want_mult = 1},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += check_case(&cases[i]);
  assert_int_equal(failures, 0);
}

/* A CUTOFF, CHARGE, or MULTIPLICITY directive in the typed &SYSTEM section,
 * or a SYSTEM.<keyword> set directive, would follow the line the renderer
 * writes for that scalar. CUTOFF is written on every deck, including the
 * default 70 Ry. CHARGE is written only when it is not 0 and MULTIPLICITY
 * only when it is greater than 1, so a default leaves the directive as the
 * only line. Both are rejected: the directive is a second value when the
 * line is written, and a value the recorded configuration does not carry
 * when the line is omitted. No overlay is set. */
struct SystemDirectiveCase {
  const char *name;
  int has_sys;
  double cutoff_ry;
  int charge;
  int mult;
  const char *directive;
  const char *set_key;
  const char *set_value;
  const char *keyword;
  const char *reason;
};

static void build_system_directive(struct capn *arena,
                                   const struct SystemDirectiveCase *c,
                                   CPMDParams_ptr *out) {
  capn_init_malloc(arena);
  capn_ptr root = capn_root(arena);
  CPMDParams_ptr ptr = new_CPMDParams(root.seg);
  struct CPMDParams view;
  memset(&view, 0, sizeof(view));
  view.task = text_of("gradient");
  int nsec = (c->has_sys ? 1 : 0) + (c->set_key ? 1 : 0);
  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, nsec);
  int idx = 0;
  if (c->has_sys) {
    struct CPMDSystemSection sys;
    memset(&sys, 0, sizeof(sys));
    sys.angstrom = 1;
    sys.cutOffRy = c->cutoff_ry;
    sys.charge = c->charge;
    sys.multiplicity = c->mult > 0 ? c->mult : 1;
    if (c->directive) {
      CPMDDirective_list dirs = new_CPMDDirective_list(root.seg, 1);
      struct CPMDDirective d;
      memset(&d, 0, sizeof(d));
      d.keyword = text_of(c->directive);
      set_CPMDDirective(&d, dirs, 0);
      sys.directives = dirs;
    }
    CPMDSystemSection_ptr sp = new_CPMDSystemSection(root.seg);
    write_CPMDSystemSection(&sys, sp);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sp;
    set_CPMDInputSection(&sec, secs, idx++);
  }
  if (c->set_key) {
    struct CPMDSetDirective setd;
    memset(&setd, 0, sizeof(setd));
    setd.key = text_of(c->set_key);
    if (c->set_value)
      setd.value = text_of(c->set_value);
    CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
    write_CPMDSetDirective(&setd, set_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_set;
    sec.set = set_ptr;
    set_CPMDInputSection(&sec, secs, idx++);
  }
  view.inputSections = secs;
  write_CPMDParams(&view, ptr);
  *out = ptr;
}

static int expect_system_directive_rejected(const struct SystemDirectiveCase *c) {
  struct capn arena;
  CPMDParams_ptr params;
  build_system_directive(&arena, c, &params);
  char input_deck[CPMDC_BLOCKS];
  char geometry_deck[CPMDC_BLOCKS];
  const double positions[6] = {0.0, 0.0, 0.0, 0.757, 0.586, 0.0};
  const int atomic_numbers[2] = {8, 1};
  int rc_input = cpmdc_params_render_input_deck_ov(params, NULL, input_deck,
                                                   sizeof(input_deck));
  char reason_input[256];
  snprintf(reason_input, sizeof(reason_input), "%s",
           cpmdc_params_render_reason());
  int rc_geometry = cpmdc_params_render_deck_with_geometry_ov(
      params, NULL, 2, positions, atomic_numbers, NULL, 0, geometry_deck,
      sizeof(geometry_deck));
  char reason_geometry[256];
  snprintf(reason_geometry, sizeof(reason_geometry), "%s",
           cpmdc_params_render_reason());
  int input_ok = rc_input != 0 && strstr(reason_input, c->reason) != NULL;
  int geometry_ok =
      rc_geometry != 0 && strstr(reason_geometry, c->reason) != NULL;
  int failures = 0;
  if (!input_ok || !geometry_ok) {
    fprintf(stderr,
            "%s: input rc %d reason \"%s\" (keyword x%d); geometry rc %d "
            "reason \"%s\" (keyword x%d); want \"%s\"\n",
            c->name, rc_input, reason_input,
            rc_input == 0 ? count_text(input_deck, c->keyword) : -1,
            rc_geometry, reason_geometry,
            rc_geometry == 0 ? count_text(geometry_deck, c->keyword) : -1,
            c->reason);
    failures = 1;
  }
  capn_free(&arena);
  return failures;
}

static void test_render_rejects_system_scalar_directive(void **state) {
  (void)state;
  const struct SystemDirectiveCase cases[] = {
      {.name = "typed CUTOFF, cutoff 50",
       .has_sys = 1, .cutoff_ry = 50.0, .mult = 1,
       .directive = "CUTOFF 90", .keyword = "CUTOFF",
       .reason = "a CUTOFF directive in the typed &SYSTEM section"},
      {.name = "typed CUTOFF, cutoff 70",
       .has_sys = 1, .cutoff_ry = 70.0, .mult = 1,
       .directive = "CUTOFF 90", .keyword = "CUTOFF",
       .reason = "a CUTOFF directive in the typed &SYSTEM section"},
      {.name = "SYSTEM.CUTOFF, cutoff 50",
       .has_sys = 1, .cutoff_ry = 50.0, .mult = 1,
       .set_key = "SYSTEM.CUTOFF", .set_value = "90", .keyword = "CUTOFF",
       .reason = "a SYSTEM.CUTOFF set directive"},
      {.name = "SYSTEM.CUTOFF, cutoff 70",
       .has_sys = 1, .cutoff_ry = 70.0, .mult = 1,
       .set_key = "SYSTEM.CUTOFF", .set_value = "90", .keyword = "CUTOFF",
       .reason = "a SYSTEM.CUTOFF set directive"},
      {.name = "SYSTEM.CUTOFF, no typed section",
       .set_key = "SYSTEM.CUTOFF", .set_value = "90", .keyword = "CUTOFF",
       .reason = "a SYSTEM.CUTOFF set directive"},
      {.name = "typed CHARGE, charge 2",
       .has_sys = 1, .cutoff_ry = 70.0, .charge = 2, .mult = 1,
       .directive = "CHARGE 4", .keyword = "CHARGE",
       .reason = "a CHARGE directive in the typed &SYSTEM section"},
      {.name = "typed CHARGE, charge 0",
       .has_sys = 1, .cutoff_ry = 70.0, .charge = 0, .mult = 1,
       .directive = "CHARGE 4", .keyword = "CHARGE",
       .reason = "a CHARGE directive in the typed &SYSTEM section"},
      {.name = "SYSTEM.CHARGE, charge 2",
       .has_sys = 1, .cutoff_ry = 70.0, .charge = 2, .mult = 1,
       .set_key = "SYSTEM.CHARGE", .set_value = "4", .keyword = "CHARGE",
       .reason = "a SYSTEM.CHARGE set directive"},
      {.name = "SYSTEM.CHARGE, charge 0",
       .has_sys = 1, .cutoff_ry = 70.0, .charge = 0, .mult = 1,
       .set_key = "SYSTEM.CHARGE", .set_value = "4", .keyword = "CHARGE",
       .reason = "a SYSTEM.CHARGE set directive"},
      {.name = "SYSTEM.CHARGE, no typed section",
       .set_key = "SYSTEM.CHARGE", .set_value = "4", .keyword = "CHARGE",
       .reason = "a SYSTEM.CHARGE set directive"},
      {.name = "typed MULTIPLICITY, multiplicity 3",
       .has_sys = 1, .cutoff_ry = 70.0, .mult = 3,
       .directive = "MULTIPLICITY 5", .keyword = "MULTIPLICITY",
       .reason = "a MULTIPLICITY directive in the typed &SYSTEM section"},
      {.name = "typed MULTIPLICITY, multiplicity 1",
       .has_sys = 1, .cutoff_ry = 70.0, .mult = 1,
       .directive = "MULTIPLICITY 5", .keyword = "MULTIPLICITY",
       .reason = "a MULTIPLICITY directive in the typed &SYSTEM section"},
      {.name = "SYSTEM.MULTIPLICITY, multiplicity 3",
       .has_sys = 1, .cutoff_ry = 70.0, .mult = 3,
       .set_key = "SYSTEM.MULTIPLICITY", .set_value = "5",
       .keyword = "MULTIPLICITY",
       .reason = "a SYSTEM.MULTIPLICITY set directive"},
      {.name = "SYSTEM.MULTIPLICITY, multiplicity 1",
       .has_sys = 1, .cutoff_ry = 70.0, .mult = 1,
       .set_key = "SYSTEM.MULTIPLICITY", .set_value = "5",
       .keyword = "MULTIPLICITY",
       .reason = "a SYSTEM.MULTIPLICITY set directive"},
      {.name = "SYSTEM.MULTIPLICITY, no typed section",
       .set_key = "SYSTEM.MULTIPLICITY", .set_value = "5",
       .keyword = "MULTIPLICITY",
       .reason = "a SYSTEM.MULTIPLICITY set directive"},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += expect_system_directive_rejected(&cases[i]);
  assert_int_equal(failures, 0);
}

static void test_typed_scalars_report(void **state) {
  (void)state;
  struct capn arena;
  CPMDParams_ptr params;
  const struct Case defaults = {.has_dft = 1,
                                .dft_functional = "BLYP",
                                .has_sys = 1,
                                .sys_cutoff_ry = 70.0};
  build(&arena, &defaults, &params);
  CPMDCTypedScalars typed;
  assert_int_equal(cpmdc_params_typed_scalars(params, &typed), 0);
  assert_int_equal(typed.functional, 0);
  assert_int_equal(typed.cutoff, 0);
  assert_int_equal(typed.charge, 0);
  assert_int_equal(typed.multiplicity, 0);
  capn_free(&arena);

  const struct Case set = {.has_dft = 1,
                           .dft_functional = "LDA",
                           .has_sys = 1,
                           .sys_cutoff_ry = 50.0,
                           .sys_charge = 1,
                           .sys_mult = 2};
  build(&arena, &set, &params);
  assert_int_equal(cpmdc_params_typed_scalars(params, &typed), 0);
  assert_int_equal(typed.functional, 1);
  assert_int_equal(typed.cutoff, 1);
  assert_int_equal(typed.charge, 1);
  assert_int_equal(typed.multiplicity, 1);
  capn_free(&arena);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_overrides_fill_typed_defaults),
      cmocka_unit_test(test_typed_value_beats_override),
      cmocka_unit_test(test_top_level_value_with_default_sections),
      cmocka_unit_test(test_render_rejects_system_scalar_directive),
      cmocka_unit_test(test_typed_scalars_report),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
