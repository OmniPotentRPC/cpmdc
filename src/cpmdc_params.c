#include "cpmdc_params.h"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static const capn_text empty_text = {0, "", NULL};

static int pointer_list_len(capn_ptr *ptr) {
  capn_resolve(ptr);
  if (ptr->type == CAPN_NULL)
    return 0;
  if (ptr->type != CAPN_PTR_LIST)
    return -1;
  return ptr->len;
}

static int struct_list_len(capn_ptr *ptr) {
  capn_resolve(ptr);
  if (ptr->type == CAPN_NULL)
    return 0;
  if (ptr->type != CAPN_LIST)
    return -1;
  return ptr->len;
}

static int list64_len(capn_list64 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST || list->p.datasz != 8)
    return -1;
  return list->p.len;
}

static int list32_len(capn_list32 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST || list->p.datasz != 4)
    return -1;
  return list->p.len;
}

static int text_eq_ci(capn_text text, const char *lit) {
  size_t n = strlen(lit);
  if (text.len != (int)n || !text.str)
    return 0;
  for (size_t i = 0; i < n; ++i) {
    char a = (char)tolower((unsigned char)text.str[i]);
    char b = (char)tolower((unsigned char)lit[i]);
    if (a != b)
      return 0;
  }
  return 1;
}

static int text_starts_ci(capn_text text, const char *prefix) {
  size_t n = strlen(prefix);
  if (text.len < (int)n || !text.str)
    return 0;
  for (size_t i = 0; i < n; ++i) {
    char a = (char)tolower((unsigned char)text.str[i]);
    char b = (char)tolower((unsigned char)prefix[i]);
    if (a != b)
      return 0;
  }
  return 1;
}

static int append_fmt(char *dst, size_t dst_size, size_t *used, const char *fmt,
                      ...) {
  if (!dst || !used || *used >= dst_size)
    return -1;
  va_list ap;
  va_start(ap, fmt);
  int n = vsnprintf(dst + *used, dst_size - *used, fmt, ap);
  va_end(ap);
  if (n < 0 || (size_t)n >= dst_size - *used)
    return -1;
  *used += (size_t)n;
  return 0;
}

static int append_text(char *dst, size_t dst_size, size_t *used, const char *s) {
  return append_fmt(dst, dst_size, used, "%s", s ? s : "");
}

static int append_capn_text(char *dst, size_t dst_size, size_t *used,
                            capn_text text) {
  if (!text.str || text.len <= 0)
    return 0;
  if (*used + (size_t)text.len >= dst_size)
    return -1;
  memcpy(dst + *used, text.str, (size_t)text.len);
  *used += (size_t)text.len;
  dst[*used] = '\0';
  return 0;
}

static int append_f64_list_line(char *dst, size_t dst_size, size_t *used,
                                capn_list64 *list, int expected_len) {
  int n = list64_len(list);
  if (n < 0)
    return -1;
  if (n != expected_len)
    return n == 0 ? 0 : -1;
  if (append_text(dst, dst_size, used, "  ") != 0)
    return -1;
  for (int i = 0; i < n; ++i) {
    if (append_fmt(dst, dst_size, used, "%s%.10g", i == 0 ? "" : " ",
                   capn_to_f64(capn_get64(*list, i))) != 0)
      return -1;
  }
  return append_text(dst, dst_size, used, "\n");
}

static int append_i32_list_line(char *dst, size_t dst_size, size_t *used,
                                capn_list32 *list, int expected_len) {
  int n = list32_len(list);
  if (n < 0)
    return -1;
  if (n != expected_len)
    return n == 0 ? 0 : -1;
  if (append_text(dst, dst_size, used, "  ") != 0)
    return -1;
  for (int i = 0; i < n; ++i) {
    if (append_fmt(dst, dst_size, used, "%s%d", i == 0 ? "" : " ",
                   (int)(int32_t)capn_get32(*list, i)) != 0)
      return -1;
  }
  return append_text(dst, dst_size, used, "\n");
}

static int append_slice(char *dst, size_t dst_size, size_t *used,
                        const char *s, size_t len) {
  if (!s || len == 0)
    return 0;
  if (*used + len >= dst_size)
    return -1;
  memcpy(dst + *used, s, len);
  *used += len;
  dst[*used] = '\0';
  return 0;
}

static int append_directives(char *dst, size_t dst_size, size_t *used,
                             CPMDDirective_list directives) {
  int n = struct_list_len(&directives.p);
  if (n < 0)
    return -1;
  for (int i = 0; i < n; ++i) {
    struct CPMDDirective d;
    get_CPMDDirective(&d, directives, i);
    if (append_text(dst, dst_size, used, " ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, d.keyword) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    int na = pointer_list_len(&d.args);
    if (na < 0)
      return -1;
    for (int j = 0; j < na; ++j) {
      capn_text arg = capn_get_text(d.args, j, empty_text);
      if (append_text(dst, dst_size, used, "  ") != 0)
        return -1;
      if (append_capn_text(dst, dst_size, used, arg) != 0)
        return -1;
      if (append_text(dst, dst_size, used, "\n") != 0)
        return -1;
    }
  }
  return 0;
}

#define CPMDC_MAX_SET_DIRECTIVES 128

/** Parsed `SECTION.KEYWORD` directive waiting to merge into rendered output. */
struct RenderSetDirective {
  /** Section name slice from the set key. */
  const char *section;
  /** Byte length of `section`. */
  size_t section_len;
  /** Keyword name slice from the set key. */
  const char *keyword;
  /** Byte length of `keyword`. */
  size_t keyword_len;
  /** Optional directive value emitted on the following deck line. */
  capn_text value;
  /** Non-zero once the directive has been emitted. */
  int rendered;
};

/** Bounded collection of parsed set directives. */
struct RenderSetList {
  /** Parsed directives in input order. */
  struct RenderSetDirective items[CPMDC_MAX_SET_DIRECTIVES];
  /** Number of valid entries in `items`. */
  int len;
};

static int ascii_upper(int c) {
  if (c >= 'a' && c <= 'z')
    return c - 'a' + 'A';
  return c;
}

static int section_name_equals(const struct RenderSetDirective *set,
                               const char *name) {
  size_t name_len = strlen(name);
  if (set->section_len != name_len)
    return 0;
  for (size_t i = 0; i < name_len; ++i) {
    if (ascii_upper((unsigned char)set->section[i]) !=
        ascii_upper((unsigned char)name[i]))
      return 0;
  }
  return 1;
}

static int set_sections_equal(const struct RenderSetDirective *a,
                              const struct RenderSetDirective *b) {
  if (a->section_len != b->section_len)
    return 0;
  for (size_t i = 0; i < a->section_len; ++i) {
    if (ascii_upper((unsigned char)a->section[i]) !=
        ascii_upper((unsigned char)b->section[i]))
      return 0;
  }
  return 1;
}

static int slice_starts_ci(const char *text, size_t text_len,
                           const char *prefix) {
  size_t prefix_len = strlen(prefix);
  if (!text || text_len < prefix_len)
    return 0;
  for (size_t i = 0; i < prefix_len; ++i) {
    if (ascii_upper((unsigned char)text[i]) !=
        ascii_upper((unsigned char)prefix[i]))
      return 0;
  }
  return 1;
}

static int init_render_set(struct RenderSetDirective *out,
                           const struct CPMDSetDirective *set) {
  if (!set->key.str || set->key.len <= 2)
    return -1;
  const char *dot = memchr(set->key.str, '.', (size_t)set->key.len);
  if (!dot || dot == set->key.str ||
      dot == set->key.str + (size_t)set->key.len - 1)
    return -1;
  out->section = set->key.str;
  out->section_len = (size_t)(dot - set->key.str);
  out->keyword = dot + 1;
  out->keyword_len = (size_t)set->key.len - out->section_len - 1u;
  out->value = set->value;
  out->rendered = 0;
  return 0;
}

static int collect_set_directives(CPMDInputSection_list sections,
                                  struct RenderSetList *sets) {
  int nsec = struct_list_len(&sections.p);
  if (nsec < 0)
    return -1;
  sets->len = 0;
  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, sections, i);
    if (sec.which != CPMDInputSection_set)
      continue;
    if (sets->len >= CPMDC_MAX_SET_DIRECTIVES)
      return -1;
    struct CPMDSetDirective body;
    read_CPMDSetDirective(&body, sec.set);
    if (init_render_set(&sets->items[sets->len], &body) != 0)
      return -1;
    ++sets->len;
  }
  return 0;
}

static int set_directives_have_prefix(struct RenderSetList *sets,
                                      const char *section,
                                      const char *prefix) {
  if (!sets)
    return 0;
  for (int i = 0; i < sets->len; ++i) {
    if (sets->items[i].rendered ||
        !section_name_equals(&sets->items[i], section))
      continue;
    if (slice_starts_ci(sets->items[i].keyword, sets->items[i].keyword_len,
                        prefix))
      return 1;
  }
  return 0;
}

static int append_set_payload(char *dst, size_t dst_size, size_t *used,
                              const struct RenderSetDirective *set) {
  if (append_text(dst, dst_size, used, " ") != 0)
    return -1;
  if (append_slice(dst, dst_size, used, set->keyword, set->keyword_len) != 0)
    return -1;
  if (append_text(dst, dst_size, used, "\n") != 0)
    return -1;
  if (set->value.str && set->value.len > 0) {
    if (append_text(dst, dst_size, used, "  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, set->value) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  return 0;
}

static int append_set_directives_for_section(char *dst, size_t dst_size,
                                             size_t *used,
                                             struct RenderSetList *sets,
                                             const char *section) {
  if (!sets)
    return 0;
  for (int i = 0; i < sets->len; ++i) {
    if (sets->items[i].rendered || !section_name_equals(&sets->items[i], section))
      continue;
    if (append_set_payload(dst, dst_size, used, &sets->items[i]) != 0)
      return -1;
    sets->items[i].rendered = 1;
  }
  return 0;
}

static capn_text selected_file_path(const struct CPMDParams *view) {
  if (view->permanentDir.str && view->permanentDir.len > 0)
    return view->permanentDir;
  if (view->scratchDir.str && view->scratchDir.len > 0)
    return view->scratchDir;
  return empty_text;
}

static int append_file_path_directive(char *dst, size_t dst_size, size_t *used,
                                      capn_text path) {
  if (!path.str || path.len <= 0)
    return 0;
  if (append_text(dst, dst_size, used, " FILEPATH\n  ") != 0)
    return -1;
  if (append_capn_text(dst, dst_size, used, path) != 0)
    return -1;
  if (path.str[path.len - 1] != '/') {
    if (append_text(dst, dst_size, used, "/") != 0)
      return -1;
  }
  return append_text(dst, dst_size, used, "\n");
}

static int append_text_keyword_arg(char *dst, size_t dst_size, size_t *used,
                                   const char *keyword, capn_text value) {
  if (!value.str || value.len <= 0)
    return 0;
  if (append_text(dst, dst_size, used, " ") != 0)
    return -1;
  if (append_text(dst, dst_size, used, keyword) != 0)
    return -1;
  if (append_text(dst, dst_size, used, " ") != 0)
    return -1;
  if (append_capn_text(dst, dst_size, used, value) != 0)
    return -1;
  return append_text(dst, dst_size, used, "\n");
}

static int render_remaining_set_sections(char *dst, size_t dst_size,
                                         size_t *used,
                                         struct RenderSetList *sets) {
  if (!sets)
    return 0;
  for (int i = 0; i < sets->len; ++i) {
    if (sets->items[i].rendered)
      continue;
    if (append_text(dst, dst_size, used, "&") != 0)
      return -1;
    if (append_slice(dst, dst_size, used, sets->items[i].section,
                     sets->items[i].section_len) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    for (int j = i; j < sets->len; ++j) {
      if (sets->items[j].rendered ||
          !set_sections_equal(&sets->items[i], &sets->items[j]))
        continue;
      if (append_set_payload(dst, dst_size, used, &sets->items[j]) != 0)
        return -1;
      sets->items[j].rendered = 1;
    }
    if (append_text(dst, dst_size, used, "&END\n\n") != 0)
      return -1;
  }
  return 0;
}

static int directives_have_prefix(CPMDDirective_list directives,
                                  const char *prefix) {
  int n = struct_list_len(&directives.p);
  if (n < 0)
    return -1;
  for (int i = 0; i < n; ++i) {
    struct CPMDDirective d;
    get_CPMDDirective(&d, directives, i);
    if (text_starts_ci(d.keyword, prefix))
      return 1;
  }
  return 0;
}

static int append_cpmd_cell(char *dst, size_t dst_size, size_t *used,
                            const double *cell, int ncell) {
  if (ncell != 6 && ncell != 9)
    return 0;
  if (append_text(dst, dst_size, used, " CELL\n ") != 0)
    return -1;
  if (ncell == 6) {
    double a = cell[0];
    double b = cell[1];
    double c = cell[2];
    double ab = cell[3];
    double ac = cell[4];
    double bc = cell[5];
    if (fabs(ab) > 1.0 || fabs(ac) > 1.0 || fabs(bc) > 1.0) {
      const double deg = 3.14159265358979323846 / 180.0;
      if (a != 0.0) {
        b /= a;
        c /= a;
      }
      ab = cos(ab * deg);
      ac = cos(ac * deg);
      bc = cos(bc * deg);
    }
    return append_fmt(dst, dst_size, used, " %.10g %.10g %.10g %.10g %.10g %.10g\n",
                      a, b, c, ab, ac, bc);
  }
  if (fabs(cell[1]) < 1e-12 && fabs(cell[2]) < 1e-12 &&
      fabs(cell[3]) < 1e-12 && fabs(cell[5]) < 1e-12 &&
      fabs(cell[6]) < 1e-12 && fabs(cell[7]) < 1e-12 && cell[0] != 0.0) {
    return append_fmt(dst, dst_size, used, " %.10g %.10g %.10g 0 0 0\n",
                      cell[0], cell[4] / cell[0], cell[8] / cell[0]);
  }
  for (int i = 0; i < ncell; ++i) {
    if (append_fmt(dst, dst_size, used, " %.10g", cell[i]) != 0)
      return -1;
  }
  return append_text(dst, dst_size, used, "\n");
}

static int render_system_section_with_cell(
    char *dst, size_t dst_size, size_t *used, struct CPMDSystemSection *sys,
    double default_cutoff, int default_charge, const double *cell_override,
    int override_ncell, struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&SYSTEM\n") != 0)
    return -1;
  int symmetry = sys->symmetry;
  int charge = sys->charge != 0 ? sys->charge : default_charge;
  double cutoff = sys->cutOffRy > 0.0 ? sys->cutOffRy : default_cutoff;
  if (append_fmt(dst, dst_size, used, " SYMMETRY\n  %d\n", symmetry) != 0)
    return -1;
  if (sys->angstrom) {
    if (append_text(dst, dst_size, used, " ANGSTROM\n") != 0)
      return -1;
  }
  int ncell = list64_len(&sys->cell);
  if (ncell < 0)
    return -1;
  if (override_ncell == 9 && cell_override) {
    if (append_cpmd_cell(dst, dst_size, used, cell_override, override_ncell) != 0)
      return -1;
  } else if (ncell == 6 || ncell == 9) {
    double cell[9] = {0};
    for (int i = 0; i < ncell; ++i)
      cell[i] = capn_to_f64(capn_get64(sys->cell, i));
    if (append_cpmd_cell(dst, dst_size, used, cell, ncell) != 0)
      return -1;
  }
  int n_reference_cell = list64_len(&sys->referenceCell);
  if (n_reference_cell < 0)
    return -1;
  if (n_reference_cell > 0) {
    if (append_text(dst, dst_size, used,
                    n_reference_cell == 9 ? " REFERENCE CELL VECTORS\n"
                                          : " REFERENCE CELL\n") != 0)
      return -1;
    if (append_f64_list_line(dst, dst_size, used, &sys->referenceCell,
                             n_reference_cell == 9 ? 9 : 6) != 0)
      return -1;
  }
  int n_classical_cell = list64_len(&sys->classicalCell);
  if (n_classical_cell < 0)
    return -1;
  if (n_classical_cell > 0) {
    if (append_text(dst, dst_size, used, " CLASSICAL CELL\n") != 0)
      return -1;
    if (append_f64_list_line(dst, dst_size, used, &sys->classicalCell, 6) != 0)
      return -1;
  }
  if (sys->isotropicCell) {
    if (append_text(dst, dst_size, used, " ISOTROPIC CELL\n") != 0)
      return -1;
  }
  if (sys->zFlexibleCell) {
    if (append_text(dst, dst_size, used, " ZFLEXIBLE CELL\n") != 0)
      return -1;
  }
  if (append_fmt(dst, dst_size, used, " CUTOFF\n  %.10g\n", cutoff) != 0)
    return -1;
  if (sys->densityCutOffRy > 0.0) {
    if (append_fmt(dst, dst_size, used, " DENSITY CUTOFF\n  %.10g\n",
                   sys->densityCutOffRy) != 0)
      return -1;
  }
  if (sys->densityCutoffNumber > 0) {
    if (append_fmt(dst, dst_size, used, " DENSITY CUTOFF NUMBER\n  %d\n",
                   sys->densityCutoffNumber) != 0)
      return -1;
  }
  if (sys->dual > 0.0) {
    if (append_fmt(dst, dst_size, used, " DUAL\n  %.10g\n", sys->dual) != 0)
      return -1;
  }
  int n_constant_cutoff = list64_len(&sys->constantCutoff);
  if (n_constant_cutoff < 0)
    return -1;
  if (n_constant_cutoff > 0) {
    if (append_text(dst, dst_size, used, " CONSTANT CUTOFF\n") != 0)
      return -1;
    if (append_f64_list_line(dst, dst_size, used, &sys->constantCutoff, 3) !=
        0)
      return -1;
  }
  if (sys->scale != 0.0 || sys->scaleCartesian) {
    if (append_text(dst, dst_size, used, " SCALE") != 0)
      return -1;
    if (sys->scaleCartesian &&
        append_text(dst, dst_size, used, " CARTESIAN") != 0)
      return -1;
    if (sys->scale != 0.0 &&
        append_fmt(dst, dst_size, used, " S=%.10g", sys->scale) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (charge != 0) {
    if (append_fmt(dst, dst_size, used, " CHARGE\n  %d\n", charge) != 0)
      return -1;
  }
  int field_has_poisson =
      sys->poissonSolver.str && sys->poissonSolver.len > 0;
  if (field_has_poisson) {
    if (append_text(dst, dst_size, used, " POISSON SOLVER ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sys->poissonSolver) != 0)
      return -1;
    if (sys->poissonParameter > 0.0) {
      if (append_text(dst, dst_size, used, " PARAMETER\n  ") != 0)
        return -1;
      if (append_fmt(dst, dst_size, used, "%.10g", sys->poissonParameter) != 0)
        return -1;
    }
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  int n_mesh = list32_len(&sys->mesh);
  if (n_mesh < 0)
    return -1;
  if (n_mesh > 0) {
    if (append_text(dst, dst_size, used, " MESH\n") != 0)
      return -1;
    if (append_i32_list_line(dst, dst_size, used, &sys->mesh, 3) != 0)
      return -1;
  }
  if (sys->doubleGrid.str && sys->doubleGrid.len > 0) {
    if (append_text(dst, dst_size, used, " DOUBLE GRID ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sys->doubleGrid) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sys->symmetrizeCoordinates) {
    if (append_text(dst, dst_size, used, " SYMMETRIZE COORDINATES\n") != 0)
      return -1;
  }
  if (sys->tesr > 0) {
    if (append_fmt(dst, dst_size, used, " TESR\n  %d\n", sys->tesr) != 0)
      return -1;
  }
  if (sys->surface.str && sys->surface.len > 0) {
    if (append_text(dst, dst_size, used, " SURFACE ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sys->surface) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sys->polymer) {
    if (append_text(dst, dst_size, used, " POLYMER\n") != 0)
      return -1;
  }
  if (sys->cluster) {
    if (append_text(dst, dst_size, used, " CLUSTER\n") != 0)
      return -1;
  }
  int has_poisson = directives_have_prefix(sys->directives, "POISSON SOLVER");
  if (has_poisson < 0)
    return -1;
  int set_has_poisson =
      set_directives_have_prefix(sets, "SYSTEM", "POISSON SOLVER");
  if (set_has_poisson < 0)
    return -1;
  if (symmetry == 0 && !field_has_poisson && !has_poisson &&
      !set_has_poisson) {
    if (append_text(dst, dst_size, used, " POISSON SOLVER HOCKNEY\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, sys->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "SYSTEM") !=
      0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_system_section(char *dst, size_t dst_size, size_t *used,
                                 struct CPMDSystemSection *sys,
                                 double default_cutoff, int default_charge,
                                 struct RenderSetList *sets) {
  return render_system_section_with_cell(dst, dst_size, used, sys,
                                         default_cutoff, default_charge, NULL, 0,
                                         sets);
}

static int render_cpmd_section(char *dst, size_t dst_size, size_t *used,
                               const struct CPMDCpmdSection *sec,
                               const char *task, capn_text file_path,
                               struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&CPMD\n") != 0)
    return -1;
  int opt_wf = sec->optimizeWavefunction;
  int do_md = sec->molecularDynamics;
  if (task && (strcmp(task, "md") == 0 || strcmp(task, "MD") == 0))
    do_md = 1;
  if (task && (strcmp(task, "energy") == 0 || strcmp(task, "gradient") == 0 ||
               strcmp(task, "optimize") == 0))
    opt_wf = 1;
  if (opt_wf && !do_md) {
    if (append_text(dst, dst_size, used, " OPTIMIZE WAVEFUNCTION\n") != 0)
      return -1;
  }
  if (sec->optimizeGeometry) {
    if (append_text(dst, dst_size, used, " OPTIMIZE GEOMETRY\n") != 0)
      return -1;
  }
  if (do_md) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS\n") != 0)
      return -1;
  }
  if (sec->isolatedMolecule) {
    if (append_text(dst, dst_size, used, " ISOLATED MOLECULE\n") != 0)
      return -1;
  }
  if (sec->molecularDynamicsCp) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS CP\n") != 0)
      return -1;
  }
  if (sec->molecularDynamicsBo) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS BO\n") != 0)
      return -1;
  }
  if (sec->molecularDynamicsEh) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS EH\n") != 0)
      return -1;
  }
  if (sec->molecularDynamicsPt) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS PT\n") != 0)
      return -1;
  }
  if (sec->molecularDynamicsClassical) {
    if (append_text(dst, dst_size, used,
                    " MOLECULAR DYNAMICS CLASSICAL\n") != 0)
      return -1;
  }
  if (sec->molecularDynamicsFile.str && sec->molecularDynamicsFile.len > 0) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS FILE\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->molecularDynamicsFile) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->convergenceOrbitals > 0.0) {
    if (append_fmt(dst, dst_size, used, " CONVERGENCE ORBITALS\n  %.10g\n",
                   sec->convergenceOrbitals) != 0)
      return -1;
  }
  if (sec->maxStep > 0) {
    if (append_fmt(dst, dst_size, used, " MAXSTEP\n  %d\n", sec->maxStep) != 0)
      return -1;
  }
  if (sec->maxIter > 0) {
    if (append_fmt(dst, dst_size, used, " MAXITER\n  %d\n", sec->maxIter) != 0)
      return -1;
  }
  if (sec->convergenceGeometry > 0.0) {
    if (append_fmt(dst, dst_size, used, " CONVERGENCE GEOMETRY\n  %.10g\n",
                   sec->convergenceGeometry) != 0)
      return -1;
  }
  if (sec->timestep > 0.0) {
    if (append_fmt(dst, dst_size, used, " TIMESTEP\n  %.10g\n", sec->timestep) !=
        0)
      return -1;
  }
  if (sec->electronMass > 0.0) {
    if (append_fmt(dst, dst_size, used, " EMASS\n  %.10g\n",
                   sec->electronMass) != 0)
      return -1;
  }
  if (sec->nose) {
    if (append_text(dst, dst_size, used, " NOSE\n") != 0)
      return -1;
  }
  if (sec->noseIons) {
    if (append_text(dst, dst_size, used, " NOSE IONS\n") != 0)
      return -1;
  }
  if (sec->noseElectrons) {
    if (append_text(dst, dst_size, used, " NOSE ELECTRONS\n") != 0)
      return -1;
  }
  if (sec->berendsen.str && sec->berendsen.len > 0) {
    if (append_text(dst, dst_size, used, " BERENDSEN\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->berendsen) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->langevin) {
    if (append_text(dst, dst_size, used, " LANGEVIN\n") != 0)
      return -1;
  }
  if (sec->annealing.str && sec->annealing.len > 0) {
    if (append_text(dst, dst_size, used, " ANNEALING\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->annealing) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->quench) {
    if (append_text(dst, dst_size, used, " QUENCH\n") != 0)
      return -1;
  }
  if (sec->rattle) {
    if (append_text(dst, dst_size, used, " RATTLE\n") != 0)
      return -1;
  }
  if (sec->shake) {
    if (append_text(dst, dst_size, used, " SHAKE\n") != 0)
      return -1;
  }
  if (sec->constraint.str && sec->constraint.len > 0) {
    if (append_text(dst, dst_size, used, " CONSTRAINT\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->constraint) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->trotter.str && sec->trotter.len > 0) {
    if (append_text(dst, dst_size, used, " TROTTER\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->trotter) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->restart) {
    if (append_text(dst, dst_size, used, " RESTART\n") != 0)
      return -1;
  }
  if (sec->printOptions.str && sec->printOptions.len > 0) {
    if (append_text(dst, dst_size, used, " PRINT\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->printOptions) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->storeOptions.str && sec->storeOptions.len > 0) {
    if (append_text(dst, dst_size, used, " STORE\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sec->storeOptions) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (sec->centerMoleculeOff) {
    if (append_text(dst, dst_size, used, " CENTER MOLECULE OFF\n") != 0)
      return -1;
  }
  if (sec->centerMoleculeOn) {
    if (append_text(dst, dst_size, used, " CENTER MOLECULE ON\n") != 0)
      return -1;
  }
  if (sec->diis) {
    if (append_text(dst, dst_size, used, " DIIS\n") != 0)
      return -1;
  }
  if (sec->odiis) {
    if (append_text(dst, dst_size, used, " ODIIS\n") != 0)
      return -1;
  }
  if (sec->pcg) {
    if (append_text(dst, dst_size, used, " PCG\n") != 0)
      return -1;
  }
  if (sec->diagonalization) {
    if (append_text(dst, dst_size, used, " DIAGONALIZATION\n") != 0)
      return -1;
  }
  if (sec->freeEnergy) {
    if (append_text(dst, dst_size, used, " FREE-ENERGY\n") != 0)
      return -1;
  }
  if (sec->_interface) {
    if (append_text(dst, dst_size, used, " INTERFACE\n") != 0)
      return -1;
  }
  if (sec->qmmm) {
    if (append_text(dst, dst_size, used, " QMMM\n") != 0)
      return -1;
  }
  if (sec->bicanonicalEnsemble) {
    if (append_text(dst, dst_size, used, " BICANONICAL ENSEMBLE\n") != 0)
      return -1;
  }
  if (sec->cdft) {
    if (append_text(dst, dst_size, used, " CDFT\n") != 0)
      return -1;
  }
  if (sec->properties) {
    if (append_text(dst, dst_size, used, " PROPERTIES\n") != 0)
      return -1;
  }
  if (append_text_keyword_arg(dst, dst_size, used, "VDW CORRECTION",
                              sec->vdwCorrection) != 0)
    return -1;
  if (append_text_keyword_arg(dst, dst_size, used, "VDW WANNIER",
                              sec->vdwWannier) != 0)
    return -1;
  if (sec->dcacp) {
    if (append_text(dst, dst_size, used, " DCACP\n") != 0)
      return -1;
  }
  if (sec->restartWavefunction) {
    if (append_text(dst, dst_size, used, " RESTART WAVEFUNCTION\n") != 0)
      return -1;
  }
  if (sec->trajectory) {
    if (append_text(dst, dst_size, used, " TRAJECTORY\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, sec->directives) != 0)
    return -1;
  int has_filepath = directives_have_prefix(sec->directives, "FILEPATH");
  if (has_filepath < 0)
    return -1;
  int set_has_filepath = set_directives_have_prefix(sets, "CPMD", "FILEPATH");
  if (set_has_filepath < 0)
    return -1;
  if (!has_filepath && !set_has_filepath) {
    if (append_file_path_directive(dst, dst_size, used, file_path) != 0)
      return -1;
  }
  if (append_set_directives_for_section(dst, dst_size, used, sets, "CPMD") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_dft_section(char *dst, size_t dst_size, size_t *used,
                              const struct CPMDDftSection *dft,
                              const char *default_functional, int multiplicity,
                              struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&DFT\n") != 0)
    return -1;
  if (dft->functional.str && dft->functional.len > 0) {
    if (append_text(dst, dst_size, used, " FUNCTIONAL ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->functional) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  } else {
    if (append_fmt(dst, dst_size, used, " FUNCTIONAL %s\n",
                   default_functional ? default_functional : "BLYP") != 0)
      return -1;
  }
  if (dft->lsd || multiplicity > 1) {
    if (append_text(dst, dst_size, used, " LSD\n") != 0)
      return -1;
  }
  if (dft->gcCutoff > 0.0) {
    if (append_fmt(dst, dst_size, used, " GC-CUTOFF\n  %.10g\n",
                   dft->gcCutoff) != 0)
      return -1;
  }
  if (dft->xcDriver.str && dft->xcDriver.len > 0) {
    if (append_text(dst, dst_size, used, " XC_DRIVER\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->xcDriver) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->libxc.str && dft->libxc.len > 0) {
    if (append_text(dst, dst_size, used, " LIBXC\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->libxc) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->lrKernel.str && dft->lrKernel.len > 0) {
    if (append_text(dst, dst_size, used, " LR KERNEL\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->lrKernel) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->refunct.str && dft->refunct.len > 0) {
    if (append_text(dst, dst_size, used, " REFUNCT\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->refunct) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->mtsHighFunc.str && dft->mtsHighFunc.len > 0) {
    if (append_text(dst, dst_size, used, " MTS_HIGH_FUNC\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->mtsHighFunc) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->mtsLowFunc.str && dft->mtsLowFunc.len > 0) {
    if (append_text(dst, dst_size, used, " MTS_LOW_FUNC\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->mtsLowFunc) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->hfx) {
    if (append_text(dst, dst_size, used, " HFX\n") != 0)
      return -1;
  }
  if (dft->hfxScreening.str && dft->hfxScreening.len > 0) {
    if (append_text(dst, dst_size, used, " HFX-SCREENING\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->hfxScreening) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->hubbard.str && dft->hubbard.len > 0) {
    if (append_text(dst, dst_size, used, " HUBBARD\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->hubbard) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->alpha > 0.0) {
    if (append_fmt(dst, dst_size, used, " ALPHA\n  %.10g\n", dft->alpha) != 0)
      return -1;
  }
  if (dft->beta > 0.0) {
    if (append_fmt(dst, dst_size, used, " BETA\n  %.10g\n", dft->beta) != 0)
      return -1;
  }
  if (dft->oldCode) {
    if (append_text(dst, dst_size, used, " OLDCODE\n") != 0)
      return -1;
  }
  if (dft->newCode) {
    if (append_text(dst, dst_size, used, " NEWCODE\n") != 0)
      return -1;
  }
  if (dft->correlation.str && dft->correlation.len > 0) {
    if (append_text(dst, dst_size, used, " CORRELATION\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->correlation) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->exchange.str && dft->exchange.len > 0) {
    if (append_text(dst, dst_size, used, " EXCHANGE\n  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->exchange) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  if (dft->becke88) {
    if (append_text(dst, dst_size, used, " BECKE88\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, dft->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "DFT") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

/* CPMD &ATOMS: *filename.psp then LMAX=S|P|D and coordinates (Angstrom). */
static const char *lmax_letter(int lmax) {
  if (lmax <= 0)
    return "S";
  if (lmax == 1)
    return "P";
  if (lmax == 2)
    return "D";
  return "F";
}

static int render_atoms_section(char *dst, size_t dst_size, size_t *used,
                                struct CPMDAtomsSection *atoms,
                                struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&ATOMS\n") != 0)
    return -1;
  int npsp = struct_list_len(&atoms->pseudopotentials.p);
  if (npsp < 0)
    return -1;
  for (int i = 0; i < npsp; ++i) {
    struct CPMDAtomsPseudopotential psp;
    get_CPMDAtomsPseudopotential(&psp, atoms->pseudopotentials, i);
    /* path is the PP filename (e.g. O_MT_BLYP.psp); element is Z symbol only. */
    if (append_text(dst, dst_size, used, "*") != 0)
      return -1;
    if (psp.path.str && psp.path.len > 0) {
      if (append_capn_text(dst, dst_size, used, psp.path) != 0)
        return -1;
    } else if (psp.element.str && psp.element.len > 0) {
      if (append_capn_text(dst, dst_size, used, psp.element) != 0)
        return -1;
      if (append_text(dst, dst_size, used, ".psp") != 0)
        return -1;
    } else {
      return -1;
    }
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    if (append_fmt(dst, dst_size, used, " LMAX=%s\n",
                   lmax_letter(psp.lmax >= 0 ? psp.lmax : 0)) != 0)
      return -1;
    /* Coordinate counts filled by geometry merge (see render_deck_with_geometry). */
    if (append_text(dst, dst_size, used, "   0\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, atoms->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "ATOMS") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_generic_section(char *dst, size_t dst_size, size_t *used,
                                  const struct CPMDGenericSection *sec,
                                  struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&") != 0)
    return -1;
  if (append_capn_text(dst, dst_size, used, sec->name) != 0)
    return -1;
  if (append_text(dst, dst_size, used, "\n") != 0)
    return -1;
  if (append_directives(dst, dst_size, used, sec->directives) != 0)
    return -1;
  if (sec->name.str && sec->name.len > 0) {
    char section[64];
    size_t n = (size_t)sec->name.len;
    if (n >= sizeof(section))
      return -1;
    memcpy(section, sec->name.str, n);
    section[n] = '\0';
    if (append_set_directives_for_section(dst, dst_size, used, sets, section) !=
        0)
      return -1;
  }
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int append_subsections(char *dst, size_t dst_size, size_t *used,
                              CPMDGenericSection_list subsections) {
  int nsub = struct_list_len(&subsections.p);
  if (nsub < 0)
    return -1;
  for (int i = 0; i < nsub; ++i) {
    struct CPMDGenericSection sub;
    get_CPMDGenericSection(&sub, subsections, i);
    if (!sub.name.str || sub.name.len <= 0)
      return -1;
    if (append_text(dst, dst_size, used, " ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sub.name) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    if (append_directives(dst, dst_size, used, sub.directives) != 0)
      return -1;
    if (append_text(dst, dst_size, used, " END ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, sub.name) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  return 0;
}

static int render_directive_section(char *dst, size_t dst_size, size_t *used,
                                    const char *section,
                                    const struct CPMDDirectiveSection *sec,
                                    struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&") != 0)
    return -1;
  if (append_text(dst, dst_size, used, section) != 0)
    return -1;
  if (append_text(dst, dst_size, used, "\n") != 0)
    return -1;
  if (append_directives(dst, dst_size, used, sec->directives) != 0)
    return -1;
  if (append_subsections(dst, dst_size, used, sec->subsections) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, section) != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_typed_directive_section(char *dst, size_t dst_size,
                                          size_t *used,
                                          const struct CPMDInputSection *sec,
                                          struct RenderSetList *sets) {
#define CPMDC_RENDER_DIRECTIVE_CASE(kind, section_name)                         \
  case CPMDInputSection_##kind: {                                               \
    struct CPMDDirectiveSection body;                                           \
    read_CPMDDirectiveSection(&body, sec->kind);                                \
    return render_directive_section(dst, dst_size, used, section_name, &body,   \
                                    sets);                                      \
  }
  switch (sec->which) {
    CPMDC_RENDER_DIRECTIVE_CASE(atom, "ATOM")
    CPMDC_RENDER_DIRECTIVE_CASE(basis, "BASIS")
    CPMDC_RENDER_DIRECTIVE_CASE(clas, "CLAS")
    CPMDC_RENDER_DIRECTIVE_CASE(eam, "EAM")
    CPMDC_RENDER_DIRECTIVE_CASE(exte, "EXTE")
    CPMDC_RENDER_DIRECTIVE_CASE(hardness, "HARDNESS")
    CPMDC_RENDER_DIRECTIVE_CASE(info, "INFO")
    CPMDC_RENDER_DIRECTIVE_CASE(linres, "LINRES")
    CPMDC_RENDER_DIRECTIVE_CASE(molstates, "MOLSTATES")
    CPMDC_RENDER_DIRECTIVE_CASE(mts, "MTS")
    CPMDC_RENDER_DIRECTIVE_CASE(nlcc, "NLCC")
    CPMDC_RENDER_DIRECTIVE_CASE(path, "PATH")
    CPMDC_RENDER_DIRECTIVE_CASE(pimd, "PIMD")
    CPMDC_RENDER_DIRECTIVE_CASE(potential, "POTENTIAL")
    CPMDC_RENDER_DIRECTIVE_CASE(prop, "PROP")
    CPMDC_RENDER_DIRECTIVE_CASE(ptddft, "PTDDFT")
    CPMDC_RENDER_DIRECTIVE_CASE(resp, "RESP")
    CPMDC_RENDER_DIRECTIVE_CASE(tddft, "TDDFT")
    CPMDC_RENDER_DIRECTIVE_CASE(vdw, "VDW")
    CPMDC_RENDER_DIRECTIVE_CASE(vectors, "VECTORS")
    CPMDC_RENDER_DIRECTIVE_CASE(wavefunction, "WAVEFUNCTION")
  default:
    return -1;
  }
#undef CPMDC_RENDER_DIRECTIVE_CASE
}

int cpmdc_params_root(const void *params_capnp, size_t params_capnp_size_bytes,
                      struct capn *arena, CPMDParams_ptr *params) {
  if (!params_capnp || params_capnp_size_bytes == 0 || !arena || !params)
    return -1;
  memset(arena, 0, sizeof(*arena));
  memset(params, 0, sizeof(*params));
  if (capn_init_mem(arena, (const uint8_t *)params_capnp,
                    params_capnp_size_bytes, 0) != 0)
    return -1;
  params->p = capn_getp(capn_root(arena), 0, 1);
  if (params->p.type != CAPN_STRUCT) {
    cpmdc_params_release(arena);
    memset(params, 0, sizeof(*params));
    return -1;
  }
  return 0;
}

void cpmdc_params_release(struct capn *arena) {
  if (!arena)
    return;
  capn_free(arena);
  memset(arena, 0, sizeof(*arena));
}

const char *cpmdc_params_text_or(capn_text text, const char *fallback) {
  if (text.str && text.len > 0)
    return text.str;
  return fallback ? fallback : "";
}

static void copy_capn_text_or(capn_text text, const char *fallback, char *dst,
                              size_t dst_size) {
  if (!dst || dst_size == 0)
    return;
  dst[0] = '\0';
  if (text.str && text.len > 0) {
    size_t n = (size_t)text.len;
    if (n >= dst_size)
      n = dst_size - 1u;
    memcpy(dst, text.str, n);
    dst[n] = '\0';
    return;
  }
  snprintf(dst, dst_size, "%s", fallback ? fallback : "");
}

int cpmdc_params_effective_config(CPMDParams_ptr params, char *functional,
                                  size_t functional_size, double *cutoff_ry,
                                  int *charge, int *multiplicity) {
  if (params.p.type == CAPN_NULL || !functional || functional_size == 0 ||
      !cutoff_ry || !charge || !multiplicity)
    return -1;

  struct CPMDParams view;
  read_CPMDParams(&view, params);
  copy_capn_text_or(view.functional, "BLYP", functional, functional_size);
  *cutoff_ry = view.cutOffRy > 0.0 ? view.cutOffRy : 70.0;
  *charge = view.charge;
  *multiplicity = view.multiplicity > 0 ? view.multiplicity : 1;

  int nsec = struct_list_len(&view.inputSections.p);
  if (nsec < 0)
    return -1;
  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, view.inputSections, i);
    if (sec.which == CPMDInputSection_system) {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      if (body.cutOffRy > 0.0)
        *cutoff_ry = body.cutOffRy;
      *charge = body.charge;
      if (body.multiplicity > 0)
        *multiplicity = body.multiplicity;
    } else if (sec.which == CPMDInputSection_dft) {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      if (body.functional.str && body.functional.len > 0)
        copy_capn_text_or(body.functional, "BLYP", functional,
                          functional_size);
    }
  }
  return 0;
}

int cpmdc_params_render_input_deck(CPMDParams_ptr params, char *dst,
                                   size_t dst_size) {
  if (params.p.type == CAPN_NULL || !dst || dst_size == 0)
    return -1;
  dst[0] = '\0';
  size_t used = 0;

  struct CPMDParams view;
  read_CPMDParams(&view, params);

  if (view.title.str && view.title.len > 0) {
    if (append_text(dst, dst_size, &used, "! ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, &used, view.title) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  int nblocks = pointer_list_len(&view.inputBlocks);
  if (nblocks < 0)
    return -1;
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(view.inputBlocks, i, empty_text);
    if (append_capn_text(dst, dst_size, &used, block) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  const char *task = cpmdc_params_text_or(view.task, "gradient");
  const char *functional = cpmdc_params_text_or(view.functional, "BLYP");
  capn_text file_path = selected_file_path(&view);
  double cutoff = view.cutOffRy > 0.0 ? view.cutOffRy : 70.0;
  int charge = view.charge;
  int mult = view.multiplicity > 0 ? view.multiplicity : 1;

  int has_system = 0, has_cpmd = 0, has_dft = 0, has_atoms = 0;
  int nsec = struct_list_len(&view.inputSections.p);
  if (nsec < 0)
    return -1;
  struct RenderSetList sets;
  if (collect_set_directives(view.inputSections, &sets) != 0)
    return -1;

  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, view.inputSections, i);
    switch (sec.which) {
    case CPMDInputSection_generic: {
      struct CPMDGenericSection body;
      read_CPMDGenericSection(&body, sec.generic);
      if (render_generic_section(dst, dst_size, &used, &body, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_system: {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      has_system = 1;
      if (render_system_section(dst, dst_size, &used, &body, cutoff, charge,
                                &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_cpmd: {
      struct CPMDCpmdSection body;
      read_CPMDCpmdSection(&body, sec.cpmd);
      has_cpmd = 1;
      if (render_cpmd_section(dst, dst_size, &used, &body, task, file_path,
                              &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_dft: {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      has_dft = 1;
      if (render_dft_section(dst, dst_size, &used, &body, functional, mult,
                             &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_atoms: {
      struct CPMDAtomsSection body;
      read_CPMDAtomsSection(&body, sec.atoms);
      has_atoms = 1;
      if (render_atoms_section(dst, dst_size, &used, &body, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_set: {
      break;
    }
    case CPMDInputSection_raw:
      if (append_capn_text(dst, dst_size, &used, sec.raw) != 0)
        return -1;
      if (append_text(dst, dst_size, &used, "\n\n") != 0)
        return -1;
      break;
    case CPMDInputSection_atom:
    case CPMDInputSection_basis:
    case CPMDInputSection_clas:
    case CPMDInputSection_eam:
    case CPMDInputSection_exte:
    case CPMDInputSection_hardness:
    case CPMDInputSection_info:
    case CPMDInputSection_linres:
    case CPMDInputSection_molstates:
    case CPMDInputSection_mts:
    case CPMDInputSection_nlcc:
    case CPMDInputSection_path:
    case CPMDInputSection_pimd:
    case CPMDInputSection_potential:
    case CPMDInputSection_prop:
    case CPMDInputSection_ptddft:
    case CPMDInputSection_resp:
    case CPMDInputSection_tddft:
    case CPMDInputSection_vdw:
    case CPMDInputSection_vectors:
    case CPMDInputSection_wavefunction:
      if (render_typed_directive_section(dst, dst_size, &used, &sec, &sets) != 0)
        return -1;
      break;
    default:
      return -1;
    }
  }

  if (!has_cpmd) {
    struct CPMDCpmdSection def;
    memset(&def, 0, sizeof(def));
    def.optimizeWavefunction = 1;
    def.convergenceOrbitals = 1.0e-6;
    if (render_cpmd_section(dst, dst_size, &used, &def, task, file_path,
                            &sets) != 0)
      return -1;
  }
  if (!has_system) {
    struct CPMDSystemSection def;
    memset(&def, 0, sizeof(def));
    def.symmetry = 0;
    def.angstrom = 1;
    def.cutOffRy = cutoff;
    def.charge = charge;
    def.multiplicity = mult;
    if (render_system_section(dst, dst_size, &used, &def, cutoff, charge,
                              &sets) != 0)
      return -1;
  }
  if (!has_dft) {
    struct CPMDDftSection def;
    memset(&def, 0, sizeof(def));
    def.lsd = mult > 1;
    if (render_dft_section(dst, dst_size, &used, &def, functional, mult,
                           &sets) != 0)
      return -1;
  }
  if (!has_atoms) {
    struct CPMDAtomsSection def;
    memset(&def, 0, sizeof(def));
    if (render_atoms_section(dst, dst_size, &used, &def, &sets) != 0)
      return -1;
  }
  if (render_remaining_set_sections(dst, dst_size, &used, &sets) != 0)
    return -1;

  return 0;
}

static int element_z(const char *sym, size_t len) {
  if (!sym || len == 0)
    return -1;
  static const char *const symbols[] = {
      "H",  "HE", "LI", "BE", "B",  "C",  "N",  "O",  "F",  "NE",
      "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR", "K",  "CA",
      "SC", "TI", "V",  "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
      "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y",  "ZR",
      "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
      "SB", "TE", "I",  "XE", "CS", "BA", "LA", "CE", "PR", "ND",
      "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
      "LU", "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG",
      "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
      "PA", "U",  "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
      "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS",
      "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG",
  };
  char normalized[3] = {'\0', '\0', '\0'};
  size_t used = 0;
  for (size_t i = 0; i < len && used < 2; ++i) {
    unsigned char c = (unsigned char)sym[i];
    if (c == '\0' || isspace(c))
      break;
    if (!isalpha(c))
      return -1;
    normalized[used++] = (char)ascii_upper(c);
  }
  if (used == 0)
    return -1;
  for (size_t i = 0; i < sizeof(symbols) / sizeof(symbols[0]); ++i) {
    if (strcmp(normalized, symbols[i]) == 0)
      return (int)i + 1;
  }
  return -1;
}

static int render_atoms_with_geometry(char *dst, size_t dst_size, size_t *used,
                                      struct CPMDAtomsSection *atoms,
                                      int n_atoms, const double *pos,
                                      const int *z,
                                      struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&ATOMS\n") != 0)
    return -1;
  int npsp = struct_list_len(&atoms->pseudopotentials.p);
  if (npsp < 0)
    return -1;
  if (npsp == 0) {
    /* Default BLYP PPs for H/O when Cap'n Proto omits atoms list. */
    static const struct {
      int z;
      const char *path;
      const char *lmax;
    } defs[] = {{8, "O_MT_BLYP.psp", "P"}, {1, "H_CVB_BLYP.psp", "S"}};
    int covered = 0;
    for (size_t d = 0; d < sizeof(defs) / sizeof(defs[0]); ++d) {
      int count = 0;
      for (int j = 0; j < n_atoms; ++j)
        if (z[j] == defs[d].z)
          ++count;
      if (count == 0)
        continue;
      covered += count;
      if (append_fmt(dst, dst_size, used, "*%s\n LMAX=%s\n   %d\n", defs[d].path,
                     defs[d].lmax, count) != 0)
        return -1;
      for (int j = 0; j < n_atoms; ++j) {
        if (z[j] != defs[d].z)
          continue;
        if (append_fmt(dst, dst_size, used, "  %.6f  %.6f  %.6f\n", pos[3 * j],
                       pos[3 * j + 1], pos[3 * j + 2]) != 0)
          return -1;
      }
    }
    if (covered != n_atoms)
      return -1;
    if (append_directives(dst, dst_size, used, atoms->directives) != 0)
      return -1;
    if (append_set_directives_for_section(dst, dst_size, used, sets, "ATOMS") !=
        0)
      return -1;
    return append_text(dst, dst_size, used, "&END\n\n");
  }
  int covered = 0;
  for (int i = 0; i < npsp; ++i) {
    struct CPMDAtomsPseudopotential psp;
    get_CPMDAtomsPseudopotential(&psp, atoms->pseudopotentials, i);
    int zz = element_z(psp.element.str, (size_t)psp.element.len);
    if (zz < 0)
      return -1;
    int count = 0;
    for (int j = 0; j < n_atoms; ++j)
      if (z[j] == zz)
        ++count;
    if (count == 0)
      continue;
    covered += count;
    if (append_text(dst, dst_size, used, "*") != 0)
      return -1;
    if (psp.path.str && psp.path.len > 0) {
      if (append_capn_text(dst, dst_size, used, psp.path) != 0)
        return -1;
    } else {
      return -1;
    }
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    if (append_fmt(dst, dst_size, used, " LMAX=%s\n   %d\n",
                   lmax_letter(psp.lmax >= 0 ? psp.lmax : (zz == 8 ? 1 : 0)),
                   count) != 0)
      return -1;
    for (int j = 0; j < n_atoms; ++j) {
      if (z[j] != zz)
        continue;
      if (append_fmt(dst, dst_size, used, "  %.6f  %.6f  %.6f\n", pos[3 * j],
                     pos[3 * j + 1], pos[3 * j + 2]) != 0)
        return -1;
    }
  }
  if (covered != n_atoms)
    return -1;
  if (append_directives(dst, dst_size, used, atoms->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "ATOMS") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

int cpmdc_params_render_deck_with_geometry(
    CPMDParams_ptr params, int n_atoms, const double *positions_ang,
    const int *atomic_numbers, const double *cell_ang, int has_cell, char *dst,
    size_t dst_size) {
  if (params.p.type == CAPN_NULL || !dst || dst_size == 0 || n_atoms <= 0 ||
      !positions_ang || !atomic_numbers)
    return -1;
  dst[0] = '\0';
  size_t used = 0;
  struct CPMDParams view;
  read_CPMDParams(&view, params);

  if (view.title.str && view.title.len > 0) {
    if (append_text(dst, dst_size, &used, "! ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, &used, view.title) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  int nblocks = pointer_list_len(&view.inputBlocks);
  if (nblocks < 0)
    return -1;
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(view.inputBlocks, i, empty_text);
    if (append_capn_text(dst, dst_size, &used, block) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  const char *task = cpmdc_params_text_or(view.task, "gradient");
  const char *functional = cpmdc_params_text_or(view.functional, "BLYP");
  capn_text file_path = selected_file_path(&view);
  double cutoff = view.cutOffRy > 0.0 ? view.cutOffRy : 70.0;
  int charge = view.charge;
  int mult = view.multiplicity > 0 ? view.multiplicity : 1;

  int has_system = 0, has_cpmd = 0, has_dft = 0, has_atoms = 0;
  int nsec = struct_list_len(&view.inputSections.p);
  if (nsec < 0)
    return -1;
  struct RenderSetList sets;
  if (collect_set_directives(view.inputSections, &sets) != 0)
    return -1;

  struct CPMDAtomsSection atoms_body;
  memset(&atoms_body, 0, sizeof(atoms_body));

  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, view.inputSections, i);
    switch (sec.which) {
    case CPMDInputSection_generic: {
      struct CPMDGenericSection body;
      read_CPMDGenericSection(&body, sec.generic);
      if (render_generic_section(dst, dst_size, &used, &body, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_system: {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      has_system = 1;
      if (render_system_section_with_cell(
              dst, dst_size, &used, &body, cutoff, charge,
              has_cell ? cell_ang : NULL, has_cell ? 9 : 0, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_cpmd: {
      struct CPMDCpmdSection body;
      read_CPMDCpmdSection(&body, sec.cpmd);
      has_cpmd = 1;
      if (render_cpmd_section(dst, dst_size, &used, &body, task, file_path,
                              &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_dft: {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      has_dft = 1;
      if (render_dft_section(dst, dst_size, &used, &body, functional, mult,
                             &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_atoms: {
      read_CPMDAtomsSection(&atoms_body, sec.atoms);
      has_atoms = 1;
      if (render_atoms_with_geometry(dst, dst_size, &used, &atoms_body, n_atoms,
                                     positions_ang, atomic_numbers, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_set: {
      break;
    }
    case CPMDInputSection_raw:
      if (append_capn_text(dst, dst_size, &used, sec.raw) != 0)
        return -1;
      if (append_text(dst, dst_size, &used, "\n\n") != 0)
        return -1;
      break;
    case CPMDInputSection_atom:
    case CPMDInputSection_basis:
    case CPMDInputSection_clas:
    case CPMDInputSection_eam:
    case CPMDInputSection_exte:
    case CPMDInputSection_hardness:
    case CPMDInputSection_info:
    case CPMDInputSection_linres:
    case CPMDInputSection_molstates:
    case CPMDInputSection_mts:
    case CPMDInputSection_nlcc:
    case CPMDInputSection_path:
    case CPMDInputSection_pimd:
    case CPMDInputSection_potential:
    case CPMDInputSection_prop:
    case CPMDInputSection_ptddft:
    case CPMDInputSection_resp:
    case CPMDInputSection_tddft:
    case CPMDInputSection_vdw:
    case CPMDInputSection_vectors:
    case CPMDInputSection_wavefunction:
      if (render_typed_directive_section(dst, dst_size, &used, &sec, &sets) != 0)
        return -1;
      break;
    default:
      return -1;
    }
  }

  if (!has_cpmd) {
    struct CPMDCpmdSection def;
    memset(&def, 0, sizeof(def));
    def.optimizeWavefunction = 1;
    def.convergenceOrbitals = 1.0e-5;
    if (render_cpmd_section(dst, dst_size, &used, &def, task, file_path,
                            &sets) != 0)
      return -1;
  }
  if (!has_system) {
    struct CPMDSystemSection def;
    memset(&def, 0, sizeof(def));
    def.symmetry = 0;
    def.angstrom = 1;
    def.cutOffRy = cutoff;
    def.charge = charge;
    def.multiplicity = mult;
    if (render_system_section_with_cell(dst, dst_size, &used, &def, cutoff,
                                        charge, has_cell ? cell_ang : NULL,
                                        has_cell ? 9 : 0, &sets) != 0)
      return -1;
  }
  if (!has_dft) {
    struct CPMDDftSection def;
    memset(&def, 0, sizeof(def));
    def.lsd = mult > 1;
    if (render_dft_section(dst, dst_size, &used, &def, functional, mult,
                           &sets) != 0)
      return -1;
  }
  if (!has_atoms) {
    if (render_atoms_with_geometry(dst, dst_size, &used, &atoms_body, n_atoms,
                                   positions_ang, atomic_numbers, &sets) != 0)
      return -1;
  }
  if (render_remaining_set_sections(dst, dst_size, &used, &sets) != 0)
    return -1;
  (void)has_cell;
  return 0;
}

static int force_input_length_factor(capn_text unit, double *factor) {
  if (!factor)
    return -1;
  if (!unit.str || unit.len <= 0 || text_eq_ci(unit, "angstrom") ||
      text_eq_ci(unit, "angstroms") || text_eq_ci(unit, "a") ||
      text_eq_ci(unit, "aa")) {
    *factor = 1.0;
    return 0;
  }
  if (text_eq_ci(unit, "bohr") || text_eq_ci(unit, "au") ||
      text_eq_ci(unit, "a.u.") || text_eq_ci(unit, "atomic")) {
    *factor = CPMDC_BOHR_TO_ANGSTROM;
    return 0;
  }
  if (text_eq_ci(unit, "nm") || text_eq_ci(unit, "nanometer")) {
    *factor = 10.0;
    return 0;
  }
  return -1;
}

static int force_input_energy_factor(capn_text unit, double *factor) {
  if (!factor)
    return -1;
  if (!unit.str || unit.len <= 0 || text_eq_ci(unit, "hartree") ||
      text_eq_ci(unit, "ha") || text_eq_ci(unit, "au") ||
      text_eq_ci(unit, "a.u.")) {
    *factor = 1.0;
    return 0;
  }
  if (text_eq_ci(unit, "ev") || text_eq_ci(unit, "electronvolt")) {
    *factor = 27.211386245988;
    return 0;
  }
  if (text_eq_ci(unit, "ry") || text_eq_ci(unit, "rydberg")) {
    *factor = 2.0;
    return 0;
  }
  if (text_eq_ci(unit, "kj/mol") || text_eq_ci(unit, "kjmol")) {
    *factor = 2625.499639479;
    return 0;
  }
  if (text_eq_ci(unit, "kcal/mol") || text_eq_ci(unit, "kcalmol")) {
    *factor = 627.5094740631;
    return 0;
  }
  return -1;
}

int cpmdc_force_input_root(const void *force_input_capnp,
                           size_t force_input_capnp_size_bytes,
                           struct capn *arena, ForceInput_ptr *force_input) {
  if (!force_input_capnp || force_input_capnp_size_bytes == 0 || !arena ||
      !force_input)
    return -1;
  memset(arena, 0, sizeof(*arena));
  memset(force_input, 0, sizeof(*force_input));
  if (capn_init_mem(arena, (const uint8_t *)force_input_capnp,
                    force_input_capnp_size_bytes, 0) != 0)
    return -1;
  force_input->p = capn_getp(capn_root(arena), 0, 1);
  if (force_input->p.type != CAPN_STRUCT) {
    cpmdc_params_release(arena);
    memset(force_input, 0, sizeof(*force_input));
    return -1;
  }
  return 0;
}

int cpmdc_force_input_atom_count(ForceInput_ptr force_input, size_t *n_atoms,
                                 int *has_cell) {
  if (force_input.p.type == CAPN_NULL || !n_atoms || !has_cell)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  int n_pos = list64_len(&view.pos);
  int n_z = list32_len(&view.atmnrs);
  int n_box = list64_len(&view.box);
  if (n_pos < 0 || n_z <= 0 || n_box < 0)
    return -1;
  if (n_pos != n_z * 3)
    return -1;
  if (n_box != 0 && n_box != 9)
    return -1;
  *n_atoms = (size_t)n_z;
  *has_cell = n_box == 9 ? 1 : 0;
  return 0;
}

int cpmdc_force_input_copy_geometry(ForceInput_ptr force_input,
                                    double *positions_ang, int *atomic_numbers,
                                    size_t atom_capacity, double *cell_ang,
                                    int *has_cell) {
  if (force_input.p.type == CAPN_NULL || !positions_ang || !atomic_numbers ||
      !cell_ang || !has_cell)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  size_t n_atoms = 0;
  int local_has_cell = 0;
  if (cpmdc_force_input_atom_count(force_input, &n_atoms, &local_has_cell) != 0)
    return -1;
  if (atom_capacity < n_atoms)
    return -1;
  double length_factor = 1.0;
  if (force_input_length_factor(view.lengthUnit, &length_factor) != 0)
    return -1;
  capn_resolve(&view.pos.p);
  capn_resolve(&view.atmnrs.p);
  capn_resolve(&view.box.p);
  for (size_t i = 0; i < n_atoms; ++i)
    atomic_numbers[i] = (int)(int32_t)capn_get32(view.atmnrs, (int)i);
  for (size_t i = 0; i < n_atoms * 3u; ++i)
    positions_ang[i] =
        capn_to_f64(capn_get64(view.pos, (int)i)) * length_factor;
  for (size_t i = 0; i < 9; ++i)
    cell_ang[i] = local_has_cell
                      ? capn_to_f64(capn_get64(view.box, (int)i)) * length_factor
                      : 0.0;
  *has_cell = local_has_cell;
  return 0;
}

int cpmdc_force_input_result_factors(ForceInput_ptr force_input,
                                     double *energy_factor,
                                     double *force_factor) {
  if (force_input.p.type == CAPN_NULL || !energy_factor || !force_factor)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  double length_factor = 1.0;
  double energy = 1.0;
  if (force_input_length_factor(view.lengthUnit, &length_factor) != 0 ||
      force_input_energy_factor(view.energyUnit, &energy) != 0)
    return -1;
  *energy_factor = energy;
  *force_factor = energy * length_factor / CPMDC_BOHR_TO_ANGSTROM;
  return 0;
}

size_t cpmdc_potential_result_flat_size(size_t force_count) {
  if (force_count > (SIZE_MAX - 32u) / 8u)
    return 0;
  return 32u + force_count * 8u;
}

int cpmdc_potential_result_write(double energy, const double *forces,
                                 size_t force_count,
                                 void *potential_result_capnp,
                                 size_t potential_result_capacity_bytes,
                                 size_t *potential_result_size_bytes) {
  if (!forces || !potential_result_capnp || !potential_result_size_bytes ||
      force_count > (size_t)INT_MAX)
    return -1;
  size_t required = cpmdc_potential_result_flat_size(force_count);
  *potential_result_size_bytes = required;
  if (required == 0 || potential_result_capacity_bytes < required)
    return -1;

  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  if (root.type == CAPN_NULL) {
    capn_free(&arena);
    return -1;
  }
  PotentialResult_ptr result = new_PotentialResult(root.seg);
  capn_list64 force_list = capn_new_list64(root.seg, (int)force_count);
  if (result.p.type == CAPN_NULL ||
      (force_count > 0 && force_list.p.type == CAPN_NULL)) {
    capn_free(&arena);
    return -1;
  }
  for (size_t i = 0; i < force_count; ++i)
    capn_set64(force_list, (int)i, capn_from_f64(forces[i]));

  struct PotentialResult view;
  view.energy = energy;
  view.forces = force_list;
  write_PotentialResult(&view, result);
  if (capn_setp(capn_root(&arena), 0, result.p) != 0) {
    capn_free(&arena);
    return -1;
  }
  int wrote = capn_write_mem(&arena, (uint8_t *)potential_result_capnp,
                             (int)potential_result_capacity_bytes, 0);
  capn_free(&arena);
  if (wrote < 0 || (size_t)wrote > potential_result_capacity_bytes)
    return -1;
  *potential_result_size_bytes = (size_t)wrote;
  return 0;
}
