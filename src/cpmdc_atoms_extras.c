#include "cpmdc_atoms_extras.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *const cpmdc_atoms_dropped_directives[] = {
    "VELOCITIES",
    "METADYNAMICS",
    NULL,
};

enum line_kind {
  LINE_NONE = 0,
  LINE_SECTION_END,
  LINE_END_BLOCK,
  LINE_CONSTRAINTS,
  LINE_ISOTOPE,
  LINE_DUMMY,
  LINE_CHANGE,
  LINE_CONFIN,
  LINE_CHARGES,
  LINE_MOVIE,
  LINE_DROPPED,
  LINE_GENERATE
};

enum block_kind {
  BLOCK_NONE = 0,
  BLOCK_KEEP,
  BLOCK_KEEP_END,
  BLOCK_DROP,
  BLOCK_SKIP
};

struct grow {
  char *data;
  int len;
  int cap;
};

static int set_err(char *err, int err_cap, const char *msg) {
  if (err && err_cap > 0)
    snprintf(err, (size_t)err_cap, "%s", msg ? msg : "");
  return -1;
}

static char fold(char c) {
  if (c >= 'a' && c <= 'z')
    return (char)(c - 32);
  return c;
}

static int ci_contains(const char *s, int n, const char *lit) {
  int m = (int)strlen(lit);
  int i, k;
  if (m <= 0 || m > n)
    return 0;
  for (i = 0; i + m <= n; ++i) {
    for (k = 0; k < m; ++k) {
      if (fold(s[i + k]) != fold(lit[k]))
        break;
    }
    if (k == m)
      return 1;
  }
  return 0;
}

static int starts_ci(const char *s, int n, const char *lit) {
  int m = (int)strlen(lit);
  int i = 0;
  int k;
  while (i < n && (s[i] == ' ' || s[i] == '\t'))
    i++;
  if (n - i < m)
    return 0;
  for (k = 0; k < m; ++k) {
    if (fold(s[i + k]) != fold(lit[k]))
      return 0;
  }
  if (i + m < n && s[i + m] != ' ' && s[i + m] != '\t')
    return 0;
  return 1;
}

static int is_dropped_name(const char *s, int n) {
  int i;
  for (i = 0; cpmdc_atoms_dropped_directives[i] != NULL; ++i) {
    if (ci_contains(s, n, cpmdc_atoms_dropped_directives[i]))
      return 1;
  }
  /* OpenCPMD also accepts the two-word spelling in ratom. */
  if (ci_contains(s, n, "META") && ci_contains(s, n, "DYNAM"))
    return 1;
  return 0;
}

static int is_comment_line(const char *s, int n) {
  int i = 0;
  while (i < n && (s[i] == ' ' || s[i] == '\t'))
    i++;
  return i < n && s[i] == '!';
}

static int classify_line(const char *s, int n) {
  if (is_comment_line(s, n))
    return LINE_NONE;
  if (starts_ci(s, n, "&END"))
    return LINE_SECTION_END;
  if (starts_ci(s, n, "END") &&
      (is_dropped_name(s, n) || ci_contains(s, n, "CONSTRAINT")))
    return LINE_END_BLOCK;
  if (is_dropped_name(s, n))
    return LINE_DROPPED;
  if (ci_contains(s, n, "CONSTRAINT"))
    return LINE_CONSTRAINTS;
  if (ci_contains(s, n, "ISOTOPE"))
    return LINE_ISOTOPE;
  if (ci_contains(s, n, "DUMMY") && ci_contains(s, n, "ATOM"))
    return LINE_DUMMY;
  if (ci_contains(s, n, "CHANGE") && ci_contains(s, n, "BOND"))
    return LINE_CHANGE;
  if (ci_contains(s, n, "CONFIN") && ci_contains(s, n, "POTENT"))
    return LINE_CONFIN;
  if (ci_contains(s, n, "ATOMIC") && ci_contains(s, n, "CHARGE"))
    return LINE_CHARGES;
  if (ci_contains(s, n, "MOVIE") && ci_contains(s, n, "TYPE"))
    return LINE_MOVIE;
  if (ci_contains(s, n, "GENER") && ci_contains(s, n, "COORD"))
    return LINE_GENERATE;
  return LINE_NONE;
}

static int find_ci(const char *s, int n, int from, const char *lit) {
  int m = (int)strlen(lit);
  int i, k;
  if (from < 0)
    from = 0;
  if (m <= 0 || from > n)
    return -1;
  for (i = from; i + m <= n; ++i) {
    for (k = 0; k < m; ++k) {
      if (fold(s[i + k]) != fold(lit[k]))
        break;
    }
    if (k == m)
      return i;
  }
  return -1;
}

static int take_line(const char *deck, int end, int pos, int *lb, int *le) {
  int i = pos;
  *lb = pos;
  while (i < end && deck[i] != '\n' && deck[i] != '\r')
    i++;
  *le = i;
  if (i < end && deck[i] == '\r')
    i++;
  if (i < end && deck[i] == '\n')
    i++;
  return i;
}

static int grow_add(struct grow *g, const char *s, int n) {
  int need;
  int cap;
  char *p;
  if (n < 0)
    return -1;
  if (n == 0)
    return 0;
  if (g->len > INT_MAX - n - 1)
    return -1;
  need = g->len + n + 1;
  if (need > g->cap) {
    cap = g->cap > 0 ? g->cap : 64;
    while (cap < need) {
      if (cap > INT_MAX / 2)
        return -1;
      cap *= 2;
    }
    p = (char *)realloc(g->data, (size_t)cap);
    if (!p)
      return -1;
    g->data = p;
    g->cap = cap;
  }
  memcpy(g->data + g->len, s, (size_t)n);
  g->len += n;
  g->data[g->len] = '\0';
  return 0;
}

static int all_digits(const char *s, int n) {
  int i;
  if (n < 1)
    return 0;
  for (i = 0; i < n; ++i) {
    if (s[i] < '0' || s[i] > '9')
      return 0;
  }
  return 1;
}

static int parse_u32(const char *s, int n, int *value) {
  int i;
  int v = 0;
  if (!all_digits(s, n) || n > 9)
    return 0;
  for (i = 0; i < n; ++i)
    v = v * 10 + (s[i] - '0');
  *value = v;
  return 1;
}

static int lone_integer_line(const char *s, int n) {
  int i = 0;
  int start;
  int tokens = 0;
  while (i < n && (s[i] == ' ' || s[i] == '\t'))
    i++;
  if (i >= n)
    return 0;
  start = i;
  while (i < n && s[i] != ' ' && s[i] != '\t')
    i++;
  if (!all_digits(s + start, i - start))
    return 0;
  tokens = 1;
  while (i < n && (s[i] == ' ' || s[i] == '\t'))
    i++;
  if (i < n)
    tokens = 2;
  return tokens == 1;
}

/* Same numeric-token rule as deck_has_real_atoms: three or more numeric
 * tokens, and a '*' or '=' line is never a coordinate line. */
static int is_coord_line(const char *s, int n) {
  int i;
  int ntok = 0;
  if (ci_contains(s, n, "*") || ci_contains(s, n, "="))
    return 0;
  i = 0;
  while (i < n) {
    int digit = 0;
    int ok = 1;
    int j;
    if (s[i] == ' ' || s[i] == '\t' || s[i] == ',') {
      i++;
      continue;
    }
    j = i;
    while (j < n && s[j] != ' ' && s[j] != '\t' && s[j] != ',') {
      if (s[j] >= '0' && s[j] <= '9')
        digit = 1;
      else if (!(s[j] == '.' || s[j] == '+' || s[j] == '-' || s[j] == 'e' ||
                 s[j] == 'E' || s[j] == 'd' || s[j] == 'D'))
        ok = 0;
      j++;
    }
    if (ok && digit)
      ntok++;
    i = j;
  }
  return ntok >= 3;
}

static int line_empty(const char *s, int n) {
  int i;
  for (i = 0; i < n; ++i) {
    if (s[i] != ' ' && s[i] != '\t' && s[i] != '\r')
      return 0;
  }
  return 1;
}

static int starts_star(const char *s, int n) {
  int i = 0;
  while (i < n && (s[i] == ' ' || s[i] == '\t'))
    i++;
  return i < n && s[i] == '*';
}

static int emit_raw(struct grow *g, const char *s, int n) {
  if (grow_add(g, s, n) != 0)
    return -1;
  return grow_add(g, "\n", 1);
}

static int emit_remapped(struct grow *g, const char *s, int n, const int *map,
                         int n_atoms, const char *directive, char *err,
                         int err_cap) {
  int i = 0;
  if (lone_integer_line(s, n))
    return emit_raw(g, s, n);
  while (i < n) {
    int start;
    int value;
    int neu;
    char num[16];
    int m;
    if (s[i] == ' ' || s[i] == '\t') {
      if (grow_add(g, s + i, 1) != 0)
        return -1;
      i++;
      continue;
    }
    start = i;
    while (i < n && s[i] != ' ' && s[i] != '\t')
      i++;
    if (all_digits(s + start, i - start) &&
        parse_u32(s + start, i - start, &value)) {
      if (value >= 1 && value <= n_atoms) {
        neu = map[value - 1];
        if (neu < 1) {
          char msg[128];
          snprintf(msg, sizeof(msg),
                   "%s atom index is outside the rebuilt order", directive);
          return set_err(err, err_cap, msg);
        }
        m = snprintf(num, sizeof(num), "%d", neu);
        if (m < 0 || grow_add(g, num, m) != 0)
          return set_err(err, err_cap, "out of memory");
        continue;
      }
        if (value == 0) {
        if (grow_add(g, s + start, i - start) != 0)
          return set_err(err, err_cap, "out of memory");
        continue;
      }
      {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "%s atom index is outside the rebuilt order", directive);
        return set_err(err, err_cap, msg);
      }
    }
    if (grow_add(g, s + start, i - start) != 0)
      return set_err(err, err_cap, "out of memory");
  }
  if (grow_add(g, "\n", 1) != 0)
    return set_err(err, err_cap, "out of memory");
  return 0;
}

static int build_old_to_new(const int *z, int n_atoms, int *map, char *err,
                            int err_cap) {
  unsigned char seen[121];
  int next = 1;
  int i, j, zz;
  memset(seen, 0, sizeof(seen));
  for (i = 0; i < n_atoms; ++i)
    map[i] = 0;
  for (i = 0; i < n_atoms; ++i) {
    zz = z[i];
    if (zz < 1 || zz > 120)
      return set_err(err, err_cap, "atomic number is outside 1..120");
    if (seen[zz])
      continue;
    seen[zz] = 1;
    for (j = 0; j < n_atoms; ++j) {
      if (z[j] != zz)
        continue;
      map[j] = next++;
    }
  }
  if (next != n_atoms + 1)
    return set_err(err, err_cap, "atomic numbers do not cover the step");
  return 0;
}

static int keep_block(struct grow *g, const char *deck, int from, int end,
                      const int *map, int n_atoms, char *err, int err_cap) {
  int pos = from;
  int block = BLOCK_NONE;
  int remap = 0;
  int skip_left = 0;
  int expect_count = 0;
  const char *name = "";
  while (pos < end) {
    int lb = 0;
    int le = 0;
    int next;
    int kind;
    const char *line;
    int n;
    int consumed;
    next = take_line(deck, end, pos, &lb, &le);
    if (next <= pos)
      break;
    line = deck + lb;
    n = le - lb;
    for (;;) {
      kind = classify_line(line, n);
      consumed = 1;
      if (kind == LINE_SECTION_END) {
        pos = end;
        break;
      }
      if (skip_left > 0 && block == BLOCK_NONE) {
        skip_left--;
        break;
      }
      if (block == BLOCK_DROP || block == BLOCK_SKIP) {
        if (block == BLOCK_DROP && kind == LINE_END_BLOCK &&
            is_dropped_name(line, n)) {
          block = BLOCK_NONE;
          break;
        }
        if (kind != LINE_NONE) {
          block = BLOCK_NONE;
          consumed = 0;
        }
        break;
      }
      if (block == BLOCK_KEEP_END) {
        if (kind == LINE_END_BLOCK && ci_contains(line, n, "CONSTRAINT")) {
          if (emit_raw(g, line, n) != 0)
            return set_err(err, err_cap, "out of memory");
          block = BLOCK_NONE;
          break;
        }
        if (remap) {
          if (emit_remapped(g, line, n, map, n_atoms, name, err, err_cap) != 0)
            return -1;
        } else if (emit_raw(g, line, n) != 0) {
          return set_err(err, err_cap, "out of memory");
        }
        break;
      }
      if (block == BLOCK_KEEP) {
        if (kind != LINE_NONE) {
          block = BLOCK_NONE;
          consumed = 0;
          break;
        }
        if (remap) {
          if (emit_remapped(g, line, n, map, n_atoms, name, err, err_cap) != 0)
            return -1;
        } else if (emit_raw(g, line, n) != 0) {
          return set_err(err, err_cap, "out of memory");
        }
        break;
      }
      if (kind == LINE_CONSTRAINTS) {
        if (emit_raw(g, line, n) != 0)
          return set_err(err, err_cap, "out of memory");
        block = BLOCK_KEEP_END;
        remap = 1;
        name = "CONSTRAINTS";
        expect_count = 0;
        break;
      }
      if (kind == LINE_ISOTOPE || kind == LINE_CONFIN || kind == LINE_CHARGES ||
          kind == LINE_MOVIE) {
        if (emit_raw(g, line, n) != 0)
          return set_err(err, err_cap, "out of memory");
        block = BLOCK_KEEP;
        remap = 0;
        name = kind == LINE_ISOTOPE ? "ISOTOPE"
               : kind == LINE_CONFIN ? "CONFINEMENT POTENTIAL"
               : kind == LINE_CHARGES ? "ATOMIC CHARGES"
                                      : "MOVIE TYPE";
        expect_count = 0;
        break;
      }
      if (kind == LINE_DUMMY || kind == LINE_CHANGE) {
        if (emit_raw(g, line, n) != 0)
          return set_err(err, err_cap, "out of memory");
        block = BLOCK_KEEP;
        remap = 1;
        name = kind == LINE_DUMMY ? "DUMMY ATOMS" : "CHANGE BONDS";
        expect_count = 0;
        break;
      }
      if (kind == LINE_DROPPED) {
        block = BLOCK_DROP;
        expect_count = 0;
        break;
      }
      if (kind == LINE_GENERATE) {
        block = BLOCK_SKIP;
        expect_count = 0;
        break;
      }
      if (kind == LINE_END_BLOCK || is_comment_line(line, n) ||
          starts_star(line, n) || line_empty(line, n)) {
        expect_count = 0;
        break;
      }
      if (ci_contains(line, n, "=")) {
        expect_count = 1;
        break;
      }
      if (expect_count && lone_integer_line(line, n)) {
        int count = 0;
        int i = 0;
        int start;
        while (i < n && (line[i] == ' ' || line[i] == '\t'))
          i++;
        start = i;
        while (i < n && line[i] != ' ' && line[i] != '\t')
          i++;
        if (!parse_u32(line + start, i - start, &count) || count < 0)
          count = 0;
        skip_left = count;
        expect_count = 0;
        break;
      }
      if (is_coord_line(line, n) || lone_integer_line(line, n)) {
        expect_count = 0;
        break;
      }
      expect_count = 0;
      break;
    }
    if (consumed)
      pos = next;
  }
  return 0;
}

int cpmdc_atoms_extras_keep(const char *deck, int deck_len, const int *z,
                            int n_atoms, char *out, int out_cap, int *out_len,
                            char *err, int err_cap) {
  struct grow g;
  int *map = NULL;
  int from;
  int rc;
  memset(&g, 0, sizeof(g));
  if (out_len)
    *out_len = 0;
  if (!deck || deck_len < 0 || n_atoms < 0 || !out_len)
    return set_err(err, err_cap, "invalid atoms extras arguments");
  if (n_atoms > 0 && !z)
    return set_err(err, err_cap, "atomic numbers are missing");
  if (n_atoms > 0) {
    map = (int *)calloc((size_t)n_atoms, sizeof(int));
    if (!map)
      return set_err(err, err_cap, "out of memory");
    if (build_old_to_new(z, n_atoms, map, err, err_cap) != 0) {
      free(map);
      return -1;
    }
  }
  from = 0;
  while (from < deck_len) {
    int at = find_ci(deck, deck_len, from, "&ATOMS");
    int end;
    if (at < 0)
      break;
    end = find_ci(deck, deck_len, at + 6, "&END");
    if (end < 0)
      break;
    rc = keep_block(&g, deck, at + 6, end, map, n_atoms, err, err_cap);
    if (rc != 0) {
      free(map);
      free(g.data);
      *out_len = 0;
      return -1;
    }
    from = end + 4;
  }
  free(map);
  if (g.len == 0) {
    if (out && out_cap > 0)
      out[0] = '\0';
    *out_len = 0;
    free(g.data);
    if (err && err_cap > 0)
      err[0] = '\0';
    return 0;
  }
  if (!out || out_cap < g.len + 1) {
    *out_len = g.len;
    free(g.data);
    return set_err(err, err_cap, "atoms extras do not fit");
  }
  memcpy(out, g.data, (size_t)g.len);
  out[g.len] = '\0';
  *out_len = g.len;
  free(g.data);
  if (err && err_cap > 0)
    err[0] = '\0';
  return 0;
}
