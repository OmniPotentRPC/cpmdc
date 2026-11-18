#define _GNU_SOURCE
#include "cpmdc_restart.h"

#include <dirent.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef CPMDC_RESTART_TEST_HOOK
/* Defined in src/cpmdc_restart.c only when this test is built. */
extern void (*cpmdc_restart_test_before_follow)(void);
static char hook_swap_path[PATH_MAX];
static int hook_swap_armed;
static int hook_swap_ok;

/* Replace one intermediate symlink with a regular file in the window
 * after realpath has reported the chain missing. */
static void swap_intermediate_for_file(void) {
  int fd;
  ssize_t n;
  if (!hook_swap_armed)
    return;
  hook_swap_armed = 0;
  if (unlink(hook_swap_path) != 0)
    return;
  fd = open(hook_swap_path, O_WRONLY | O_CREAT | O_EXCL, 0644);
  if (fd < 0)
    return;
  n = write(fd, "KEEP", 4);
  if (n == 4 && close(fd) == 0)
    hook_swap_ok = 1;
  else
    close(fd);
}
#endif

static int fail(const char *msg) {
  fprintf(stderr, "test_cpmd_restart: %s\n", msg);
  return 1;
}

static int add_bytes(unsigned char **buf, size_t *len, size_t *cap, const void *src,
                     size_t n) {
  if (*len + n > *cap) {
    size_t ncap = *cap ? *cap * 2 : 256;
    while (ncap < *len + n)
      ncap *= 2;
    unsigned char *grown = realloc(*buf, ncap);
    if (!grown)
      return -1;
    *buf = grown;
    *cap = ncap;
  }
  if (n)
    memcpy(*buf + *len, src, n);
  *len += n;
  return 0;
}

static int add_rec(unsigned char **buf, size_t *len, size_t *cap, const void *src,
                   size_t n) {
  int32_t marker = (int32_t)n;
  if (add_bytes(buf, len, cap, &marker, 4) != 0)
    return -1;
  if (add_bytes(buf, len, cap, src, n) != 0)
    return -1;
  return add_bytes(buf, len, cap, &marker, 4);
}

static int add_i32(unsigned char **buf, size_t *len, size_t *cap, int32_t v) {
  return add_rec(buf, len, cap, &v, 4);
}

static int add_count(unsigned char **buf, size_t *len, size_t *cap, int32_t count,
                     int with_fpos) {
  unsigned char raw[12];
  int64_t fpos = 0;
  memcpy(raw, &count, 4);
  if (!with_fpos)
    return add_rec(buf, len, cap, raw, 4);
  memcpy(raw + 4, &fpos, 8);
  return add_rec(buf, len, cap, raw, 12);
}

static int build_sample(unsigned char **out, size_t *nout, int stream, int with_fpos) {
  unsigned char *buf = NULL;
  size_t len = 0, cap = 0;
  char header[80];
  int32_t sym[2];
  double cell[6] = {14.0, 1, 1, 0, 0, 0};
  int32_t nsp = 2;
  int32_t na[2] = {2, 1};
  double xyz[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
  double vel[9] = {0.1, 0, 0, 0, 0.2, 0, 0, 0, 0.3};
  unsigned char cut[36];
  int32_t states[6] = {8, 0, 100, 80, 20, 10};
  unsigned char wave_a[] = {1, 2, 3, 4, 5};
  unsigned char wave_b[] = {9, 9, 9};
  double ecut = 70.0, cdual = 4.0;
  int32_t dual_flag = 0, nel = 32, nr1 = 20, nr2 = 20, nr3 = 20;
  memset(header, ' ', sizeof header);
  memcpy(header, " STUTTGART VERSION 3.0", 22);
  if (stream)
    memcpy(header + 22, " STREAM", 7);
  sym[0] = 1;
  sym[1] = 0;
  memset(cut, 0, sizeof cut);
  memcpy(cut, &ecut, 8);
  memcpy(cut + 8, &cdual, 8);
  memcpy(cut + 16, &dual_flag, 4);
  memcpy(cut + 20, &nel, 4);
  memcpy(cut + 24, &nr1, 4);
  memcpy(cut + 28, &nr2, 4);
  memcpy(cut + 32, &nr3, 4);
  if (add_rec(&buf, &len, &cap, header, 80) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 2, with_fpos) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, sym, 8) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, cell, 48) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 2, with_fpos) != 0)
    return -1;
  if (add_i32(&buf, &len, &cap, nsp) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, na, 8) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 3, with_fpos) != 0)
    return -1;
  for (int a = 0; a < 3; a++) {
    if (add_rec(&buf, &len, &cap, xyz + 3 * a, 24) != 0)
      return -1;
  }
  if (add_count(&buf, &len, &cap, 3, with_fpos) != 0)
    return -1;
  for (int a = 0; a < 3; a++) {
    if (add_rec(&buf, &len, &cap, vel + 3 * a, 24) != 0)
      return -1;
  }
  if (add_count(&buf, &len, &cap, 0, with_fpos) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 1, with_fpos) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, cut, 36) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 1, with_fpos) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, states, 24) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, -2, with_fpos) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, wave_a, sizeof wave_a) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, wave_b, sizeof wave_b) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 0, with_fpos) != 0)
    return -1;
  *out = buf;
  *nout = len;
  return 0;
}

static int check_sample(const unsigned char *bytes, size_t n, int expect_stream) {
  char err[256];
  CPMDCRestart *file;
  double cell[6];
  double xyz[9];
  double vel[9];
  const int *na = NULL;
  int nsp = 0, ibrav = 0, indpg = 0;
  int nel = 0, n1 = 0, n2 = 0, n3 = 0, ns = 0, nk = 0, dual_flag = -1;
  double ecut = 0, cdual = 0;
  int sec9 = 0;
  void *rewritten = NULL;
  size_t rn = 0;
  double moved[9] = {2, 0, 0, 3, 0, 0, 4, 1, 0};
  const unsigned char *wave;
  err[0] = '\0';
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  if (!file)
    return fail(err[0] ? err : "read failed");
  if (cpmdc_restart_is_stream(file) != expect_stream) {
    cpmdc_restart_free(file);
    return fail("stream flag");
  }
  if (cpmdc_restart_ncoords(file) != 3) {
    cpmdc_restart_free(file);
    return fail("ncoords");
  }
  if (cpmdc_restart_cell(file, &ibrav, &indpg, cell) != 0 || ibrav != 1 || cell[0] != 14.0) {
    cpmdc_restart_free(file);
    return fail("cell");
  }
  if (cpmdc_restart_species(file, &nsp, &na) != 0 || nsp != 2 || na[0] != 2 || na[1] != 1) {
    cpmdc_restart_free(file);
    return fail("species");
  }
  {
    const char *species[] = {"H", "O"};
    const char *atoms[3] = {NULL, NULL, NULL};
    if (cpmdc_restart_expand_species_symbols(file, species, 2, atoms, 3) != 0 ||
        strcmp(atoms[0], "H") != 0 || strcmp(atoms[1], "H") != 0 ||
        strcmp(atoms[2], "O") != 0) {
      cpmdc_restart_free(file);
      return fail("species symbol expansion");
    }
  }
  if (cpmdc_restart_coordinates(file, xyz, 9) != 0 || xyz[3] != 1.0) {
    cpmdc_restart_free(file);
    return fail("coordinates");
  }
  if (cpmdc_restart_velocities(file, vel, 9) != 0 || vel[4] != 0.2) {
    cpmdc_restart_free(file);
    return fail("velocities");
  }
  if (cpmdc_restart_cutoff(file, &ecut, &cdual, &dual_flag, &nel, &n1, &n2, &n3) != 0 ||
      ecut != 70.0 || cdual != 4.0 || dual_flag != 0 || nel != 32 || n1 != 20) {
    cpmdc_restart_free(file);
    return fail("cutoff");
  }
  if (cpmdc_restart_states(file, &ns, &nk, NULL, NULL, NULL, NULL) != 0 || ns != 8 || nk != 0) {
    cpmdc_restart_free(file);
    return fail("states");
  }
  if (cpmdc_restart_section(file, 9, &sec9) != 0 || sec9 != -2) {
    cpmdc_restart_free(file);
    return fail("section 9 count");
  }
  if (cpmdc_restart_initial_coordinates(file, xyz, 3) == 0) {
    cpmdc_restart_free(file);
    return fail("empty section 6 accepted a coordinate triple");
  }
  if (cpmdc_restart_set_coordinates(file, moved, 6) == 0) {
    cpmdc_restart_free(file);
    return fail("short coordinate write should fail");
  }
  if (cpmdc_restart_set_coordinates(file, moved, 9) != 0) {
    cpmdc_restart_free(file);
    return fail("coordinate write");
  }
  if (cpmdc_restart_write_mem(file, &rewritten, &rn, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    return fail("rewrite");
  }
  cpmdc_restart_free(file);
  wave = memmem(rewritten, rn, "\x01\x02\x03\x04\x05", 5);
  if (!wave || !memmem(rewritten, rn, "\x09\x09\x09", 3)) {
    free(rewritten);
    return fail("wavefunction bytes were not preserved");
  }
  file = cpmdc_restart_read_mem(rewritten, rn, err, sizeof err);
  free(rewritten);
  if (!file)
    return fail(err[0] ? err : "reread failed");
  if (cpmdc_restart_coordinates(file, xyz, 9) != 0 || xyz[0] != 2.0 || xyz[6] != 4.0 ||
      xyz[8] != 0.0) {
    cpmdc_restart_free(file);
    return fail("patched coordinates did not round-trip");
  }
  if (cpmdc_restart_velocities(file, vel, 9) != 0 || vel[8] != 0.3) {
    cpmdc_restart_free(file);
    return fail("velocities changed during a coordinate patch");
  }
  cpmdc_restart_free(file);
  return 0;
}

static int check_split_header(void) {
  unsigned char *buf = NULL;
  size_t len = 0, cap = 0;
  char part1[10];
  char part2[70];
  int32_t m1 = -10, m2 = 70;
  int32_t count = 0;
  char err[256];
  CPMDCRestart *file;
  memset(part1, ' ', sizeof part1);
  memset(part2, ' ', sizeof part2);
  memcpy(part1, " STUTTGART", 10);
  memcpy(part2, " VERSION 3.0", 12);
  if (add_bytes(&buf, &len, &cap, &m1, 4) || add_bytes(&buf, &len, &cap, part1, 10) ||
      add_bytes(&buf, &len, &cap, &m1, 4) || add_bytes(&buf, &len, &cap, &m2, 4) ||
      add_bytes(&buf, &len, &cap, part2, 70) || add_bytes(&buf, &len, &cap, &m2, 4))
    return fail("split header build");
  if (add_rec(&buf, &len, &cap, &count, 4) != 0)
    return fail("split trailer");
  file = cpmdc_restart_read_mem(buf, len, err, sizeof err);
  free(buf);
  if (!file)
    return fail(err[0] ? err : "split header was rejected");
  if (!strstr(cpmdc_restart_header(file), "STUTTGART VERSION 3.0")) {
    cpmdc_restart_free(file);
    return fail("split header text");
  }
  /* 0.4.0 callers spell the type cpmdc_restart; it must stay the same type. */
  cpmdc_restart *old_name = file;
  cpmdc_restart_free(old_name);
  return 0;
}

static int check_bad_magic(void) {
  unsigned char *buf = NULL;
  size_t len = 0, cap = 0;
  char header[80];
  char err[256];
  memset(header, ' ', sizeof header);
  memcpy(header, " NOT A RESTART", 14);
  if (add_rec(&buf, &len, &cap, header, 80) != 0)
    return fail("bad magic build");
  if (cpmdc_restart_read_mem(buf, len, err, sizeof err) != NULL) {
    free(buf);
    return fail("bad magic was accepted");
  }
  free(buf);
  if (!strstr(err, "STUTTGART"))
    return fail("bad magic error text");
  return 0;
}

static int check_geo1(const char *path) {
  char err[256];
  CPMDCRestart *file;
  double cell[6];
  double taui[21];
  const int *na = NULL;
  int nsp = 0, ibrav = 0, count6 = 0, sec4 = -1;
  int nel = 0, n1 = 0, dual = -1, nstates = 0, nk = -1;
  double ecut = 0, cdual = 0;
  err[0] = '\0';
  file = cpmdc_restart_read_path(path, err, sizeof err);
  if (!file)
    return fail(err[0] ? err : "geo1 read failed");
  if (cpmdc_restart_is_stream(file)) {
    cpmdc_restart_free(file);
    return fail("geo1 header is not marked STREAM");
  }
  if (cpmdc_restart_species(file, &nsp, &na) != 0 || nsp != 2 || na[0] != 3 || na[1] != 4) {
    cpmdc_restart_free(file);
    return fail("geo1 species");
  }
  if (cpmdc_restart_section(file, 4, &sec4) != 0 || sec4 != 0 || cpmdc_restart_ncoords(file) != 0) {
    cpmdc_restart_free(file);
    return fail("geo1 section 4 is empty");
  }
  if (cpmdc_restart_section(file, 6, &count6) != 0 || count6 != 7) {
    cpmdc_restart_free(file);
    return fail("geo1 initial-coordinate count");
  }
  if (cpmdc_restart_initial_coordinates(file, taui, 21) != 0) {
    cpmdc_restart_free(file);
    return fail("geo1 initial coordinates");
  }
  if (cpmdc_restart_cell(file, &ibrav, NULL, cell) != 0 || ibrav != 1 || cell[0] < 10.0) {
    cpmdc_restart_free(file);
    return fail("geo1 cell");
  }
  if (cpmdc_restart_cutoff(file, &ecut, &cdual, &dual, &nel, &n1, NULL, NULL) != 0 || ecut != 70.0 ||
      cdual != 4.0 || dual != 0 || nel != 32 || n1 != 144) {
    cpmdc_restart_free(file);
    return fail("geo1 cutoff");
  }
  if (cpmdc_restart_states(file, &nstates, &nk, NULL, NULL, NULL, NULL) != 0 || nstates != 16 ||
      nk != 0) {
    cpmdc_restart_free(file);
    return fail("geo1 states");
  }
  cpmdc_restart_free(file);
  return 0;
}

/* The docs examples run cpmdc-restart on this sample: three atoms with
 * coordinates and velocities, which the geo1 fixture does not carry. */
static int write_sample(const char *path) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("build");
  FILE *fp = fopen(path, "wb");
  if (fp == NULL) {
    free(bytes);
    return fail("open sample");
  }
  size_t wrote = fwrite(bytes, 1, n, fp);
  free(bytes);
  if (fclose(fp) != 0 || wrote != n)
    return fail("write sample");
  return 0;
}

/* TMPDIR first, as POSIX asks; P_tmpdir is the C library's own default. */
static const char *temp_root(void) {
  const char *root = getenv("TMPDIR");
  if (root && root[0])
    return root;
#ifdef P_tmpdir
  return P_tmpdir;
#else
  return ".";
#endif
}

static int read_whole(const char *path, unsigned char **out, size_t *len) {
  FILE *fp = fopen(path, "rb");
  long size;
  unsigned char *buf;
  if (!fp)
    return -1;
  if (fseek(fp, 0, SEEK_END) != 0 || (size = ftell(fp)) < 0 || fseek(fp, 0, SEEK_SET) != 0) {
    fclose(fp);
    return -1;
  }
  buf = malloc(size ? (size_t)size : 1u);
  if (!buf || fread(buf, 1, (size_t)size, fp) != (size_t)size) {
    free(buf);
    fclose(fp);
    return -1;
  }
  fclose(fp);
  *out = buf;
  *len = (size_t)size;
  return 0;
}

/* Every name in `dir` other than `.`, `..` and `keep` is a leftover. */
static int leftovers(const char *dir, const char *keep, int *kept) {
  DIR *listing = opendir(dir);
  struct dirent *ent;
  int extra = 0;
  *kept = 0;
  if (!listing)
    return -1;
  while ((ent = readdir(listing)) != NULL) {
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
      continue;
    if (keep && strcmp(ent->d_name, keep) == 0)
      (*kept)++;
    else
      extra++;
  }
  closedir(listing);
  return extra;
}

static int atomic_case(const char *dir) {
  unsigned char *bytes = NULL;
  unsigned char *disk = NULL;
  void *rewritten = NULL;
  size_t n = 0, have = 0, rn = 0;
  char path[PATH_MAX];
  char err[256];
  CPMDCRestart *file;
  CPMDCRestart *back;
  FILE *marker;
  int kept = 0;
  if (snprintf(path, sizeof path, "%s/RESTART", dir) >= (int)sizeof path)
    return fail("restart path too long");
  marker = fopen(path, "wb");
  if (!marker)
    return fail("seed restart");
  if (fwrite("OLDOLDOLDOLD", 1, 12, marker) != 12) {
    fclose(marker);
    return fail("seed bytes");
  }
  fclose(marker);
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("atomic sample");
  err[0] = '\0';
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  if (!file)
    return fail(err[0] ? err : "atomic sample rejected");
  if (cpmdc_restart_write_path(file, path, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    return fail(err[0] ? err : "atomic write failed");
  }
  if (cpmdc_restart_write_mem(file, &rewritten, &rn, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    return fail("atomic memory image");
  }
  cpmdc_restart_free(file);
  back = cpmdc_restart_read_path(path, err, sizeof err);
  if (!back) {
    free(rewritten);
    return fail(err[0] ? err : "atomic read-back failed");
  }
  cpmdc_restart_free(back);
  if (read_whole(path, &disk, &have) != 0) {
    free(rewritten);
    return fail("read published restart");
  }
  if (have != rn || memcmp(disk, rewritten, rn) != 0) {
    free(disk);
    free(rewritten);
    return fail("published restart is not the finished image");
  }
  free(disk);
  free(rewritten);
  if (leftovers(dir, "RESTART", &kept) != 0 || kept != 1)
    return fail("restart directory does not hold one RESTART and nothing else");

  if (snprintf(path, sizeof path, "%s/missing/RESTART", dir) >= (int)sizeof path)
    return fail("restart path too long");
  bytes = NULL;
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("rebuild for missing directory");
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  if (!file)
    return fail("missing-directory sample");
  err[0] = '\0';
  if (cpmdc_restart_write_path(file, path, err, sizeof err) == 0) {
    cpmdc_restart_free(file);
    return fail("write into a missing directory succeeded");
  }
  cpmdc_restart_free(file);
  if (err[0] == '\0')
    return fail("a failed write left no diagnostic");
  if (access(path, F_OK) == 0)
    return fail("a failed write published a restart");
  return 0;
}

/* A rename onto a directory fails after the bytes are synced. The
 * temporary file goes and the directory stays. */
static int rename_failure_case(const char *dir) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  char sub[PATH_MAX];
  char path[PATH_MAX];
  char err[256];
  CPMDCRestart *file;
  int kept = 0, extra;
  if (snprintf(sub, sizeof sub, "%s/rename", dir) >= (int)sizeof sub ||
      snprintf(path, sizeof path, "%s/RESTART.dir", sub) >= (int)sizeof path)
    return fail("restart path too long");
  if (mkdir(sub, 0700) != 0 || mkdir(path, 0700) != 0)
    return fail("make the blocking directory");
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("rename-failure sample");
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  if (!file)
    return fail("rename-failure sample rejected");
  err[0] = '\0';
  if (cpmdc_restart_write_path(file, path, err, sizeof err) == 0) {
    cpmdc_restart_free(file);
    return fail("a rename onto a directory succeeded");
  }
  cpmdc_restart_free(file);
  if (strstr(err, "cannot replace the RESTART file") == NULL)
    return fail(err[0] ? err : "a failed rename left no diagnostic");
  extra = leftovers(sub, "RESTART.dir", &kept);
  rmdir(path);
  if (rmdir(sub) != 0 || extra != 0 || kept != 1)
    return fail("a failed rename left its temporary file behind");
  return 0;
}

/* An existing file keeps its permissions; only a new file takes 0666 under
 * the umask, matching fopen(path, "wb"). */
static int mode_case(const char *dir) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  char path[PATH_MAX];
  char err[256];
  CPMDCRestart *file;
  struct stat st;
  FILE *seed;
  mode_t old;
  int rc;
  if (snprintf(path, sizeof path, "%s/RESTART", dir) >= (int)sizeof path)
    return fail("restart path too long");
  seed = fopen(path, "wb");
  if (!seed || fclose(seed) != 0 || chmod(path, 0600) != 0)
    return fail("seed 0600 restart");
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("mode sample");
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  if (!file)
    return fail("mode sample rejected");
  old = umask(022);
  rc = cpmdc_restart_write_path(file, path, err, sizeof err);
  umask(old);
  cpmdc_restart_free(file);
  if (rc != 0)
    return fail(err);
  if (stat(path, &st) != 0 || (st.st_mode & 0777) != 0600)
    return fail("published restart did not preserve mode 0600");
  if (unlink(path) != 0)
    return fail("remove mode-preservation restart");
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("new-file mode sample");
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  if (!file)
    return fail("new-file mode sample rejected");
  old = umask(022);
  rc = cpmdc_restart_write_path(file, path, err, sizeof err);
  umask(old);
  cpmdc_restart_free(file);
  if (rc != 0)
    return fail(err);
  if (stat(path, &st) != 0 || (st.st_mode & 0777) != 0644)
    return fail("published restart does not take 0666 under the umask");
  return 0;
}

/* Publish one sample at `path` under umask 022. `*image` is the bytes that
 * must land on disk. */
static int write_restart_at(const char *path, void **image, size_t *image_n, char *err,
                            size_t err_cap) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  CPMDCRestart *file;
  mode_t old;
  int rc;
  *image = NULL;
  *image_n = 0;
  err[0] = '\0';
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return -1;
  file = cpmdc_restart_read_mem(bytes, n, err, err_cap);
  free(bytes);
  if (!file)
    return -1;
  if (cpmdc_restart_write_mem(file, image, image_n, err, err_cap) != 0) {
    cpmdc_restart_free(file);
    return -1;
  }
  old = umask(022);
  rc = cpmdc_restart_write_path(file, path, err, err_cap);
  umask(old);
  cpmdc_restart_free(file);
  return rc;
}

static int published_matches(const char *path, const void *image, size_t image_n, mode_t mode) {
  unsigned char *disk = NULL;
  size_t have = 0;
  struct stat st;
  int ok;
  if (stat(path, &st) != 0 || (st.st_mode & 0777) != mode)
    return 0;
  if (read_whole(path, &disk, &have) != 0)
    return 0;
  ok = have == image_n && memcmp(disk, image, image_n) == 0;
  free(disk);
  return ok;
}

/* RESTART -> RESTART.1 with RESTART.1 absent. realpath fails here; the
 * target is created beside the link. This case fails on the publisher
 * that stops at that realpath error. */
static int dangling_relative_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], target[PATH_MAX], linkpath[PATH_MAX], err[256];
  char linktext[PATH_MAX];
  struct stat st;
  ssize_t linklen;
  int kept = 0;
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/dangling-rel", dir) >= (int)sizeof sub ||
      snprintf(target, sizeof target, "%s/RESTART.1", sub) >= (int)sizeof target ||
      snprintf(linkpath, sizeof linkpath, "%s/RESTART", sub) >= (int)sizeof linkpath)
    return fail("dangling relative path too long");
  if (mkdir(sub, 0700) != 0 || symlink("RESTART.1", linkpath) != 0)
    return fail("make dangling relative restart link");
  if (write_restart_at(linkpath, &image, &image_n, err, sizeof err) != 0)
    goto done;
  if (lstat(linkpath, &st) != 0 || !S_ISLNK(st.st_mode))
    goto done;
  linklen = readlink(linkpath, linktext, sizeof linktext - 1);
  if (linklen < 0)
    goto done;
  linktext[linklen] = '\0';
  if (strcmp(linktext, "RESTART.1") != 0 || !published_matches(target, image, image_n, 0644))
    goto done;
  if (leftovers(sub, "RESTART", &kept) != 1 || kept != 1)
    goto done;
  rc = 0;
done:
  free(image);
  unlink(linkpath);
  unlink(target);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0 : fail(err[0] ? err : "dangling relative restart was not created");
}

/* An absolute dangling target is created at that path, not beside the
 * link under a new name. This case fails on the publisher that rejects
 * every symlink realpath cannot resolve. */
static int dangling_absolute_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], target[PATH_MAX], linkpath[PATH_MAX], err[256];
  char linktext[PATH_MAX];
  struct stat st;
  ssize_t linklen;
  int kept = 0;
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/dangling-abs", dir) >= (int)sizeof sub ||
      snprintf(target, sizeof target, "%s/RESTART.1", sub) >= (int)sizeof target ||
      snprintf(linkpath, sizeof linkpath, "%s/RESTART", sub) >= (int)sizeof linkpath)
    return fail("dangling absolute path too long");
  if (mkdir(sub, 0700) != 0 || symlink(target, linkpath) != 0)
    return fail("make dangling absolute restart link");
  if (write_restart_at(linkpath, &image, &image_n, err, sizeof err) != 0)
    goto done;
  if (lstat(linkpath, &st) != 0 || !S_ISLNK(st.st_mode))
    goto done;
  linklen = readlink(linkpath, linktext, sizeof linktext - 1);
  if (linklen < 0 || (size_t)linklen >= sizeof linktext)
    goto done;
  linktext[linklen] = '\0';
  if (strcmp(linktext, target) != 0 || !published_matches(target, image, image_n, 0644))
    goto done;
  if (leftovers(sub, "RESTART", &kept) != 1 || kept != 1)
    goto done;
  rc = 0;
done:
  free(image);
  unlink(linkpath);
  unlink(target);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0 : fail(err[0] ? err : "dangling absolute restart was not created");
}

/* True exactly when `path` is a symbolic link whose literal target text
 * is `target`, unchanged by whatever the publisher did at the far end of
 * the chain. */
static int is_symlink_to(const char *path, const char *target) {
  struct stat st;
  char text[PATH_MAX];
  ssize_t n;
  if (lstat(path, &st) != 0 || !S_ISLNK(st.st_mode))
    return 0;
  n = readlink(path, text, sizeof text - 1);
  if (n < 0)
    return 0;
  text[n] = '\0';
  return strcmp(text, target) == 0;
}

/* RESTART -> RESTART.1 -> RESTART.2, with RESTART.2 absent. realpath fails
 * on the whole chain, and the previous fallback read only the first hop,
 * so it replaced RESTART.1 with a new regular file instead of following
 * it to RESTART.2. Both links must stay exactly as they were, and the
 * published file must land at the chain's actual end. */
static int dangling_chain_2hop_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], link1[PATH_MAX], link2[PATH_MAX], target[PATH_MAX], err[256];
  int kept = 0, extra;
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/dangling-chain2", dir) >= (int)sizeof sub ||
      snprintf(link1, sizeof link1, "%s/RESTART", sub) >= (int)sizeof link1 ||
      snprintf(link2, sizeof link2, "%s/RESTART.1", sub) >= (int)sizeof link2 ||
      snprintf(target, sizeof target, "%s/RESTART.2", sub) >= (int)sizeof target)
    return fail("two-hop dangling chain path too long");
  if (mkdir(sub, 0700) != 0 || symlink("RESTART.1", link1) != 0 ||
      symlink("RESTART.2", link2) != 0)
    return fail("make two-hop dangling restart chain");
  if (write_restart_at(link1, &image, &image_n, err, sizeof err) != 0)
    goto done;
  if (!is_symlink_to(link1, "RESTART.1") || !is_symlink_to(link2, "RESTART.2"))
    goto done;
  if (!published_matches(target, image, image_n, 0644))
    goto done;
  extra = leftovers(sub, "RESTART", &kept);
  if (extra != 2 || kept != 1)
    goto done;
  rc = 0;
done:
  free(image);
  unlink(link1);
  unlink(link2);
  unlink(target);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0
                 : fail(err[0] ? err : "two-hop dangling chain was not resolved to its end");
}

/* RESTART -> RESTART.1 -> RESTART.2 -> RESTART.3, with RESTART.3 absent.
 * One more hop than the two-hop case, to show the walk is not bounded to a
 * fixed small depth. All three links survive; only RESTART.3 is new. */
static int dangling_chain_3hop_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], link1[PATH_MAX], link2[PATH_MAX], link3[PATH_MAX];
  char target[PATH_MAX], err[256];
  int kept = 0, extra;
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/dangling-chain3", dir) >= (int)sizeof sub ||
      snprintf(link1, sizeof link1, "%s/RESTART", sub) >= (int)sizeof link1 ||
      snprintf(link2, sizeof link2, "%s/RESTART.1", sub) >= (int)sizeof link2 ||
      snprintf(link3, sizeof link3, "%s/RESTART.2", sub) >= (int)sizeof link3 ||
      snprintf(target, sizeof target, "%s/RESTART.3", sub) >= (int)sizeof target)
    return fail("three-hop dangling chain path too long");
  if (mkdir(sub, 0700) != 0 || symlink("RESTART.1", link1) != 0 ||
      symlink("RESTART.2", link2) != 0 || symlink("RESTART.3", link3) != 0)
    return fail("make three-hop dangling restart chain");
  if (write_restart_at(link1, &image, &image_n, err, sizeof err) != 0)
    goto done;
  if (!is_symlink_to(link1, "RESTART.1") || !is_symlink_to(link2, "RESTART.2") ||
      !is_symlink_to(link3, "RESTART.3"))
    goto done;
  if (!published_matches(target, image, image_n, 0644))
    goto done;
  extra = leftovers(sub, "RESTART", &kept);
  if (extra != 3 || kept != 1)
    goto done;
  rc = 0;
done:
  free(image);
  unlink(link1);
  unlink(link2);
  unlink(link3);
  unlink(target);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0
                 : fail(err[0] ? err : "three-hop dangling chain was not resolved to its end");
}

/* RESTART -> RESTART.1 -> RESTART, a loop of two links rather than one
 * link pointing at itself. The write must fail with a clear diagnostic and
 * leave both links exactly as they were. */
static int dangling_chain_loop_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], link1[PATH_MAX], link2[PATH_MAX], err[256];
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/dangling-chain-loop", dir) >= (int)sizeof sub ||
      snprintf(link1, sizeof link1, "%s/RESTART", sub) >= (int)sizeof link1 ||
      snprintf(link2, sizeof link2, "%s/RESTART.1", sub) >= (int)sizeof link2)
    return fail("two-link loop path too long");
  if (mkdir(sub, 0700) != 0 || symlink("RESTART.1", link1) != 0 ||
      symlink("RESTART", link2) != 0)
    return fail("make a two-link restart symlink loop");
  if (write_restart_at(link1, &image, &image_n, err, sizeof err) == 0)
    goto done;
  if (strstr(err, "cannot resolve the RESTART symlink") == NULL)
    goto done;
  if (!is_symlink_to(link1, "RESTART.1") || !is_symlink_to(link2, "RESTART"))
    goto done;
  rc = 0;
done:
  free(image);
  unlink(link1);
  unlink(link2);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0 : fail(err[0] ? err : "two-link symlink loop was not rejected");
}

#ifdef CPMDC_RESTART_TEST_HOOK
/* RESTART -> RESTART.1 -> RESTART.2, with RESTART.2 absent. After realpath
 * reports ENOENT the hook replaces the intermediate symlink RESTART.1 with
 * a regular file. The publisher used to return that file as the dangling
 * destination and rename over it. The write must fail, the file must still
 * hold the swapped bytes, and RESTART must still name RESTART.1. */
static int swapped_intermediate_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], link1[PATH_MAX], link2[PATH_MAX], final[PATH_MAX], err[256];
  char bytes[8];
  struct stat st;
  int fd, kept = 0, extra, rc = 1;
  ssize_t n;
  if (snprintf(sub, sizeof sub, "%s/swap-mid", dir) >= (int)sizeof sub ||
      snprintf(link1, sizeof link1, "%s/RESTART", sub) >= (int)sizeof link1 ||
      snprintf(link2, sizeof link2, "%s/RESTART.1", sub) >= (int)sizeof link2 ||
      snprintf(final, sizeof final, "%s/RESTART.2", sub) >= (int)sizeof final ||
      snprintf(hook_swap_path, sizeof hook_swap_path, "%s", link2) >= (int)sizeof hook_swap_path)
    return fail("swapped intermediate path too long");
  if (mkdir(sub, 0700) != 0 || symlink("RESTART.1", link1) != 0 || symlink("RESTART.2", link2) != 0)
    return fail("make the chain a swap can break");
  hook_swap_armed = 1;
  hook_swap_ok = 0;
  cpmdc_restart_test_before_follow = swap_intermediate_for_file;
  if (write_restart_at(link1, &image, &image_n, err, sizeof err) == 0)
    goto done;
  cpmdc_restart_test_before_follow = NULL;
  hook_swap_armed = 0;
  if (!hook_swap_ok)
    goto done;
  if (strstr(err, "the RESTART symlink chain changed before publication") == NULL)
    goto done;
  if (!is_symlink_to(link1, "RESTART.1"))
    goto done;
  if (lstat(link2, &st) != 0 || !S_ISREG(st.st_mode))
    goto done;
  fd = open(link2, O_RDONLY);
  if (fd < 0)
    goto done;
  n = read(fd, bytes, sizeof bytes);
  close(fd);
  if (n != 4 || memcmp(bytes, "KEEP", 4) != 0)
    goto done;
  if (lstat(final, &st) == 0)
    goto done;
  extra = leftovers(sub, "RESTART", &kept);
  if (extra != 1 || kept != 1)
    goto done;
  rc = 0;
done:
  cpmdc_restart_test_before_follow = NULL;
  hook_swap_armed = 0;
  free(image);
  unlink(link1);
  unlink(link2);
  unlink(final);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0 : fail(err[0] ? err : "a swapped intermediate was not refused");
}
#endif

/* A link that points at itself is not a place to create a file. */
static int symlink_loop_case(const char *dir) {
  void *image = NULL;
  size_t image_n = 0;
  char sub[PATH_MAX], linkpath[PATH_MAX], err[256];
  char linktext[64];
  struct stat st;
  ssize_t linklen;
  int kept = 0;
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/symlink-loop", dir) >= (int)sizeof sub ||
      snprintf(linkpath, sizeof linkpath, "%s/RESTART", sub) >= (int)sizeof linkpath)
    return fail("symlink loop path too long");
  if (mkdir(sub, 0700) != 0 || symlink("RESTART", linkpath) != 0)
    return fail("make restart symlink loop");
  if (write_restart_at(linkpath, &image, &image_n, err, sizeof err) == 0)
    goto done;
  if (strstr(err, "cannot resolve the RESTART symlink") == NULL)
    goto done;
  if (lstat(linkpath, &st) != 0 || !S_ISLNK(st.st_mode))
    goto done;
  linklen = readlink(linkpath, linktext, sizeof linktext - 1);
  if (linklen < 0)
    goto done;
  linktext[linklen] = '\0';
  if (strcmp(linktext, "RESTART") != 0 || leftovers(sub, "RESTART", &kept) != 0 || kept != 1)
    goto done;
  rc = 0;
done:
  free(image);
  unlink(linkpath);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0 : fail(err[0] ? err : "symlink loop was not rejected");
}

static int symlink_case(const char *dir) {
  unsigned char *bytes = NULL;
  unsigned char *disk = NULL;
  void *rewritten = NULL;
  size_t n = 0, have = 0, rn = 0;
  char sub[PATH_MAX], target[PATH_MAX], linkpath[PATH_MAX], err[256];
  char linktext[PATH_MAX];
  struct stat st;
  CPMDCRestart *file;
  FILE *seed;
  ssize_t linklen;
  int rc = 1;
  if (snprintf(sub, sizeof sub, "%s/symlink", dir) >= (int)sizeof sub ||
      snprintf(target, sizeof target, "%s/RESTART.1", sub) >= (int)sizeof target ||
      snprintf(linkpath, sizeof linkpath, "%s/RESTART", sub) >= (int)sizeof linkpath)
    return fail("symlink restart path too long");
  if (mkdir(sub, 0700) != 0)
    return fail("make symlink restart directory");
  seed = fopen(target, "wb");
  if (!seed || fwrite("STALE", 1, 5, seed) != 5 || fclose(seed) != 0)
    goto done;
  if (symlink("RESTART.1", linkpath) != 0)
    goto done;
  if (build_sample(&bytes, &n, 0, 1) != 0)
    goto done;
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  bytes = NULL;
  if (!file)
    goto done;
  if (cpmdc_restart_write_mem(file, &rewritten, &rn, err, sizeof err) != 0 ||
      cpmdc_restart_write_path(file, linkpath, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    goto done;
  }
  cpmdc_restart_free(file);
  if (lstat(linkpath, &st) != 0 || !S_ISLNK(st.st_mode))
    goto done;
  linklen = readlink(linkpath, linktext, sizeof linktext - 1);
  if (linklen < 0)
    goto done;
  linktext[linklen] = '\0';
  if (strcmp(linktext, "RESTART.1") != 0 ||
      read_whole(target, &disk, &have) != 0 ||
      have != rn || memcmp(disk, rewritten, rn) != 0)
    goto done;
  rc = 0;
done:
  free(bytes);
  free(disk);
  free(rewritten);
  unlink(linkpath);
  unlink(target);
  if (rmdir(sub) != 0 && rc == 0)
    rc = 1;
  return rc == 0 ? 0 : fail("symlink restart was not updated through its target");
}

static int long_name_case(const char *dir) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  char err[256];
  CPMDCRestart *file = NULL;
  long name_max = pathconf(dir, _PC_NAME_MAX);
  size_t dir_len = strlen(dir);
  char *path;
  int rc = 1;
  if (name_max < 32)
    name_max = NAME_MAX;
  if ((size_t)name_max + dir_len + 2 > PATH_MAX)
    name_max = (long)(PATH_MAX - dir_len - 2);
  path = malloc(dir_len + (size_t)name_max + 2);
  if (!path)
    return fail("allocate long restart path");
  memcpy(path, dir, dir_len);
  path[dir_len] = '/';
  memset(path + dir_len + 1, 'R', (size_t)name_max);
  path[dir_len + 1 + (size_t)name_max] = '\0';
  if (build_sample(&bytes, &n, 0, 1) != 0)
    goto done;
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  bytes = NULL;
  if (!file || cpmdc_restart_write_path(file, path, err, sizeof err) != 0)
    goto done;
  rc = access(path, F_OK) == 0 ? 0 : 1;
done:
  cpmdc_restart_free(file);
  free(bytes);
  unlink(path);
  free(path);
  return rc == 0 ? 0 : fail("NAME_MAX restart name could not be published");
}

/* A bare name lands in, and syncs, the current directory. */
static int bare_name_case(const char *dir) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  char here[PATH_MAX];
  char err[256];
  CPMDCRestart *file;
  int rc, kept = 0;
  if (!getcwd(here, sizeof here))
    return fail("getcwd");
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("bare-name sample");
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  free(bytes);
  if (!file)
    return fail("bare-name sample rejected");
  if (chdir(dir) != 0) {
    cpmdc_restart_free(file);
    return fail("chdir into the restart directory");
  }
  rc = cpmdc_restart_write_path(file, "RESTART", err, sizeof err);
  cpmdc_restart_free(file);
  if (chdir(here) != 0)
    return fail("chdir back");
  if (rc != 0)
    return fail(err);
  if (leftovers(dir, "RESTART", &kept) != 0 || kept != 1)
    return fail("a bare-name write left more than RESTART");
  return 0;
}

static int check_atomic_replace(void) {
  char dir[PATH_MAX];
  char path[PATH_MAX];
  int rc;
  if (snprintf(dir, sizeof dir, "%s/cpmdc-restart-XXXXXX", temp_root()) >= (int)sizeof dir)
    return fail("temporary directory name too long");
  if (!mkdtemp(dir))
    return fail("restart directory");
  rc = atomic_case(dir);
  if (rc == 0)
    rc = rename_failure_case(dir);
  if (rc == 0)
    rc = mode_case(dir);
  if (rc == 0)
    rc = symlink_case(dir);
  if (rc == 0)
    rc = dangling_relative_case(dir);
  if (rc == 0)
    rc = dangling_absolute_case(dir);
  if (rc == 0)
    rc = dangling_chain_2hop_case(dir);
  if (rc == 0)
    rc = dangling_chain_3hop_case(dir);
  if (rc == 0)
    rc = dangling_chain_loop_case(dir);
#ifdef CPMDC_RESTART_TEST_HOOK
  if (rc == 0)
    rc = swapped_intermediate_case(dir);
#endif
  if (rc == 0)
    rc = symlink_loop_case(dir);
  if (rc == 0)
    rc = long_name_case(dir);
  if (rc == 0)
    rc = bare_name_case(dir);
  if (snprintf(path, sizeof path, "%s/RESTART", dir) < (int)sizeof path)
    unlink(path);
  if (rmdir(dir) != 0 && rc == 0)
    rc = fail("restart directory not empty after the test");
  return rc;
}

int main(int argc, char **argv) {
  unsigned char *bytes = NULL;
  size_t n = 0;
  if (argc == 3 && strcmp(argv[1], "--write-sample") == 0)
    return write_sample(argv[2]);
  if (build_sample(&bytes, &n, 0, 1) != 0)
    return fail("build");
  if (check_sample(bytes, n, 0) != 0) {
    free(bytes);
    return 1;
  }
  free(bytes);
  bytes = NULL;
  if (build_sample(&bytes, &n, 1, 1) != 0)
    return fail("build stream");
  if (check_sample(bytes, n, 1) != 0) {
    free(bytes);
    return 1;
  }
  free(bytes);
  bytes = NULL;
  if (build_sample(&bytes, &n, 0, 0) != 0)
    return fail("build classic");
  if (check_sample(bytes, n, 0) != 0) {
    free(bytes);
    return 1;
  }
  free(bytes);
  if (check_split_header() != 0)
    return 1;
  if (check_bad_magic() != 0)
    return 1;
  if (check_atomic_replace() != 0)
    return 1;
  if (argc > 1 && check_geo1(argv[1]) != 0)
    return 1;
  return 0;
}
