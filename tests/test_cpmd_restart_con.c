#define _GNU_SOURCE
#include "cpmdc_restart_con.h"

#include "readcon-core.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* The original bohr bits do not survive the trip through angstrom.
 * The inverse divide is not exact for every input.
 * RoundTrip preserves the f64 the multiply produced.
 */
static const char k_bohr_note[] =
    "original bohr bits do not survive the angstrom trip; the inverse divide "
    "is not exact for every input; RoundTrip preserves the f64 the multiply produced";

static int fail(const char *msg) {
  fprintf(stderr, "test_cpmd_restart_con: %s\n", msg);
  return 1;
}

static int add_bytes(unsigned char **buf, size_t *len, size_t *cap, const void *src,
                     size_t n) {
  if (*len + n < *len)
    return -1;
  if (*len + n > *cap) {
    size_t ncap = *cap ? *cap * 2 : 256;
    while (ncap < *len + n) {
      if (ncap > (size_t)-1 / 2)
        return -1;
      ncap *= 2;
    }
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

static int add_count(unsigned char **buf, size_t *len, size_t *cap, int32_t count) {
  unsigned char raw[12];
  int64_t fpos = 0;
  memcpy(raw, &count, 4);
  memcpy(raw + 4, &fpos, 8);
  return add_rec(buf, len, cap, raw, 12);
}

struct payload {
  const unsigned char *data;
  size_t len;
};

static int split_records(const unsigned char *buf, size_t n, struct payload **out,
                         int *nout) {
  size_t off = 0;
  struct payload *recs = NULL;
  int nrec = 0;
  int cap = 0;
  *out = NULL;
  *nout = 0;
  while (off < n) {
    unsigned char *payload = NULL;
    size_t plen = 0;
    size_t pcap = 0;
    for (;;) {
      int32_t lead;
      uint32_t slen;
      int cont;
      int32_t trail;
      if (n - off < 4) {
        free(payload);
        goto fail;
      }
      memcpy(&lead, buf + off, 4);
      off += 4;
      cont = lead < 0;
      if (lead == INT32_MIN) {
        free(payload);
        goto fail;
      }
      slen = cont ? (uint32_t)(-lead) : (uint32_t)lead;
      if ((size_t)slen > n - off) {
        free(payload);
        goto fail;
      }
      if (add_bytes(&payload, &plen, &pcap, buf + off, slen) != 0) {
        free(payload);
        goto fail;
      }
      off += slen;
      if (n - off < 4) {
        free(payload);
        goto fail;
      }
      memcpy(&trail, buf + off, 4);
      off += 4;
      if (trail != lead) {
        free(payload);
        goto fail;
      }
      if (!cont)
        break;
    }
    if (nrec == cap) {
      int ncap = cap ? cap * 2 : 32;
      struct payload *grown = realloc(recs, (size_t)ncap * sizeof(*recs));
      if (!grown) {
        free(payload);
        goto fail;
      }
      recs = grown;
      cap = ncap;
    }
    recs[nrec].data = payload;
    recs[nrec].len = plen;
    nrec++;
  }
  *out = recs;
  *nout = nrec;
  return 0;
fail:
  for (int i = 0; i < nrec; i++)
    free((void *)recs[i].data);
  free(recs);
  return -1;
}

static void free_records(struct payload *recs, int nrec) {
  if (!recs)
    return;
  for (int i = 0; i < nrec; i++)
    free((void *)recs[i].data);
  free(recs);
}

static int index_sections(const struct payload *recs, int nrec, int *start, int *count) {
  size_t i;
  for (int s = 0; s <= 99; s++) {
    start[s] = -1;
    count[s] = 0;
  }
  if (nrec < 1 || !memmem(recs[0].data, recs[0].len, "STUTTGART VERSION 3.0", 21))
    return -1;
  start[1] = 0;
  count[1] = 1;
  i = 1;
  for (int s = 2; s <= 99 && i < (size_t)nrec; s++) {
    int32_t raw;
    int64_t nfollow;
    if (recs[i].len < 4)
      return -1;
    memcpy(&raw, recs[i].data, 4);
    if (raw == INT32_MIN)
      return -1;
    nfollow = raw < 0 ? -(int64_t)raw : (int64_t)raw;
    if (nfollow > (int64_t)(nrec - (int)i - 1))
      return -1;
    start[s] = (int)i;
    count[s] = (int)(1 + nfollow);
    i += (size_t)(1 + nfollow);
  }
  return i == (size_t)nrec ? 0 : -1;
}

static int record_is_patched(const int *start, const int *count, int index) {
  if (start[2] >= 0 && index == start[2] + 2)
    return 1;
  if (start[4] >= 0 && index > start[4] && index < start[4] + count[4])
    return 1;
  if (start[5] >= 0 && index > start[5] && index < start[5] + count[5])
    return 1;
  return 0;
}

static int read_file(const char *path, unsigned char **out, size_t *nout) {
  FILE *fp = fopen(path, "rb");
  unsigned char *buf = NULL;
  size_t n = 0;
  size_t cap = 0;
  unsigned char chunk[4096];
  if (!fp)
    return -1;
  for (;;) {
    size_t got = fread(chunk, 1, sizeof chunk, fp);
    if (got && add_bytes(&buf, &n, &cap, chunk, got) != 0) {
      free(buf);
      fclose(fp);
      return -1;
    }
    if (got < sizeof chunk) {
      if (ferror(fp)) {
        free(buf);
        fclose(fp);
        return -1;
      }
      break;
    }
  }
  fclose(fp);
  *out = buf;
  *nout = n;
  return 0;
}

static int write_temp(const unsigned char *bytes, size_t n, const char *pattern, char *path,
                      size_t path_cap) {
  int fd;
  size_t wrote = 0;
  if (strlen(pattern) + 1 > path_cap)
    return -1;
  memcpy(path, pattern, strlen(pattern) + 1);
  fd = mkstemps(path, 4);
  if (fd < 0)
    return -1;
  while (wrote < n) {
    ssize_t w = write(fd, bytes + wrote, n - wrote);
    if (w < 0) {
      close(fd);
      return -1;
    }
    wrote += (size_t)w;
  }
  if (close(fd) != 0)
    return -1;
  return 0;
}

static int bits_differ(const void *a, const void *b, size_t n) {
  return memcmp(a, b, n) != 0;
}

static int build_water(unsigned char **out, size_t *nout, const double xyz[9],
                       const double vel[9], const double cell[6], const int32_t na[2]) {
  unsigned char *buf = NULL;
  size_t len = 0;
  size_t cap = 0;
  char header[80];
  int32_t sym[2];
  int32_t nsp = 2;
  double taui[9] = {9, 9, 9, 8, 8, 8, 7, 7, 7};
  unsigned char cut[36];
  int32_t states[6] = {8, 0, 100, 80, 20, 10};
  char wave_a[] = "WAVEFUNC-COEFF-A";
  char wave_b[] = "WAVEFUNC-COEFF-B";
  double ecut = 70.0;
  double cdual = 4.0;
  int32_t dual_flag = 0;
  int32_t nel = 8;
  int32_t nr1 = 20;
  int32_t nr2 = 20;
  int32_t nr3 = 20;
  memset(header, ' ', sizeof header);
  memcpy(header, " STUTTGART VERSION 3.0", 22);
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
  if (add_count(&buf, &len, &cap, 2) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, sym, 8) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, cell, 48) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 2) != 0)
    return -1;
  if (add_i32(&buf, &len, &cap, nsp) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, na, 8) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 3) != 0)
    return -1;
  for (int a = 0; a < 3; a++) {
    if (add_rec(&buf, &len, &cap, xyz + 3 * a, 24) != 0)
      return -1;
  }
  if (add_count(&buf, &len, &cap, 3) != 0)
    return -1;
  for (int a = 0; a < 3; a++) {
    if (add_rec(&buf, &len, &cap, vel + 3 * a, 24) != 0)
      return -1;
  }
  if (add_count(&buf, &len, &cap, 3) != 0)
    return -1;
  for (int a = 0; a < 3; a++) {
    if (add_rec(&buf, &len, &cap, taui + 3 * a, 24) != 0)
      return -1;
  }
  if (add_count(&buf, &len, &cap, 1) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, cut, 36) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 1) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, states, 24) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, -2) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, wave_a, sizeof wave_a) != 0)
    return -1;
  if (add_rec(&buf, &len, &cap, wave_b, sizeof wave_b) != 0)
    return -1;
  if (add_count(&buf, &len, &cap, 0) != 0)
    return -1;
  *out = buf;
  *nout = len;
  return 0;
}

static int con_matches_multiply(const char *path, const double original[9],
                                const double original_vel[9], const uint64_t expect_ids[3]) {
  struct RKRConFrame *frame = rkr_read_first_frame(path);
  const double *xyz;
  const double *vel;
  uint64_t ids[3];
  uintptr_t n = 0;
  uintptr_t nvel = 0;
  unsigned char *text = NULL;
  size_t text_n = 0;
  if (!frame)
    return fail("read exported CON");
  xyz = rkr_frame_xyz_f64(frame, &n);
  if (!xyz || n != 3) {
    free_rkr_frame(frame);
    return fail("exported CON atom count");
  }
  if (rkr_frame_copy_atom_ids(frame, ids, 3) != RKR_STATUS_SUCCESS) {
    free_rkr_frame(frame);
    return fail("exported CON atom ids");
  }
  for (int i = 0; i < 3; i++) {
    if (ids[i] != expect_ids[i]) {
      free_rkr_frame(frame);
      return fail("atom_id column does not match this two-species water fixture");
    }
    for (int k = 0; k < 3; k++) {
      double produced = cpmdc_restart_bohr_to_angstrom(original[3 * (int)ids[i] + k]);
      if (memcmp(xyz + 3 * i + k, &produced, sizeof produced) != 0) {
        free_rkr_frame(frame);
        return fail("RoundTrip did not preserve the f64 the multiply produced");
      }
    }
  }
  vel = rkr_frame_velocities_f64(frame, &nvel);
  if (!vel || nvel != 3) {
    free_rkr_frame(frame);
    return fail("exported CON is missing velocities");
  }
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      double produced =
          cpmdc_restart_bohr_to_angstrom(original_vel[3 * (int)ids[i] + k]);
      if (memcmp(vel + 3 * i + k, &produced, sizeof produced) != 0) {
        free_rkr_frame(frame);
        return fail("RoundTrip did not preserve the f64 the multiply produced");
      }
    }
  }
  free_rkr_frame(frame);
  if (read_file(path, &text, &text_n) != 0)
    return fail("reread CON");
  if (!memmem(text, text_n, "(amu*angstrom^2/eV)^0.5", 23)) {
    free(text);
    return fail("velocity export is missing the natural time label");
  }
  if (!memmem(text, text_n, "angstrom", 8) || memmem(text, text_n, "WAVEFUNC", 8)) {
    free(text);
    return fail("CON units or a coefficient leaked into the resume file");
  }
  free(text);
  return 0;
}

static int untouched_records_match(const unsigned char *before, size_t nb,
                                   const unsigned char *after, size_t na,
                                   const double expect_cell[6], const double expect_xyz[9],
                                   const double expect_vel[9]) {
  struct payload *left = NULL;
  struct payload *right = NULL;
  int nl = 0;
  int nr = 0;
  int start[100];
  int count[100];
  int xyz_i = 0;
  int vel_i = 0;
  if (split_records(before, nb, &left, &nl) != 0 || split_records(after, na, &right, &nr) != 0) {
    free_records(left, nl);
    free_records(right, nr);
    return fail("record split");
  }
  if (nl != nr || index_sections(left, nl, start, count) != 0) {
    free_records(left, nl);
    free_records(right, nr);
    return fail("section index");
  }
  for (int i = 0; i < nl; i++) {
    if (left[i].len != right[i].len) {
      free_records(left, nl);
      free_records(right, nr);
      return fail("record length changed");
    }
    if (!record_is_patched(start, count, i)) {
      if (memcmp(left[i].data, right[i].data, left[i].len) != 0) {
        free_records(left, nl);
        free_records(right, nr);
        return fail("a non-coordinate, non-velocity, non-cell record changed");
      }
      continue;
    }
    if (start[2] >= 0 && i == start[2] + 2) {
      if (left[i].len != 48 || memcmp(right[i].data, expect_cell, 48) != 0) {
        free_records(left, nl);
        free_records(right, nr);
        return fail("cell record is not the converted celldm");
      }
    } else if (start[4] >= 0 && i > start[4] && i < start[4] + count[4]) {
      if (left[i].len != 24 || memcmp(right[i].data, expect_xyz + 3 * xyz_i, 24) != 0) {
        free_records(left, nl);
        free_records(right, nr);
        return fail("coordinate record is not the scattered bohr value");
      }
      xyz_i++;
    } else if (start[5] >= 0 && i > start[5] && i < start[5] + count[5]) {
      if (left[i].len != 24 || memcmp(right[i].data, expect_vel + 3 * vel_i, 24) != 0) {
        free_records(left, nl);
        free_records(right, nr);
        return fail("velocity record is not the scattered bohr value");
      }
      vel_i++;
    }
  }
  free_records(left, nl);
  free_records(right, nr);
  if (xyz_i != 3 || vel_i != 3)
    return fail("patched record count");
  if (!memmem(after, na, "WAVEFUNC-COEFF-A", 17) || !memmem(after, na, "WAVEFUNC-COEFF-B", 17))
    return fail("coefficient records were dropped");
  return 0;
}

static int round_trip(const double original[9], double out[9]) {
  for (int i = 0; i < 9; i++)
    out[i] = cpmdc_restart_angstrom_to_bohr(cpmdc_restart_bohr_to_angstrom(original[i]));
  return 0;
}

static int check_fixture(const unsigned char *bytes, size_t n, const char *const *symbols,
                         const uint64_t expect_ids[3], const double xyz[9], const double vel[9],
                         const double cell[6], const char *golden) {
  char err[256];
  char con1[] = "/tmp/cpmdc-water-XXXXXX.con";
  char con2[] = "/tmp/cpmdc-water-XXXXXX.con";
  char restart_path[] = "/tmp/cpmdc-restart-XXXXXX.bin";
  int con1_fd;
  int con2_fd;
  cpmdc_restart *file;
  cpmdc_restart *again;
  unsigned char *first = NULL;
  unsigned char *second = NULL;
  unsigned char *rewritten = NULL;
  unsigned char *golden_bytes = NULL;
  size_t first_n = 0;
  size_t second_n = 0;
  size_t rewritten_n = 0;
  size_t golden_n = 0;
  double round_xyz[9];
  double round_vel[9];
  double lengths[3];
  double angles[3];
  double round_cell[6];
  double got_cell[6];
  double got_xyz[9];
  double got_taui[9];
  void *image = NULL;
  size_t image_n = 0;
  err[0] = '\0';
  con1_fd = mkstemps(con1, 4);
  con2_fd = mkstemps(con2, 4);
  if (con1_fd < 0 || con2_fd < 0)
    return fail("temp CON");
  close(con1_fd);
  close(con2_fd);
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  if (!file)
    return fail(err[0] ? err : "read water");
  if (cpmdc_restart_con_export(file, con1, symbols, 3, err, sizeof err) != 0 ||
      cpmdc_restart_con_export(file, con2, symbols, 3, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    return fail(err[0] ? err : "export");
  }
  cpmdc_restart_free(file);
  if (read_file(con1, &first, &first_n) != 0 || read_file(con2, &second, &second_n) != 0)
    return fail("read export");
  if (first_n != second_n || memcmp(first, second, first_n) != 0) {
    free(first);
    free(second);
    return fail("a second export does not match the first CON");
  }
  free(second);
  if (golden) {
    if (read_file(golden, &golden_bytes, &golden_n) != 0) {
      FILE *fp = fopen("/tmp/cpmdc-restart-water.con", "wb");
      if (fp) {
        fwrite(first, 1, first_n, fp);
        fclose(fp);
      }
      fprintf(stderr, "golden missing, export bytes %zu path %s\n", first_n, con1);
      free(first);
      return fail("golden CON is missing");
    }
    if (golden_n != first_n || memcmp(golden_bytes, first, first_n) != 0) {
      FILE *fp = fopen("/tmp/cpmdc-restart-water.con", "wb");
      if (fp) {
        fwrite(first, 1, first_n, fp);
        fclose(fp);
      }
      fprintf(stderr, "golden bytes %zu export bytes %zu path %s\n", golden_n, first_n, con1);
      free(golden_bytes);
      free(first);
      return fail("export is not byte-exact against the golden CON");
    }
    free(golden_bytes);
  }
  free(first);
  if (con_matches_multiply(con1, xyz, vel, expect_ids) != 0)
    return 1;
  round_trip(xyz, round_xyz);
  round_trip(vel, round_vel);
  if (!bits_differ(xyz, round_xyz, sizeof round_xyz) &&
      !bits_differ(vel, round_vel, sizeof round_vel))
    return fail(k_bohr_note);
  if (cpmdc_restart_celldm_to_con(cell, lengths, angles) != 0 ||
      cpmdc_restart_con_to_celldm(lengths, angles, round_cell) != 0)
    return fail("cell conversion");
  file = cpmdc_restart_read_mem(bytes, n, err, sizeof err);
  if (!file)
    return fail("reread water");
  if (cpmdc_restart_con_import(file, con1, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    return fail(err[0] ? err : "import");
  }
  if (cpmdc_restart_coordinates(file, got_xyz, 9) != 0 ||
      memcmp(got_xyz, round_xyz, sizeof round_xyz) != 0) {
    cpmdc_restart_free(file);
    return fail("imported coordinates are not the scattered divide");
  }
  if (cpmdc_restart_cell(file, NULL, NULL, got_cell) != 0 ||
      memcmp(got_cell, round_cell, sizeof round_cell) != 0) {
    cpmdc_restart_free(file);
    return fail("imported cell is not cpmdc_restart_set_cell of the converted celldm");
  }
  if (cpmdc_restart_initial_coordinates(file, got_taui, 9) != 0 || got_taui[0] != 9.0 ||
      got_taui[8] != 7.0) {
    cpmdc_restart_free(file);
    return fail("initial coordinates changed");
  }
  if (cpmdc_restart_write_mem(file, &image, &image_n, err, sizeof err) != 0) {
    cpmdc_restart_free(file);
    return fail("write imported image");
  }
  cpmdc_restart_free(file);
  if (untouched_records_match(bytes, n, image, image_n, round_cell, round_xyz, round_vel) != 0) {
    free(image);
    return 1;
  }
  free(image);
  if (write_temp(bytes, n, "/tmp/cpmdc-restart-XXXXXX.bin", restart_path,
                 sizeof restart_path) != 0)
    return fail("temp restart");
  again = cpmdc_restart_read_path(restart_path, err, sizeof err);
  if (!again)
    return fail("path read");
  if (cpmdc_restart_con_import(again, con1, err, sizeof err) != 0 ||
      cpmdc_restart_write_path(again, restart_path, err, sizeof err) != 0) {
    cpmdc_restart_free(again);
    return fail(err[0] ? err : "import write");
  }
  cpmdc_restart_free(again);
  if (read_file(restart_path, &rewritten, &rewritten_n) != 0)
    return fail("reread imported restart");
  if (untouched_records_match(bytes, n, rewritten, rewritten_n, round_cell, round_xyz,
                              round_vel) != 0) {
    free(rewritten);
    return 1;
  }
  free(rewritten);
  unlink(con1);
  unlink(con2);
  unlink(restart_path);
  return 0;
}

static int check_cli(const char *tool, const unsigned char *bytes, size_t n,
                     const char *library_con) {
  char restart_path[] = "/tmp/cpmdc-cli-XXXXXX.bin";
  char out_con[] = "/tmp/cpmdc-cli-XXXXXX.con";
  char out_restart[] = "/tmp/cpmdc-cli-XXXXXX.bin";
  char cmd[1024];
  unsigned char *cli = NULL;
  unsigned char *lib = NULL;
  size_t cli_n = 0;
  size_t lib_n = 0;
  int status;
  int restart_fd = mkstemps(restart_path, 4);
  int con_fd = mkstemps(out_con, 4);
  int out_fd = mkstemps(out_restart, 4);
  if (restart_fd < 0 || con_fd < 0 || out_fd < 0)
    return fail("cli temp");
  close(restart_fd);
  close(con_fd);
  close(out_fd);
  {
    FILE *fp = fopen(restart_path, "wb");
    if (!fp || fwrite(bytes, 1, n, fp) != n) {
      if (fp)
        fclose(fp);
      return fail("cli restart write");
    }
    fclose(fp);
  }
  snprintf(cmd, sizeof cmd, "'%s' export-con '%s' '%s' --species O,H", tool, restart_path,
           out_con);
  status = system(cmd);
  if (status != 0)
    return fail("export-con command failed");
  if (read_file(out_con, &cli, &cli_n) != 0 || read_file(library_con, &lib, &lib_n) != 0)
    return fail("cli con read");
  if (cli_n != lib_n || memcmp(cli, lib, cli_n) != 0) {
    free(cli);
    free(lib);
    return fail("export-con bytes differ from the library writer");
  }
  free(cli);
  free(lib);
  snprintf(cmd, sizeof cmd, "'%s' import-con '%s' '%s' '%s'", tool, restart_path, out_con,
           out_restart);
  status = system(cmd);
  if (status != 0)
    return fail("import-con command failed");
  unlink(restart_path);
  unlink(out_con);
  unlink(out_restart);
  return 0;
}

int main(int argc, char **argv) {
  unsigned char *mixed = NULL;
  unsigned char *ordered = NULL;
  size_t mixed_n = 0;
  size_t ordered_n = 0;
  const char *golden = argc > 1 ? argv[1] : NULL;
  const char *tool = argc > 2 ? argv[2] : NULL;
  /* Water stored H, O, H. build groups the two H atoms, so file order is H, H, O. */
  const double mixed_xyz[9] = {0.123456789012345, 0.20, 0.30, 1.11, 2.22,
                               3.33,              4.11, 5.22, 6.33};
  const double mixed_vel[9] = {0.0077, -0.002, 0.0011, 0.013, 0.0, 0.0003, -0.008, 0.0015, 0.004};
  const double cell[6] = {10.123456789012345, 1.0, 1.25, 0.0, 0.0, 0.5};
  const int32_t mixed_na[2] = {2, 1};
  const char *mixed_symbols[] = {"H", "O", "H"};
  const uint64_t mixed_ids[3] = {0, 2, 1};
  const double ordered_xyz[9] = {0.0, 0.0, 0.123456789012345, 1.4301, -0.9876543210987654,
                                 0.02, -1.1, 0.204, 0.3333333333333333};
  const double ordered_vel[9] = {0.0, 0.0, 0.0077, 0.013, -0.002, 0.0011, -0.008, 3.857142857142857,
                                 0.0003};
  const int32_t ordered_na[2] = {1, 2};
  const char *ordered_species[] = {"O", "H"};
  const char *ordered_symbols[3] = {NULL, NULL, NULL};
  const uint64_t ordered_ids[3] = {0, 1, 2};
  char err[256];
  char ordered_con[] = "/tmp/cpmdc-ordered-XXXXXX.con";
  cpmdc_restart *ordered_file;
  if (!golden)
    return fail("golden path");
  if (build_water(&mixed, &mixed_n, mixed_xyz, mixed_vel, cell, mixed_na) != 0)
    return fail("build mixed water");
  if (check_fixture(mixed, mixed_n, mixed_symbols, mixed_ids, mixed_xyz, mixed_vel, cell,
                    golden) != 0) {
    free(mixed);
    return 1;
  }
  if (build_water(&ordered, &ordered_n, ordered_xyz, ordered_vel, cell, ordered_na) != 0) {
    free(mixed);
    return fail("build ordered water");
  }
  ordered_file = cpmdc_restart_read_mem(ordered, ordered_n, err, sizeof err);
  if (!ordered_file ||
      cpmdc_restart_expand_species_symbols(ordered_file, ordered_species, 2, ordered_symbols,
                                           3) != 0) {
    cpmdc_restart_free(ordered_file);
    free(mixed);
    free(ordered);
    return fail("expand O,H");
  }
  if (check_fixture(ordered, ordered_n, ordered_symbols, ordered_ids, ordered_xyz, ordered_vel,
                    cell, NULL) != 0) {
    cpmdc_restart_free(ordered_file);
    free(mixed);
    free(ordered);
    return 1;
  }
  {
    int ordered_fd = mkstemps(ordered_con, 4);
    if (ordered_fd < 0) {
      cpmdc_restart_free(ordered_file);
      free(mixed);
      free(ordered);
      return fail("ordered CON");
    }
    close(ordered_fd);
  }
  if (cpmdc_restart_con_export(ordered_file, ordered_con, ordered_symbols, 3, err, sizeof err) !=
      0) {
    cpmdc_restart_free(ordered_file);
    free(mixed);
    free(ordered);
    return fail(err[0] ? err : "ordered export");
  }
  cpmdc_restart_free(ordered_file);
  if (tool && check_cli(tool, ordered, ordered_n, ordered_con) != 0) {
    free(mixed);
    free(ordered);
    return 1;
  }
  free(mixed);
  free(ordered);
  unlink(ordered_con);
  return 0;
}
