/* cpmdc_session_set_constant_cutoff and the CONSTANT CUTOFF override of the
 * deck renderer. The line goes into the &SYSTEM section once, as A, sigma,
 * and E0 in Ry, with A = 50 x the deck cutoff and sigma = E0 / 10 when they
 * are not given. A request that cannot be placed fails with a reason: E0 +
 * 2 sigma above the deck cutoff, a raw &SYSTEM block, or a message that
 * already sets system.constantCutoff. */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_failed = 0;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n);
  if (buf && fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    buf = NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static void check(int ok, const char *what, const char *detail) {
  if (ok)
    return;
  fprintf(stderr, "FAIL %s\n%s\n", what, detail ? detail : "");
  g_failed = 1;
}

static int count_token(const char *deck, const char *token) {
  int n = 0;
  size_t len = strlen(token);
  for (const char *p = deck; (p = strstr(p, token)) != NULL; p += len)
    n++;
  return n;
}

/* The CONSTANT CUTOFF line lies between &SYSTEM and the &END after it. */
static int in_system_section(const char *deck, const char *needle) {
  const char *sys = strstr(deck, "&SYSTEM");
  const char *hit = strstr(deck, needle);
  if (!sys || !hit || hit < sys)
    return 0;
  const char *end = strstr(sys, "&END");
  return end && hit < end;
}

struct Message {
  unsigned char *bytes;
  size_t size;
  struct capn arena;
  CPMDParams_ptr root;
};

static int open_message(const char *path, struct Message *msg) {
  msg->bytes = read_file(path, &msg->size);
  if (!msg->bytes) {
    fprintf(stderr, "failed to read %s\n", path);
    return -1;
  }
  if (cpmdc_params_root(msg->bytes, msg->size, &msg->arena, &msg->root) != 0) {
    fprintf(stderr, "failed to decode %s\n", path);
    free(msg->bytes);
    return -1;
  }
  return 0;
}

static void close_message(struct Message *msg) {
  cpmdc_params_release(&msg->arena);
  free(msg->bytes);
}

static int render(const struct Message *msg, double cutoff, double height,
                  double width, char *deck, size_t cap) {
  CPMDCScalarOverrides ov;
  memset(&ov, 0, sizeof(ov));
  ov.constant_cutoff_ry = cutoff;
  ov.constant_cutoff_height_ry = height;
  ov.constant_cutoff_width_ry = width;
  deck[0] = '\0';
  return cpmdc_params_render_input_deck_ov(msg->root, &ov, deck, cap);
}

static void structured_system(const struct Message *msg) {
  char deck[CPMDC_BLOCKS];

  check(render(msg, 0.0, 0.0, 0.0, deck, sizeof(deck)) == 0 &&
            strstr(deck, "CONSTANT CUTOFF") == NULL,
        "no request leaves the deck as it is", deck);

  check(render(msg, 30.0, 0.0, 0.0, deck, sizeof(deck)) == 0, "defaults",
        cpmdc_params_render_reason());
  check(in_system_section(deck, " CONSTANT CUTOFF\n  2000 3 30\n"),
        "defaults are 50 x the 40 Ry cutoff and E0 / 10", deck);
  check(count_token(deck, "CONSTANT CUTOFF") == 1, "one line", deck);

  check(render(msg, 30.0, 500.0, 2.0, deck, sizeof(deck)) == 0, "explicit",
        cpmdc_params_render_reason());
  check(in_system_section(deck, " CONSTANT CUTOFF\n  500 2 30\n"),
        "explicit height and width", deck);

  check(render(msg, 34.0, 0.0, 3.0, deck, sizeof(deck)) == 0,
        "E0 + 2 sigma at the cutoff", cpmdc_params_render_reason());

  check(render(msg, 38.0, 0.0, 0.0, deck, sizeof(deck)) != 0,
        "E0 + 2 sigma above the cutoff fails", deck);
  check(strstr(cpmdc_params_render_reason(),
               "exceeds the plane-wave cutoff 40 Ry") != NULL,
        "the reason names the cutoff", cpmdc_params_render_reason());

  double pos[9] = {0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572,
                   -0.4692};
  int z[3] = {8, 1, 1};
  double cube[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  CPMDCScalarOverrides ov;
  memset(&ov, 0, sizeof(ov));
  ov.constant_cutoff_ry = 30.0;
  check(cpmdc_params_render_deck_with_geometry_ov(msg->root, &ov, 3, pos, z,
                                                  cube, 1, deck,
                                                  sizeof(deck)) == 0,
        "geometry render", cpmdc_params_render_reason());
  check(in_system_section(deck, " CONSTANT CUTOFF\n  2000 3 30\n") &&
            count_token(deck, "CONSTANT CUTOFF") == 1,
        "the geometry render carries the line once", deck);
}

static void generated_system(const struct Message *msg) {
  char deck[CPMDC_BLOCKS];
  check(render(msg, 50.0, 0.0, 0.0, deck, sizeof(deck)) == 0,
        "generated &SYSTEM", cpmdc_params_render_reason());
  check(in_system_section(deck, " CONSTANT CUTOFF\n  3500 5 50\n"),
        "the generated &SYSTEM takes the 70 Ry defaults", deck);
}

static void raw_system(const struct Message *msg) {
  char deck[CPMDC_BLOCKS];
  check(render(msg, 50.0, 0.0, 0.0, deck, sizeof(deck)) != 0,
        "a raw &SYSTEM block fails", deck);
  check(strstr(cpmdc_params_render_reason(), "raw text") != NULL,
        "the reason names the raw block", cpmdc_params_render_reason());
}

static void listed_system(const struct Message *msg) {
  char deck[CPMDC_BLOCKS];
  check(render(msg, 0.0, 0.0, 0.0, deck, sizeof(deck)) == 0 &&
            in_system_section(deck, " CONSTANT CUTOFF\n") &&
            count_token(deck, "CONSTANT CUTOFF") == 1,
        "system.constantCutoff renders alone", deck);
  check(render(msg, 50.0, 0.0, 0.0, deck, sizeof(deck)) != 0,
        "a request over system.constantCutoff fails", deck);
  check(strstr(cpmdc_params_render_reason(), "system.constantCutoff") != NULL,
        "the reason names the field", cpmdc_params_render_reason());
}

static void session_abi(const struct Message *msg) {
  check(cpmdc_session_set_constant_cutoff(NULL, 30.0, 0.0, 0.0) == -1,
        "a NULL session fails", cpmdc_last_error());
  CPMDCSession *session = cpmdc_session_create(msg->bytes, msg->size);
  check(session != NULL, "session_create", cpmdc_last_error());
  if (!session)
    return;
  check(cpmdc_session_set_constant_cutoff(session, 30.0, 0.0, 0.0) == 0,
        "set 30 Ry", cpmdc_last_error());
  check(cpmdc_session_set_constant_cutoff(session, NAN, 0.0, 0.0) == -1,
        "a NaN fails", cpmdc_last_error());
  check(cpmdc_session_set_constant_cutoff(session, 38.0, 0.0, 0.0) == -1,
        "38 Ry fails on the 40 Ry deck", cpmdc_last_error());
  check(strstr(cpmdc_last_error(), "exceeds the plane-wave cutoff") != NULL,
        "last_error gives the reason", cpmdc_last_error());
  check(cpmdc_session_set_constant_cutoff(session, 0.0, 0.0, 0.0) == 0,
        "0 turns it off", cpmdc_last_error());
  cpmdc_session_destroy(session);
}

int main(int argc, char **argv) {
  if (argc < 5) {
    fprintf(stderr,
            "usage: %s water_periodic.bin scalar.bin raw_periodic.bin "
            "system_controls.bin\n",
            argv[0]);
    return 2;
  }
  struct Message msg;
  if (open_message(argv[1], &msg) != 0)
    return 1;
  structured_system(&msg);
  session_abi(&msg);
  close_message(&msg);

  if (open_message(argv[2], &msg) != 0)
    return 1;
  generated_system(&msg);
  close_message(&msg);

  if (open_message(argv[3], &msg) != 0)
    return 1;
  raw_system(&msg);
  close_message(&msg);

  if (open_message(argv[4], &msg) != 0)
    return 1;
  listed_system(&msg);
  close_message(&msg);

  cpmdc_finalize();
  if (g_failed)
    return 1;
  printf("constant cutoff deck: ok\n");
  return 0;
}
