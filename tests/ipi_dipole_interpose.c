/* Test stand-in for cpmdc_session_last_dipole.

 * CPMDC_TEST_DIPOLE selects the value the i-PI client observes:
 * "once" returns a finite dipole on the first call and reports none after
 * that, "nan" returns a non-finite dipole, and any other non-empty value
 * returns the same finite dipole on every call. An unset variable reports
 * that no dipole was recomputed.
 */
#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int valid;
  double values[3];
} Dipole;

int cpmdc_session_last_dipole(const void *session, Dipole *out) {
  const char *spec;
  static int calls = 0;

  (void)session;
  if (out == NULL)
    return -1;
  memset(out, 0, sizeof(*out));
  spec = getenv("CPMDC_TEST_DIPOLE");
  if (spec == NULL || spec[0] == '\0')
    return -1;
  calls += 1;
  if (strcmp(spec, "once") == 0 && calls > 1)
    return -1;
  out->valid = 1;
  if (strcmp(spec, "nan") == 0) {
    out->values[0] = nan("");
    out->values[1] = 0.0;
    out->values[2] = 0.0;
    return 0;
  }
  out->values[0] = 1.0;
  out->values[1] = -2.5;
  out->values[2] = 0.0;
  return 0;
}
