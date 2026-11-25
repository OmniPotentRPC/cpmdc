/* Preloaded in front of libcpmdc so cpmdc-ipi-c sees a dipole the reference
 * evaluator never computes. CPMDC_TEST_DIPOLE holds three reals that strtod
 * reads, nan and inf included; unset, the evaluation has no dipole. */
#include "cpmdc.h"

#include <stdlib.h>
#include <string.h>

int cpmdc_session_last_dipole(const CPMDCSession *session, CPMDCDipole *out) {
  const char *text = getenv("CPMDC_TEST_DIPOLE");
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!session || !text)
    return -1;
  char *at = (char *)text;
  for (int i = 0; i < 3; ++i) {
    char *stop = NULL;
    out->values[i] = strtod(at, &stop);
    if (stop == at) {
      memset(out, 0, sizeof(*out));
      return -1;
    }
    at = stop;
  }
  out->valid = 1;
  return 0;
}
