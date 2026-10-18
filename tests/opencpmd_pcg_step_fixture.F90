MODULE native_pcg_step_fixture
  USE iso_c_binding
  IMPLICIT NONE
  INTEGER, PARAMETER :: real_8=c_double
CONTAINS
#include "native_scaled_pcg_step.inc"
  FUNCTION native_step(f0,f1,a) RESULT(value) BIND(C)
    REAL(c_double), VALUE :: f0,f1,a
    REAL(c_double) :: value
    value=scaled_pcg_step(f0,f1,a)
  END FUNCTION native_step
END MODULE native_pcg_step_fixture
