SUBROUTINE cpmdc_test_park_partial_orbitals(next_key) &
    BIND(C, NAME='cpmdc_test_park_partial_orbitals')
  USE, INTRINSIC :: iso_c_binding, ONLY: c_long_long
  USE rwfopt_utils, ONLY: embed_park_partial_orbital_slot
  IMPLICIT NONE
  INTEGER(c_long_long), INTENT(IN), VALUE :: next_key

  CALL embed_park_partial_orbital_slot(next_key)
END SUBROUTINE
