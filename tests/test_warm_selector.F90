! SPDX-License-Identifier: MIT
PROGRAM test_warm_selector
  USE warm_selector
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  IMPLICIT NONE
  REAL(real_8) :: cell(9), geo(3), dist
  INTEGER :: source, parsed_limit
  LOGICAL :: parsed_valid
  INTEGER(embed_key_kind) :: donor
  COMPLEX(real_8) :: orbitals(1,1,1)
  cell = [4D0,0D0,0D0, 0D0,5D0,0D0, 0D0,0D0,6D0]
  CALL embed_select_orbital_slot(11_embed_key_kind)
  geo = [.1D0,0D0,0D0]
  CALL choose(.TRUE., .5D0)
  IF (source /= -2) ERROR STOP 1
  orbitals = CMPLX(11D0,0D0,real_8)
  CALL embed_save_orbitals(orbitals)
  CALL embed_select_orbital_slot(12_embed_key_kind)
  geo(1) = .4D0
  CALL choose(.TRUE., .5D0)
  orbitals = CMPLX(12D0,0D0,real_8)
  CALL embed_save_orbitals(orbitals)
  CALL embed_select_orbital_slot(13_embed_key_kind)
  geo(1) = 3.9D0
  CALL choose(.TRUE., .5D0)
  IF (source /= 1 .OR. donor /= 11_embed_key_kind) ERROR STOP 2
  IF (ABS(dist-.2D0) > 1D-12) ERROR STOP 3
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-11D0) > 1D-12) ERROR STOP 4
  ! Boundary-model and cell changes invalidate every parked state.
  CALL choose(.FALSE., .5D0)
  IF (source /= -2 .OR. embed_orbital_slot_count() /= 0) ERROR STOP 5
  orbitals = 13
  CALL embed_save_orbitals(orbitals)
  geo(1) = .1D0
  CALL choose(.FALSE., .5D0)
  IF (source /= -1 .OR. ABS(dist-3.8D0) > 1D-12) ERROR STOP 6
  CALL embed_save_orbitals(orbitals)
  cell(1) = 4.01D0
  CALL choose(.FALSE., -1D0)
  IF (source /= -2 .OR. embed_have_orbitals) ERROR STOP 7
  ! A skew-cell closest image is not componentwise fractional wrapping.
  cell(1) = 4D0
  cell(4) = .4D0
  geo = 0
  CALL choose(.TRUE., 4D0)
  orbitals = 15
  CALL embed_save_orbitals(orbitals)
  geo = [-5.375D0,-2.5D0,-3.125D0]
  CALL choose(.TRUE., 3.95D0)
  IF (source /= 0 .OR. ABS(dist-SQRT(15.46625D0)) > 1D-12) ERROR STOP 8
  ! Unset limit keeps all parked slots; positive limit evicts oldest inactive
  ! slots in park order while protecting the active and incoming keys.
  CALL embed_reset_warm_orbitals()
  CALL embed_parse_orbital_slot_limit(parsed_limit, parsed_valid)
  IF (.NOT. parsed_valid .OR. parsed_limit /= 0) ERROR STOP 9
  cell = [20D0,0D0,0D0, 0D0,20D0,0D0, 0D0,0D0,20D0]
  CALL embed_set_orbital_slot_limit(2)
  CALL park_bead(101_embed_key_kind, 0.0D0, 101D0)
  IF (embed_orbital_slot_count() /= 0) ERROR STOP 10
  CALL park_bead(102_embed_key_kind, 2.0D0, 102D0)
  IF (embed_orbital_slot_count() /= 1) ERROR STOP 11
  CALL park_bead(103_embed_key_kind, 4.0D0, 103D0)
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 12
  CALL park_bead(104_embed_key_kind, 6.0D0, 104D0)
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 13
  CALL park_bead(105_embed_key_kind, 8.0D0, 105D0)
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 14
  ! Parked slots are now {103, 104} and active is 105; selecting 103 must
  ! restore 103 without evicting it while parking 105 -> parked {104, 105}.
  CALL embed_select_orbital_slot(103_embed_key_kind)
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 15
  geo = [4.0D0,0D0,0D0]
  CALL choose(.FALSE., 0.1D0)
  IF (source /= 0) ERROR STOP 16
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-103D0) > 1D-12) ERROR STOP 17
  ! Oldest bead 101 was evicted earlier, so revisiting 101 falls back to guess.
  CALL embed_select_orbital_slot(101_embed_key_kind)
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 18
  geo = [0.0D0,0D0,0D0]
  CALL choose(.FALSE., 0.1D0)
  IF (source /= -1) ERROR STOP 19
  ! Shrinking the limit to 1 trims parked slots to 1 immediately, and
  ! ping-ponging between 2 active keys never evicts the incoming key.
  CALL park_bead(201_embed_key_kind, 1.0D0, 201D0)
  CALL park_bead(202_embed_key_kind, 3.0D0, 202D0)
  CALL embed_set_orbital_slot_limit(1)
  IF (embed_orbital_slot_count() /= 1) ERROR STOP 20
  CALL embed_select_orbital_slot(201_embed_key_kind)
  IF (embed_orbital_slot_count() /= 1) ERROR STOP 21
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-201D0) > 1D-12) ERROR STOP 22
  CALL embed_set_orbital_slot_limit(-1)
CONTAINS
  SUBROUTINE park_bead(k, x, val)
    INTEGER(embed_key_kind), INTENT(IN) :: k
    REAL(real_8), INTENT(IN) :: x, val
    CALL embed_select_orbital_slot(k)
    geo = [x, 0D0, 0D0]
    CALL choose(.FALSE., 0.1D0)
    orbitals = CMPLX(val, 0D0, real_8)
    CALL embed_save_orbitals(orbitals)
  END SUBROUTINE
  SUBROUTINE choose(periodic, limit)
    LOGICAL, INTENT(IN) :: periodic
    REAL(real_8), INTENT(IN) :: limit
    CALL embed_choose_warm_start_metric_v2(geo,cell,periodic,limit,dist,source,donor,metric)
  END SUBROUTINE
  FUNCTION metric(a,b,box,periodic) RESULT(d)
    REAL(real_8), INTENT(IN) :: a(:),b(:),box(9)
    LOGICAL, INTENT(IN) :: periodic
    REAL(real_8) :: d
    INTERFACE
      FUNCTION cpmdc_geometry_rmsd(n,a,b,box,periodic) BIND(C) RESULT(d)
        IMPORT c_int,c_double
        INTEGER(c_int), VALUE :: n,periodic
        REAL(c_double), INTENT(IN) :: a(*),b(*),box(*)
        REAL(c_double) :: d
      END FUNCTION
    END INTERFACE
    d=cpmdc_geometry_rmsd(SIZE(a)/3,a,b,box,MERGE(1,0,periodic))
  END FUNCTION
END PROGRAM
