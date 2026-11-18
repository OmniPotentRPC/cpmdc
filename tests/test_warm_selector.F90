! SPDX-License-Identifier: MIT
MODULE test_warm_selector_metric
  USE warm_selector, ONLY: real_8
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  IMPLICIT NONE
CONTAINS
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
END MODULE

PROGRAM test_warm_selector
  USE warm_selector
  USE test_warm_selector_metric, ONLY: metric
  IMPLICIT NONE
  REAL(real_8) :: cell(9), geo(3), dist
  INTEGER :: source, nslots
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
  ! Copy duplicates c0 and geometry. An occupied slot is replaced. The
  ! live store and the selected key stay put. A key that never converged
  ! is refused and adds no slot.
  CALL embed_reset_warm_orbitals()
  CALL embed_select_orbital_slot(21_embed_key_kind)
  geo = [.1D0, 0D0, 0D0]
  CALL choose(.TRUE., .5D0)
  orbitals = CMPLX(21D0, 0D0, real_8)
  CALL embed_save_orbitals(orbitals)
  CALL embed_select_orbital_slot(22_embed_key_kind)
  geo = [.4D0, 0D0, 0D0]
  CALL choose(.TRUE., .5D0)
  orbitals = CMPLX(22D0, 0D0, real_8)
  CALL embed_save_orbitals(orbitals)
  CALL embed_select_orbital_slot(23_embed_key_kind)
  geo = [1D0, 0D0, 0D0]
  CALL choose(.TRUE., .5D0)
  orbitals = CMPLX(23D0, 0D0, real_8)
  CALL embed_save_orbitals(orbitals)
  nslots = embed_orbital_slot_count()
  IF (nslots /= 2) ERROR STOP 9
  IF (.NOT. embed_copy_orbital_slot(21_embed_key_kind, 22_embed_key_kind)) &
       ERROR STOP 10
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 11
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(23D0,0D0,real_8)) > 1D-12) ERROR STOP 12
  IF (embed_copy_orbital_slot(99_embed_key_kind, 30_embed_key_kind)) ERROR STOP 13
  IF (embed_orbital_slot_count() /= 2) ERROR STOP 14
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(23D0,0D0,real_8)) > 1D-12) ERROR STOP 15
  CALL embed_select_orbital_slot(22_embed_key_kind)
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(21D0,0D0,real_8)) > 1D-12) ERROR STOP 16
  geo = [.1D0, 0D0, 0D0]
  CALL choose(.TRUE., .05D0)
  IF (source /= 0 .OR. ABS(dist) > 1D-12) ERROR STOP 17
  CALL embed_select_orbital_slot(21_embed_key_kind)
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(21D0,0D0,real_8)) > 1D-12) ERROR STOP 18
  ! The selected key's orbitals live in the store, so a copy onto it
  ! replaces that store and leaves the key selected.
  CALL embed_reset_warm_orbitals()
  CALL embed_select_orbital_slot(31_embed_key_kind)
  geo = [.2D0, 0D0, 0D0]
  CALL choose(.TRUE., .5D0)
  orbitals = CMPLX(31D0, 0D0, real_8)
  CALL embed_save_orbitals(orbitals)
  CALL embed_select_orbital_slot(32_embed_key_kind)
  geo = [.8D0, 0D0, 0D0]
  CALL choose(.TRUE., .5D0)
  orbitals = CMPLX(32D0, 0D0, real_8)
  CALL embed_save_orbitals(orbitals)
  IF (.NOT. embed_copy_orbital_slot(31_embed_key_kind, 32_embed_key_kind)) &
       ERROR STOP 19
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(31D0,0D0,real_8)) > 1D-12) ERROR STOP 20
  CALL embed_select_orbital_slot(32_embed_key_kind)
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(31D0,0D0,real_8)) > 1D-12) ERROR STOP 21
  CALL embed_select_orbital_slot(31_embed_key_kind)
  orbitals = 0
  CALL embed_restore_orbitals(orbitals)
  IF (ABS(orbitals(1,1,1)-CMPLX(31D0,0D0,real_8)) > 1D-12) ERROR STOP 22
CONTAINS
  SUBROUTINE choose(periodic, limit)
    LOGICAL, INTENT(IN) :: periodic
    REAL(real_8), INTENT(IN) :: limit
    CALL embed_choose_warm_start_metric_v2(geo,cell,periodic,limit,dist,source,donor,metric)
  END SUBROUTINE
END PROGRAM
