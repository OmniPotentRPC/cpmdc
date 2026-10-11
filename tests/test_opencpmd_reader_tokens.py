#!/usr/bin/env python3
"""Extractor checks for reader-scoped OpenCPMD keyword snapshots."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from opencpmd_reader_tokens import extract_tree, tested_literals  # noqa: E402


VDW = """\
MODULE vdwin_utils
CONTAINS
  SUBROUTINE vdwin
    INTEGER :: ierr, iunit
    CHARACTER(len=80) :: line
    ierr = inscan(iunit,'&VDW')
    ! keyword_contains(line,'ONLYCOMMENT')
    IF (keyword_contains(line,'EMPIRICAL',and='CORRECTION')) THEN
       CALL empvdwin(iunit)
    ENDIF
  END SUBROUTINE vdwin

  SUBROUTINE empvdwin(iunit)
    INTEGER :: iunit
    CHARACTER(len=80) :: line
    IF (keyword_contains(line,'S6GRIMME')) THEN
       CONTINUE
    ENDIF
    WRITE(output_unit,*) 'GRIMME is not a keyword here'
  END SUBROUTINE empvdwin
END MODULE vdwin_utils
"""

CPMD = """\
MODULE control_utils
CONTAINS
  SUBROUTINE read_cpmd
    INTEGER :: ierr, iunit
    CHARACTER(len=80) :: line
    ierr = inscan(iunit,'&CPMD')
    IF (keyword_contains(line,'LINALG',and='NEWORTHO')) THEN
       CONTINUE
    ENDIF
  END SUBROUTINE read_cpmd
END MODULE control_utils
"""

VECTORS = """\
MODULE cplngs_utils
CONTAINS
  SUBROUTINE cplrdvec
    INTEGER :: ierr, iunit
    CHARACTER(len=80) :: line
    ierr = inscan(iunit,'&VECTORS')
    ! NEWORTHO belongs to &CPMD, not to this reader
    tfdens = (INDEX(line,'FIXED').NE.0)
  END SUBROUTINE cplrdvec
END MODULE cplngs_utils
"""


class ReaderTokenTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        src = root / "src"
        src.mkdir()
        (src / "vdwin_utils.mod.F90").write_text(VDW, encoding="utf-8")
        (src / "control_utils.mod.F90").write_text(CPMD, encoding="utf-8")
        (src / "cplngs_utils.mod.F90").write_text(VECTORS, encoding="utf-8")
        self.snapshot = extract_tree(root, sha="0123456789abcdef")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_s6grimme_under_vdw_and_grimme_absent(self) -> None:
        literals = tested_literals(self.snapshot, "VDW")
        self.assertIn("S6GRIMME", literals)
        self.assertNotIn("GRIMME", literals)

    def test_comment_keyword_is_not_a_token(self) -> None:
        self.assertNotIn("ONLYCOMMENT", tested_literals(self.snapshot, "VDW"))

    def test_cpmd_keyword_is_not_attributed_to_vectors(self) -> None:
        self.assertIn("NEWORTHO", tested_literals(self.snapshot, "CPMD"))
        self.assertNotIn("NEWORTHO", tested_literals(self.snapshot, "VECTORS"))
        self.assertIn("FIXED", tested_literals(self.snapshot, "VECTORS"))


if __name__ == "__main__":
    raise SystemExit(unittest.main())
