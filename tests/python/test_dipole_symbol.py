"""Old and new libcpmdc builds at the same ABI.

An older library has every symbol the loader requires except
``cpmdc_session_last_dipole``. It still evaluates energy and forces.
Asking ASE for the dipole names the missing symbol. A library that
exports the symbol returns the dipole from the last evaluation.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError

from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Potential, Session
from cpmdc.ase import CPMDCCalculator

PARAMS = b"shim-params"
POSITION = [[0.0, 0.0, 0.0]]
DIPOLE = (0.1, 0.2, 0.3)

_SHIM = r"""
#include <stddef.h>
#include <string.h>

typedef struct {
  int ok;
  double energy_h;
  char message[512];
} CPMDCResult;

typedef struct {
  int valid;
  double values[3];
} CPMDCDipole;

static int live = 1;

void *cpmdc_session_create(const void *buf, size_t n) {
  (void)buf;
  return n ? (void *)&live : NULL;
}

void cpmdc_session_destroy(void *session) { (void)session; }

static CPMDCResult ok_result(void) {
  CPMDCResult result;
  memset(&result, 0, sizeof result);
  result.ok = 1;
  result.energy_h = 1.25;
  return result;
}

CPMDCResult cpmdc_session_energy_forces(void *session, int n, const double *pos,
                                       const int *numbers, double *forces) {
  (void)session;
  (void)pos;
  (void)numbers;
  if (forces) {
    for (int i = 0; i < n * 3; ++i)
      forces[i] = 0.1;
  }
  return ok_result();
}

CPMDCResult cpmdc_session_calculate_forces(void *session, const void *message,
                                          size_t message_len, double *forces,
                                          size_t nforces) {
  (void)session;
  (void)message;
  (void)message_len;
  if (forces) {
    for (size_t i = 0; i < nforces; ++i)
      forces[i] = 0.1;
  }
  return ok_result();
}

int cpmdc_session_calculate_au(void *session, int n, const double *pos,
                               const int *numbers, const double *cell,
                               double *energy, double *forces, double *virial) {
  (void)session;
  (void)pos;
  (void)numbers;
  (void)cell;
  if (energy)
    *energy = 1.25;
  if (forces) {
    for (int i = 0; i < n * 3; ++i)
      forces[i] = 0.1;
  }
  if (virial) {
    for (int i = 0; i < 9; ++i)
      virial[i] = 0.0;
  }
  return 0;
}

int cpmdc_session_select_orbitals(void *session, long long key) {
  (void)session;
  (void)key;
  return 0;
}

int cpmdc_adopt_calculator_comm(void *session, size_t n, int rank) {
  (void)session;
  (void)n;
  (void)rank;
  return 0;
}

int cpmdc_session_last_stress(void *session, void *out) {
  (void)session;
  if (out)
    memset(out, 0, 80);
  return -1;
}

const char *cpmdc_last_error(void) { return ""; }

const char *cpmdc_version(void) { return "cpmdc/shim"; }

int cpmdc_abi_version(void) { return 0; }

int cpmdc_available(void) { return 1; }

#ifdef CPMDC_SHIM_DIPOLE
int cpmdc_session_last_dipole(void *session, CPMDCDipole *out) {
  (void)session;
  if (!out)
    return -1;
  out->valid = 1;
  out->values[0] = 0.1;
  out->values[1] = 0.2;
  out->values[2] = 0.3;
  return 0;
}
#endif
"""


def _library(tmp_path: Path, *, dipole: bool) -> Path:
    gcc = shutil.which("gcc")
    if gcc is None:
        pytest.skip("gcc is required to build the dipole shim")
    source = tmp_path / ("new.c" if dipole else "old.c")
    shared = tmp_path / ("libnew.so" if dipole else "libold.so")
    source.write_text(_SHIM, encoding="utf-8")
    command = [gcc, "-shared", "-fPIC", "-o", str(shared), str(source)]
    if dipole:
        command.insert(1, "-DCPMDC_SHIM_DIPOLE")
    subprocess.run(command, check=True)
    return shared


def test_an_older_library_loads_and_refuses_a_requested_dipole(tmp_path: Path) -> None:
    library = _library(tmp_path, dipole=False)
    with Session(PARAMS, library=library) as session:
        assert session.abi_version() == 0
        assert session.available() is True
        assert session.dipole_possible() is False
        done = session.energy_forces(POSITION, [1])
        assert done.energy_h == 1.25
        assert done.dipole_e_bohr is None
        assert session.dipole_possible() is False
        energy, forces, virial = session.calculate_au([0.0, 0.0, 0.0], [1])
        assert energy == 1.25
        assert list(forces) == pytest.approx([0.1, 0.1, 0.1])
        assert list(virial) == pytest.approx([0.0] * 9)
        assert session._missing_dipole_symbol() is not None
    with Potential(PARAMS, library=library) as pot:
        result = pot.evaluate(POSITION, [1])
        assert "dipole" not in result
        assert pot.caps() == {"stress": True, "dipole": False, "energies": False}
    atoms = Atoms("H", positions=POSITION)
    with CPMDCCalculator(PARAMS, library=library) as calc:
        atoms.calc = calc
        assert atoms.get_potential_energy() == pytest.approx(1.25 * HARTREE_TO_EV)
        with pytest.raises(PropertyNotImplementedError, match="cpmdc_session_last_dipole"):
            atoms.get_dipole_moment()


def test_a_library_with_the_dipole_symbol_returns_it(tmp_path: Path) -> None:
    library = _library(tmp_path, dipole=True)
    with Session(PARAMS, library=library) as session:
        assert session.dipole_possible() is True
        done = session.energy_forces(POSITION, [1])
        assert done.dipole_e_bohr == pytest.approx(DIPOLE)
        assert session.dipole_possible() is True
        assert session._missing_dipole_symbol() is None
    with Potential(PARAMS, library=library) as pot:
        assert pot.caps()["dipole"] is True
        result = pot.evaluate(POSITION, [1])
        assert result["dipole"] == pytest.approx([value * BOHR_TO_ANGSTROM for value in DIPOLE])
        assert pot.caps() == {"stress": True, "dipole": True, "energies": False}
    atoms = Atoms("H", positions=POSITION)
    with CPMDCCalculator(PARAMS, library=library) as calc:
        atoms.calc = calc
        np.testing.assert_allclose(
            atoms.get_dipole_moment(),
            np.asarray(DIPOLE) * BOHR_TO_ANGSTROM,
        )
