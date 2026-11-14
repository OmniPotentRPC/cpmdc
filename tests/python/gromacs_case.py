"""Build, run, and read one GROMACS system on cpmdc for test_gromacs.py.

GROMACS is the double-precision ``gmx_d`` named by ``CPMDC_GMX`` or found on
``PATH``, configured with ``-DGMX_CP2K=ON`` against cpmdc's ``libcp2k.pc``.
:func:`write_system` writes a .gro, a .top of flexible TIP3P waters, an
index file, an .mdp, and the ``-qmi`` input; :func:`run` calls grompp and
mdrun; :func:`read_trr` and :func:`energies` read the results at full
precision. :func:`ewald_qm_mm` is an independent NumPy evaluation of the
QM-MM Coulomb sum of mechanical embedding.
"""

from __future__ import annotations

import math
import os
import shlex
import shutil
import struct
import subprocess
from pathlib import Path
from typing import Any, Sequence

import numpy as np

# Two waters, nm; the same molecules as the LAMMPS cases, inside every box.
WATERS_NM = [
    [0.10, 0.10, 0.11173],
    [0.10, 0.17572, 0.05308],
    [0.10, 0.02428, 0.05308],
    [0.50, 0.40, 0.41],
    [0.56, 0.47, 0.42],
    [0.43, 0.44, 0.36],
]
SYMBOLS = "OHHOHH"
TIP3P = {"O": -0.834, "H": 0.417}

TOPOLOGY_HEAD = """[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1 2 yes 0.5 0.8333

[ atomtypes ]
; name at.num mass charge ptype sigma epsilon
OW 8 15.9994 0.0 A 0.315061 0.636386
HW 1 1.008 0.0 A 0.0 0.0
OQ 8 15.9994 0.0 A 0.0 0.0
HQ 1 1.008 0.0 A 0.0 0.0
HL 1 1.008 0.0 A 0.2 0.2
"""

MOLECULE = """
[ moleculetype ]
{name} 2

[ atoms ]
1 {o} 1 {name} OW 1 -0.834 15.9994
2 {h} 1 {name} HW1 1 0.417 1.008
3 {h} 1 {name} HW2 1 0.417 1.008

[ bonds ]
1 2 1 0.09572 502416.0
1 3 1 0.09572 502416.0

[ angles ]
2 1 3 1 104.52 628.02
"""


def gmx_problem() -> tuple[str | None, str]:
    """The gmx_d to run, or None and why it cannot run."""
    named = os.environ.get("CPMDC_GMX")
    gmx = named or shutil.which("gmx_d")
    if not gmx:
        return None, "no GROMACS: set CPMDC_GMX or put gmx_d (GMX_DOUBLE, GMX_CP2K) on PATH"
    try:
        done = subprocess.run(
            [gmx, "-quiet", "--version"], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, f"{gmx} --version failed: {exc}"
    out = done.stdout + done.stderr
    if done.returncode != 0:
        return None, f"{gmx} --version failed: {out.strip()[-300:]}"
    if "Precision:           double" not in out and "Precision: double" not in out:
        return None, f"{gmx} is not a double-precision build (GMX_DOUBLE=ON)"
    return gmx, ""


def gmx_mpi_problem() -> tuple[list[str] | None, str | None, str]:
    """The MPI launcher and gmx_mpi_d to run, or Nones and why they cannot run.

    ``CPMDC_GMX_MPI`` names a gmx_mpi_d built with ``-DGMX_MPI=ON`` against
    the MPI that libcpmdc_gromacs links; ``CPMDC_MPIEXEC`` is its launcher, a
    command line such as ``mpiexec --oversubscribe`` (default ``mpiexec`` on
    PATH).
    """
    gmx = os.environ.get("CPMDC_GMX_MPI") or shutil.which("gmx_mpi_d")
    if not gmx:
        return None, None, "no MPI GROMACS: set CPMDC_GMX_MPI or put gmx_mpi_d on PATH"
    launcher = shlex.split(os.environ.get("CPMDC_MPIEXEC", "")) or [shutil.which("mpiexec")]
    if not launcher[0]:
        return None, None, "no MPI launcher: set CPMDC_MPIEXEC or put mpiexec on PATH"
    try:
        done = subprocess.run(
            [gmx, "-quiet", "--version"], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, None, f"{gmx} --version failed: {exc}"
    out = done.stdout + done.stderr
    if done.returncode != 0:
        return None, None, f"{gmx} --version failed: {out.strip()[-300:]}"
    if "Precision:           double" not in out:
        return None, None, f"{gmx} is not a double-precision build (GMX_DOUBLE=ON)"
    if "MPI library:         MPI" not in out:
        return None, None, f"{gmx} is not built with an external MPI library (GMX_MPI=ON)"
    return launcher, gmx, ""


def water_lattice(n: int, spacing: float, seed: int = 7) -> np.ndarray:
    """n^3 rigidly placed TIP3P-shaped waters on a cubic lattice, nm."""
    rng = np.random.default_rng(seed)
    shape = np.array([[0.0, 0.0, 0.0], [0.09572, 0.0, 0.0],
                      [0.09572 * math.cos(math.radians(104.52)),
                       0.09572 * math.sin(math.radians(104.52)), 0.0]])
    out = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
                centre = (np.array([i, j, k]) + 0.5) * spacing
                out.extend(centre + shape @ q.T)
    return np.array(out)


def gro_box(rows: np.ndarray) -> str:
    """The .gro box line: v1x v2y v3z v1y v1z v2x v2z v3x v3y."""
    r = np.asarray(rows, dtype=float)
    values = [r[0, 0], r[1, 1], r[2, 2], r[0, 1], r[0, 2], r[1, 0], r[1, 2], r[2, 0], r[2, 1]]
    if not any(values[3:]):
        values = values[:3]
    return " ".join(f"{v:.10f}" for v in values)


def write_system(
    workdir: Path,
    *,
    positions: Sequence[Sequence[float]],
    box: Sequence[Sequence[float]],
    qm_molecules: Sequence[int],
    qmi: str,
    mdp: dict[str, Any],
    qm_lj: bool = False,
) -> None:
    """Waters ``i`` in ``qm_molecules`` are QM.

    QM waters take the zero-LJ types OQ and HQ, so that every force on them
    comes from cpmdc, unless ``qm_lj`` gives them TIP3P's oxygen and a
    repulsive hydrogen HL: the reference evaluator does not bind the QM
    hydrogens, and a free TIP3P hydrogen falls into an MM oxygen.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    x = np.asarray(positions, dtype=float)
    nmol = len(x) // 3
    lines = ["waters for cpmdc", str(len(x))]
    for i, row in enumerate(x):
        mol, k = divmod(i, 3)
        name = ("OW", "HW1", "HW2")[k]
        lines.append(f"{mol + 1:5d}{'SOL':<5s}{name:>5s}{i + 1:5d}" + "".join(f"{v:15.10f}" for v in row))
    lines.append(gro_box(np.asarray(box)))
    (workdir / "conf.gro").write_text("\n".join(lines) + "\n")

    top = [TOPOLOGY_HEAD]
    top.append(MOLECULE.format(name="SOL", o="OW", h="HW"))
    top.append(MOLECULE.format(name="QMW", o="OW" if qm_lj else "OQ", h="HL" if qm_lj else "HQ"))
    top.append("\n[ system ]\nwaters\n\n[ molecules ]")
    for mol in range(nmol):
        top.append(f"{'QMW' if mol in qm_molecules else 'SOL'} 1")
    (workdir / "topol.top").write_text("\n".join(top) + "\n")

    qm_atoms = [3 * m + k + 1 for m in sorted(qm_molecules) for k in range(3)]
    everything = list(range(1, len(x) + 1))
    (workdir / "index.ndx").write_text(
        "[ System ]\n"
        + " ".join(map(str, everything))
        + "\n[ QM ]\n"
        + " ".join(map(str, qm_atoms))
        + "\n"
    )
    settings = {
        "integrator": "md",
        "nsteps": 0,
        "dt": 0.00025,
        "cutoff-scheme": "Verlet",
        "coulombtype": "PME",
        "rcoulomb": 0.5,
        "rvdw": 0.5,
        "nstxout": 1,
        "nstvout": 1,
        "nstfout": 1,
        "nstenergy": 1,
        "nstcalcenergy": 1,
        "tcoupl": "no",
        "pcoupl": "no",
        "constraints": "none",
        # The reference evaluator is a well anchored at the origin: the net
        # force is not zero, and removing the centre-of-mass motion would
        # take energy out of an NVE run.
        "comm-mode": "None",
        "qmmm-cp2k-active": "true",
        "qmmm-cp2k-qmgroup": "QM",
        "qmmm-cp2k-qmmethod": "INPUT",
    }
    settings.update(mdp)
    (workdir / "md.mdp").write_text("".join(f"{k} = {v}\n" for k, v in settings.items()))
    (workdir / "cpmdc.inp").write_text(qmi)


def run(
    gmx: str,
    workdir: Path,
    *,
    check: bool = True,
    with_qmi: bool = True,
    launcher: Sequence[str] | None = None,
    ranks: int = 1,
    mdrun_args: Sequence[str] = (),
    env_extra: dict[str, str] | None = None,
    timeout: float | None = None,
) -> subprocess.CompletedProcess:
    """grompp, then mdrun; returns mdrun's process.

    Without ``launcher``, mdrun runs one thread-MPI rank. With it, ``gmx`` is
    an MPI build and mdrun runs as ``launcher -n ranks gmx mdrun``.
    """
    env = dict(os.environ, GMX_MAXBACKUP="-1", **(env_extra or {}))
    qmi = ["-qmi", "cpmdc.inp"] if with_qmi else []
    grompp = subprocess.run(
        [gmx, "-quiet", "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top",
         "-n", "index.ndx", *qmi, "-o", "md.tpr", "-maxwarn", "10"],
        cwd=workdir, capture_output=True, text=True, env=env,
    )
    if grompp.returncode != 0:
        raise RuntimeError(f"grompp failed:\n{grompp.stdout}\n{grompp.stderr}")
    if launcher is None:
        command = [gmx, "-quiet", "mdrun", "-deffnm", "md", "-ntmpi", "1"]
    else:
        command = [*launcher, "-n", str(ranks), gmx, "-quiet", "mdrun", "-deffnm", "md"]
    command += ["-ntomp", "1", "-nb", "cpu", *mdrun_args]
    mdrun = subprocess.run(
        command, cwd=workdir, capture_output=True, text=True, env=env, timeout=timeout
    )
    if check and mdrun.returncode != 0:
        raise RuntimeError(f"mdrun failed:\n{mdrun.stdout}\n{mdrun.stderr}")
    return mdrun


def read_trr(path: Path) -> list[dict[str, Any]]:
    """Frames of a .trr: step, time, box rows, x, v, f (nm, ps, kJ/mol/nm)."""
    data = path.read_bytes()
    pos = 0

    def ints(n: int) -> tuple[int, ...]:
        nonlocal pos
        out = struct.unpack_from(f">{n}i", data, pos)
        pos += 4 * n
        return out

    frames = []
    while pos < len(data):
        magic, _slen, length = ints(3)
        if magic != 1993:
            raise ValueError(f"{path}: not a trr frame at byte {pos}")
        pos += (length + 3) // 4 * 4
        (_ir, _e, box_size, vir_size, pres_size, _top, _sym, x_size, v_size, f_size,
         natoms, step, _nre) = ints(13)
        size = 8 if box_size == 72 or x_size == 24 * natoms or f_size == 24 * natoms else 4
        fmt = ">d" if size == 8 else ">f"

        def reals(n: int) -> np.ndarray:
            nonlocal pos
            out = np.frombuffer(data, dtype=fmt, count=n, offset=pos).astype(float)
            pos += size * n
            return out

        time = float(reals(1)[0])
        reals(1)
        frame: dict[str, Any] = {"step": step, "time": time, "double": size == 8}
        if box_size:
            frame["box"] = reals(9).reshape(3, 3)
        if vir_size:
            reals(9)
        if pres_size:
            reals(9)
        for key, nbytes in (("x", x_size), ("v", v_size), ("f", f_size)):
            if nbytes:
                frame[key] = reals(3 * natoms).reshape(natoms, 3)
        frames.append(frame)
    return frames


def energies(gmx: str, workdir: Path, *terms: str) -> np.ndarray:
    """Columns time and ``terms`` of md.edr, printed by gmx energy -dp."""
    done = subprocess.run(
        [gmx, "-quiet", "energy", "-f", "md.edr", "-o", "energy.xvg", "-dp"],
        cwd=workdir, input="\n".join(terms) + "\n\n", capture_output=True, text=True,
    )
    if done.returncode != 0:
        raise RuntimeError(f"gmx energy failed:\n{done.stdout}\n{done.stderr}")
    rows = [
        [float(v) for v in line.split()]
        for line in (workdir / "energy.xvg").read_text().splitlines()
        if line and line[0] not in "#@"
    ]
    return np.array(rows)


# --- QM-MM Coulomb -------------------------------------------------------


def ewald_alpha(rc: float, tol: float) -> float:
    """GROMACS' calc_ewaldcoeff_q: erfc(alpha rc) = tol by bisection."""
    beta, i = 5.0, 0
    while True:
        i += 1
        beta *= 2.0
        if math.erfc(beta * rc) <= tol:
            break
    low, high = 0.0, beta
    for _ in range(i + 60):
        beta = 0.5 * (low + high)
        if math.erfc(beta * rc) > tol:
            low = beta
        else:
            high = beta
    return beta


def default_rcut(rows: np.ndarray, bohr_to_nm: float) -> float:
    volume = abs(np.linalg.det(rows))
    heights = [volume / np.linalg.norm(np.cross(rows[(a + 1) % 3], rows[(a + 2) % 3])) for a in range(3)]
    return min(1.0 / bohr_to_nm, 0.5 * min(heights))


def ewald_qm_mm(xq, qq, xm, qm, rows, rc, tol):
    """Ewald QM-MM Coulomb energy and forces, atomic units, every k and image.

    The real-space sum takes every image pair closer than ``rc``; the
    reciprocal sum every k with 0 < |k| <= 2 alpha sqrt(-ln tol) over the
    full sphere; the background term is -pi Q_QM Q_MM / (V alpha^2).
    """
    xq, xm = np.asarray(xq, float), np.asarray(xm, float)
    qq, qm = np.asarray(qq, float), np.asarray(qm, float)
    rows = np.asarray(rows, float)
    alpha = ewald_alpha(rc, tol)
    kmax = 2.0 * alpha * math.sqrt(-math.log(tol))
    volume = abs(np.linalg.det(rows))
    inv = np.linalg.inv(rows)
    fq, fm = np.zeros_like(xq), np.zeros_like(xm)
    energy = 0.0
    heights = [volume / np.linalg.norm(np.cross(rows[(a + 1) % 3], rows[(a + 2) % 3])) for a in range(3)]
    reach = [int(math.ceil(rc / h)) + 2 for h in heights]
    shifts = np.array(
        [[i, j, k] for i in range(-reach[0], reach[0] + 1)
         for j in range(-reach[1], reach[1] + 1) for k in range(-reach[2], reach[2] + 1)],
        dtype=float,
    ) @ rows
    for i in range(len(xq)):
        for j in range(len(xm)):
            d0 = xq[i] - xm[j]
            d0 = d0 - np.rint(d0 @ inv) @ rows
            d = d0 + shifts
            r = np.linalg.norm(d, axis=1)
            keep = r < rc
            d, r = d[keep], r[keep]
            e = np.array([math.erfc(alpha * v) for v in r]) / r
            energy += qq[i] * qm[j] * e.sum()
            g = qq[i] * qm[j] * (e + 2.0 / math.sqrt(math.pi) * alpha * np.exp(-(alpha * r) ** 2)) / r**2
            fq[i] += (g[:, None] * d).sum(axis=0)
            fm[j] -= (g[:, None] * d).sum(axis=0)
    recip = 2.0 * math.pi * inv.T
    nk = [int(math.floor(kmax * np.linalg.norm(rows[a]) / (2.0 * math.pi))) for a in range(3)]
    grid = np.array(
        [[i, j, k] for i in range(-nk[0], nk[0] + 1)
         for j in range(-nk[1], nk[1] + 1) for k in range(-nk[2], nk[2] + 1)
         if (i, j, k) != (0, 0, 0)],
        dtype=float,
    )
    kv = grid @ recip
    k2 = (kv**2).sum(axis=1)
    kv, k2 = kv[k2 <= kmax**2], k2[k2 <= kmax**2]
    f = np.exp(-k2 / (4.0 * alpha**2)) / k2
    pq = np.exp(1j * xq @ kv.T)
    pm = np.exp(1j * xm @ kv.T)
    sq, sm = qq @ pq, qm @ pm
    pref = 4.0 * math.pi / volume
    energy += pref * float(np.sum(f * np.real(sq * np.conj(sm))))
    fq += pref * (qq[:, None] * np.imag(pq * np.conj(sm)[None, :]) * f[None, :]) @ kv
    fm += pref * (qm[:, None] * np.imag(pm * np.conj(sq)[None, :]) * f[None, :]) @ kv
    energy -= math.pi * qq.sum() * qm.sum() / (volume * alpha**2)
    return energy, fq, fm


def direct_qm_mm(xq, qq, xm, qm):
    """Plain Coulomb sum between QM and MM charges, atomic units."""
    xq, xm = np.asarray(xq, float), np.asarray(xm, float)
    d = xq[:, None, :] - xm[None, :, :]
    r = np.linalg.norm(d, axis=2)
    qij = np.outer(qq, qm)
    energy = float((qij / r).sum())
    g = qij / r**3
    fq = (g[:, :, None] * d).sum(axis=1)
    fm = -(g[:, :, None] * d).sum(axis=0)
    return energy, fq, fm
