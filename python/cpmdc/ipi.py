"""i-PI socket client for one cpmdc session.

``cpmdc-ipi`` connects to a server that speaks the i-PI wire protocol:
i-PI itself, ASE's ``SocketIOCalculator`` and ``SocketServer``, or any other
server that drives clients such as LAMMPS ``fix ipi``. The protocol is the
one described in the i-PI documentation, "Distributed execution",
"Communication protocol" (https://docs.ipi-code.org/distributed.html), with
the byte layout of i-PI's ``ipi/interfaces/sockets.py`` and
``drivers/py/driver.py`` and of ASE's ``ase/calculators/socketio.py``
(``IPIProtocol`` and ``SocketClient.irun_rank0``):

* a message header is 12 ASCII bytes, padded with spaces;
* integers are 4-byte ``int32`` and reals 8-byte ``float64`` in host order;
* ``STATUS`` is answered with ``NEEDINIT``, ``READY`` or ``HAVEDATA``;
* ``INIT`` carries a bead index, a byte count, and that many bytes;
* ``POSDATA`` carries the cell matrix h (lattice vectors as columns) and its
  inverse, each as nine reals in row-major order, then the atom count and
  the positions;
* ``GETFORCE`` is answered with ``FORCEREADY``, the energy, the atom count,
  the forces, the virial as nine reals laid out as the cell, and a
  length-prefixed extra string;
* ``EXIT`` ends the run.

The extra string is ``{"dipole":[x,y,z]}``, the molecular dipole in
electron Bohr, when the library recomputed a finite dipole for the step,
and empty otherwise. i-PI's ``parse_extra`` (``ipi/interfaces/utils.py``)
decodes that JSON, and its atomic unit of ``electric-dipole`` is electron
Bohr. ASE's ``SocketIOCalculator`` reads the bytes and ignores them.

Every quantity on the wire is in atomic units: Bohr, Hartree, Hartree/Bohr.
The library computes in those units, so the driver passes them through. A
position goes to :class:`cpmdc.Session` as Angstrom through
``BOHR_TO_ANGSTROM``, which the library divides out again. The virial is
``-dE/d(strain)`` in Hartree, the tensor that ASE's ``SocketClient`` sends as
``-volume * stress`` and that i-PI's ASE driver sends as ``vir_calc``. It is
the library stress times the cell volume in Bohr^3, because the library
tensor has the OpenCPMD sign ``-(1/V) dE/d(strain)``.

The protocol carries no species. The driver takes them once, at start, and
a step with another atom count stops it with :class:`ProtocolError`. One
:class:`cpmdc.Session` serves every step, so wavefunction extrapolation
sees the trajectory in order.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import socket
import struct
import sys
import time
from typing import IO, Any, Sequence

from cpmdc import BOHR_TO_ANGSTROM, CPMDCError, Session, _elements
from cpmdc._rank import _FULL_RANK, _NONFINITE, _VOLUME_UNDERFLOW, _cell_state, _cell_volume
from cpmdc._elements import species_numbers

__all__ = [
    "DEFAULT_SOCKET_PREFIX",
    "HEADER_BYTES",
    "Driver",
    "ProtocolError",
    "connect",
    "main",
    "species_numbers",
    "unix_socket_path",
]

HEADER_BYTES = 12

# i-PI's default sockets_prefix and ASE's actualunixsocketname both prepend
# this to the socket name.
DEFAULT_SOCKET_PREFIX = os.path.join(os.sep, "tmp", "ipi_")

class ProtocolError(CPMDCError):
    """The server sent a message this driver cannot answer."""


def unix_socket_path(name: str, prefix: str | None = None) -> str:
    """Path of the UNIX socket for an i-PI address.

    A name that contains a path separator is used as it is. Any other name
    gets ``prefix``, which defaults to the prefix i-PI and ASE use.
    """
    if os.sep in name:
        return name
    return (DEFAULT_SOCKET_PREFIX if prefix is None else prefix) + name


def connect(
    *,
    unix: str | None = None,
    host: str = "localhost",
    port: int | None = None,
    prefix: str | None = None,
    wait: float = 0.0,
) -> socket.socket:
    """Open the socket to the server.

    Give ``unix`` for a UNIX-domain socket, or ``port`` (and ``host``) for
    TCP. A server that is not listening yet is retried for ``wait`` seconds.
    ``wait`` must be zero or positive.
    """
    if (unix is None) == (port is None):
        raise ValueError("give exactly one of unix and port")
    if not math.isfinite(wait):
        raise ValueError(f"wait must be finite, got {wait}")
    if wait < 0.0:
        raise ValueError(f"wait must be zero or positive, got {wait}")
    where = unix_socket_path(unix, prefix) if unix is not None else f"{host}:{port}"
    deadline = time.monotonic() + wait
    while True:
        try:
            if unix is not None:
                sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                try:
                    sock.connect(where)
                except OSError:
                    sock.close()
                    raise
            else:
                sock = socket.create_connection((host, port))
                sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            return sock
        except (FileNotFoundError, ConnectionRefusedError) as exc:
            if time.monotonic() >= deadline:
                raise CPMDCError(f"cannot connect to {where}: {exc}") from exc
            time.sleep(0.05)


def _header(name: str) -> bytes:
    return name.encode("ascii").ljust(HEADER_BYTES)


_INIT_WHITESPACE = " \t\n\r"
_ASCII_INT_RE = re.compile(r"-?[0-9]+", re.ASCII)


def _batch_size(init: bytes) -> int:
    """Return 1, the only batch size this driver evaluates.

    A missing ``batch_size`` is 1. A present value that is not the
    integer 1 raises :class:`ProtocolError` instead of being read as 1.
    """
    try:
        text = init.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ProtocolError("INIT is not UTF-8") from exc
    found: int | None = None
    for token in text.split(","):
        key, sep, value = token.partition(":")
        if not sep or key.strip(_INIT_WHITESPACE) != "batch_size":
            continue
        if found is not None:
            raise ProtocolError("INIT batch_size appears more than once")
        raw = value.strip(_INIT_WHITESPACE)
        if not _ASCII_INT_RE.fullmatch(raw):
            raise ProtocolError(f"INIT batch_size is not an integer: {raw!r}")
        batch = int(raw, 10)
        if batch != 1:
            raise ProtocolError(
                f"the server asked for batch_size {batch}; cpmdc-ipi evaluates "
                "one structure per request, so set batch_size to 1"
            )
        found = batch
    return 1 if found is None else found


class _Wire:
    def __init__(self, sock: socket.socket) -> None:
        self.sock = sock

    def _read(self, nbytes: int, what: str, at_boundary: bool = False) -> bytes | None:
        chunks: list[bytes] = []
        remaining = nbytes
        while remaining > 0:
            chunk = self.sock.recv(remaining)
            if not chunk:
                if at_boundary and remaining == nbytes:
                    return None
                raise ProtocolError(
                    f"the server closed the connection in the middle of {what}"
                )
            chunks.append(chunk)
            remaining -= len(chunk)
        return b"".join(chunks)

    def header(self) -> str | None:
        raw = self._read(HEADER_BYTES, "a message header", at_boundary=True)
        if raw is None:
            return None
        try:
            return raw.decode("ascii").strip()
        except UnicodeDecodeError as exc:
            raise ProtocolError(f"the server sent a header that is not ASCII: {raw!r}") from exc

    def reals(self, count: int, what: str) -> tuple[float, ...]:
        raw = self._read(8 * count, what)
        assert raw is not None
        return struct.unpack(f"={count}d", raw)

    def int32(self, what: str) -> int:
        raw = self._read(4, what)
        assert raw is not None
        return int(struct.unpack("=i", raw)[0])

    def raw(self, nbytes: int, what: str) -> bytes:
        raw = self._read(nbytes, what)
        assert raw is not None
        return raw


class Driver:
    """Answer i-PI requests from one :class:`cpmdc.Session`.

    ``params`` is a ``CPMDParams`` path or its bytes, as for
    :class:`cpmdc.Session`. ``atomic_numbers`` are the atomic numbers in the
    order the server sends positions. ``cluster=True`` ignores the cell the server
    sends, so the deck sets the box; otherwise an all-zero cell (what ASE
    sends for a structure without a cell) is a cluster and any other cell is
    periodic. ``log`` receives one line per protocol message.
    """

    def __init__(
        self,
        params: bytes | str | os.PathLike[str],
        atomic_numbers: Sequence[int],
        *,
        library: str | os.PathLike[str] | None = None,
        cluster: bool = False,
        log: IO[str] | None = None,
    ) -> None:
        self.atomic_numbers = [int(number) for number in atomic_numbers]
        if not self.atomic_numbers:
            raise ValueError("atomic_numbers must name at least one atom")
        for number in self.atomic_numbers:
            if not 1 <= number <= len(_elements._SYMBOLS):
                raise ValueError(
                    f"atomic number {number} is outside 1..{len(_elements._SYMBOLS)}"
                )
        self.cluster = cluster
        self.steps = 0
        self._log = log
        self._warned_virial = False
        self._warned_dipole = False
        self._session = Session(params, library=library)

    def _say(self, text: str) -> None:
        if self._log is not None:
            print(f"cpmdc-ipi: {text}", file=self._log, flush=True)

    def _cell_rows(self, h: Sequence[float]) -> list[list[float]] | None:
        if self.cluster or not any(h):
            return None
        bohr = [[h[3 * i + k] for i in range(3)] for k in range(3)]
        # cpmdc-ipi-c classifies these same Bohr values with the same routine.
        state = _cell_state(bohr)
        if state == _VOLUME_UNDERFLOW:
            raise ProtocolError("the server sent a cell whose volume is too small for binary64")
        if state not in (_FULL_RANK, _NONFINITE):
            raise ProtocolError(
                "the server sent a cell without three independent vectors; send an "
                "all-zero cell or start cpmdc-ipi with --cluster for an isolated system"
            )
        return [[value * BOHR_TO_ANGSTROM for value in row] for row in bohr]

    def _check_count(self, n_atoms: int) -> None:
        if n_atoms != len(self.atomic_numbers):
            raise ProtocolError(
                f"the server sent {n_atoms} atoms, but this driver was started for "
                f"{len(self.atomic_numbers)}; the i-PI protocol carries no species, so a "
                "changed atom count or composition needs a new cpmdc-ipi process "
                "with matching --species"
            )

    def compute(
        self, positions: Sequence[float], cell: Sequence[float]
    ) -> tuple[float, list[float], list[float]]:
        """Energy, forces, and virial for one step, in atomic units.

        ``positions`` are 3N Bohr. ``cell`` is the wire's nine reals: h,
        lattice vectors as columns, row-major, in Bohr; the arguments come
        in the order of :meth:`cpmdc.Session.calculate_au`, not the wire's. The result is the
        energy in Hartree, 3N forces in Hartree/Bohr, and the row-major
        virial ``-dE/d(strain)`` in Hartree, zero for a cluster.
        """
        energy, forces, virial, _ = self._evaluate(positions, cell)
        return energy, forces, virial

    def _extra(self, dipole: Sequence[float] | None) -> bytes:
        """The FORCEREADY extra string for one dipole, or empty bytes.

        i-PI's ``parse_extra`` uses ``json.loads``, which accepts ``NaN``
        and ``Infinity``. A non-finite dipole is omitted, the same as a
        step that computed none, so the energy, forces, and virial still
        go out and the server does not receive a non-finite dipole.
        """
        if dipole is None or len(dipole) != 3:
            return b""
        if not all(map(math.isfinite, dipole)):
            if not self._warned_dipole:
                self._warned_dipole = True
                print(
                    "cpmdc-ipi: warning: the library returned a dipole that is not "
                    "finite; the extra string is empty for such steps",
                    file=sys.stderr,
                    flush=True,
                )
            return b""
        text = json.dumps(
            {"dipole": [float(value) for value in dipole]},
            separators=(",", ":"),
            allow_nan=False,
        )
        return text.encode("ascii")

    def _evaluate(
        self, positions: Sequence[float], cell: Sequence[float]
    ) -> tuple[float, list[float], list[float], bytes]:
        if len(positions) % 3:
            raise ProtocolError("positions must have three components per atom")
        n_atoms = len(positions) // 3
        self._check_count(n_atoms)
        rows = self._cell_rows(cell)
        if not all(map(math.isfinite, positions)):
            raise ValueError("positions must be finite")
        if rows is not None and not all(map(math.isfinite, cell)):
            raise ValueError("cell must be finite")
        volume = 0.0
        if rows is not None:
            # The rows cpmdc_session_calculate_au() takes, in Bohr, as cpmdc-ipi-c sends them.
            volume = _cell_volume([[cell[3 * i + k] for i in range(3)] for k in range(3)])
            if not math.isfinite(volume):
                raise CPMDCError("cell volume is too large for binary64")
        coords = [
            [positions[3 * i + k] * BOHR_TO_ANGSTROM for k in range(3)] for i in range(n_atoms)
        ]
        done = self._session.energy_forces(coords, self.atomic_numbers, rows)
        forces = [component for row in done.forces_h_bohr for component in row]
        virial = [0.0] * 9
        if rows is not None and done.stress_h_bohr3 is not None:
            virial = [volume * component for component in done.stress_h_bohr3]
            if not all(map(math.isfinite, virial)):
                raise CPMDCError("virial is not finite")
        elif rows is not None and not self._warned_virial:
            self._warned_virial = True
            print(
                "cpmdc-ipi: warning: the library computed no stress for this periodic "
                "cell (isolated deck or CPMDC_STRESS=0); the virial sent is zero",
                file=sys.stderr,
                flush=True,
            )
        extra = self._extra(done.dipole_e_bohr)
        self.steps += 1
        return float(done.energy_h), forces, virial, extra

    def serve(self, sock: socket.socket) -> int:
        """Answer the server on ``sock`` until ``EXIT`` and return the step count.

        A server that closes the connection between messages ends the run
        like ``EXIT``: ASE's ``SocketIOCalculator.close`` does that. A
        connection closed in the middle of a message, a message out of
        order, an unknown header, or an ``INIT`` whose ``batch_size`` is
        not the integer 1 raises :class:`ProtocolError`.
        """
        wire = _Wire(sock)
        initialised = False
        pending: tuple[float, list[float], list[float], bytes] | None = None
        while True:
            head = wire.header()
            if head is None:
                self._say("the server closed the connection")
                return self.steps
            self._say(f"received {head}")
            if head == "EXIT":
                return self.steps
            if head == "STATUS":
                if pending is not None:
                    reply = "HAVEDATA"
                elif initialised:
                    reply = "READY"
                else:
                    reply = "NEEDINIT"
                sock.sendall(_header(reply))
                self._say(f"sent {reply}")
            elif head == "INIT":
                bead = wire.int32("INIT")
                nbytes = wire.int32("INIT")
                if nbytes < 0:
                    raise ProtocolError(f"INIT announced {nbytes} bytes")
                text = wire.raw(nbytes, "INIT")
                _batch_size(text)
                self._say(f"bead {bead}, {nbytes} INIT bytes")
                initialised = True
            elif head == "POSDATA":
                if pending is not None:
                    raise ProtocolError("POSDATA arrived before GETFORCE collected the last step")
                cell = wire.reals(9, "POSDATA")
                wire.reals(9, "POSDATA")
                n_atoms = wire.int32("POSDATA")
                if n_atoms <= 0:
                    raise ProtocolError(f"POSDATA announced {n_atoms} atoms")
                self._check_count(n_atoms)
                positions = wire.reals(3 * n_atoms, "POSDATA")
                pending = self._evaluate(positions, cell)
                initialised = True
            elif head == "GETFORCE":
                if pending is None:
                    raise ProtocolError("GETFORCE arrived before POSDATA")
                energy, forces, virial, extra = pending
                n_atoms = len(forces) // 3
                transposed = [virial[3 * j + i] for i in range(3) for j in range(3)]
                sock.sendall(
                    _header("FORCEREADY")
                    + struct.pack("=d", energy)
                    + struct.pack("=i", n_atoms)
                    + struct.pack(f"={3 * n_atoms}d", *forces)
                    + struct.pack("=9d", *transposed)
                    + struct.pack("=i", len(extra))
                    + extra
                )
                self._say("sent FORCEREADY")
                pending = None
            else:
                raise ProtocolError(f"unknown message {head!r} from the server")

    def close(self) -> None:
        """Close the session."""
        self._session.close()

    def __enter__(self) -> Driver:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def _structure_numbers(path: str) -> list[int]:
    try:
        from ase.io import read
    except ImportError as exc:
        raise CPMDCError("--structure needs ASE: pip install 'cpmdc[ase]'") from exc
    atoms: Any = read(path)
    return [int(number) for number in atoms.get_atomic_numbers()]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cpmdc-ipi",
        description=(
            "Serve cpmdc energies, forces, and virials to an i-PI protocol server "
            "(i-PI, ASE SocketIOCalculator or SocketServer) over a UNIX or TCP socket."
        ),
    )
    where = parser.add_mutually_exclusive_group(required=True)
    where.add_argument("--unix", metavar="NAME", help="UNIX socket address, as in i-PI and ASE")
    where.add_argument("--port", type=int, help="TCP port of the server")
    parser.add_argument("--host", default="localhost", help="TCP host (default: localhost)")
    parser.add_argument("--params", required=True, help="serialized CPMDParams file")
    species = parser.add_mutually_exclusive_group(required=True)
    species.add_argument(
        "--species", help="element symbols in server order: 'O H H', 'O,H,H', or 'OH2'"
    )
    species.add_argument("--structure", help="structure file ASE can read; only species are used")
    parser.add_argument("--library", help="path of libcpmdc (default: the cpmdc search order)")
    parser.add_argument(
        "--cluster", action="store_true", help="ignore the server's cell; the deck sets the box"
    )
    parser.add_argument(
        "--socket-prefix",
        default=None,
        help="prefix for a UNIX socket name (default: i-PI's and ASE's)",
    )
    parser.add_argument(
        "--wait", type=float, default=0.0, help="seconds to retry until the server listens"
    )
    parser.add_argument("--verbose", action="store_true", help="log every protocol message")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run ``cpmdc-ipi`` and return its exit status."""
    args = _parser().parse_args(argv)
    try:
        if args.port is not None and not 0 < args.port < 65536:
            raise ValueError(f"port {args.port} is outside 1..65535")
        numbers = (
            species_numbers(args.species)
            if args.species is not None
            else _structure_numbers(args.structure)
        )
        log = sys.stderr if args.verbose else None
        with Driver(
            args.params, numbers, library=args.library, cluster=args.cluster, log=log
        ) as driver:
            sock = connect(
                unix=args.unix,
                host=args.host,
                port=args.port,
                prefix=args.socket_prefix,
                wait=args.wait,
            )
            with sock:
                steps = driver.serve(sock)
        if args.verbose:
            print(f"cpmdc-ipi: {steps} steps", file=sys.stderr)
        return 0
    except (CPMDCError, ValueError, OSError) as exc:
        print(f"cpmdc-ipi: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
