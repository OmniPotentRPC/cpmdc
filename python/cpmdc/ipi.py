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
import os
import re
import socket
import struct
import sys
import time
from typing import IO, Any, Sequence

from cpmdc import BOHR_TO_ANGSTROM, CpmdcError, Session

__all__ = [
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

_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es "
    "Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
_NUMBER = {symbol: index + 1 for index, symbol in enumerate(_SYMBOLS)}
_FORMULA_PART = re.compile(r"([A-Z][a-z]?)(\d*)")


class ProtocolError(CpmdcError):
    """The server sent a message this driver cannot answer."""


def species_numbers(text: str) -> list[int]:
    """Atomic numbers from element symbols in server order.

    Symbols may be separated by spaces or commas (``"O H H"``), run
    together with counts (``"OH2"``), or given as atomic numbers
    (``"8 1 1"``).
    """
    numbers: list[int] = []
    for token in re.split(r"[\s,]+", text.strip()):
        if not token:
            continue
        if token.isdigit():
            number = int(token)
            if not 1 <= number <= len(_SYMBOLS):
                raise ValueError(f"atomic number {number} is outside 1..{len(_SYMBOLS)}")
            numbers.append(number)
            continue
        parts = _FORMULA_PART.findall(token)
        if "".join(symbol + count for symbol, count in parts) != token:
            raise ValueError(f"cannot read species {token!r}")
        for symbol, count in parts:
            if symbol not in _NUMBER:
                raise ValueError(f"unknown element {symbol!r}")
            numbers.extend([_NUMBER[symbol]] * (int(count) if count else 1))
    if not numbers:
        raise ValueError("species must name at least one atom")
    return numbers


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
    """
    if (unix is None) == (port is None):
        raise ValueError("give exactly one of unix and port")
    where = unix_socket_path(unix, prefix) if unix is not None else f"{host}:{port}"
    deadline = time.monotonic() + max(wait, 0.0)
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
                raise CpmdcError(f"cannot connect to {where}: {exc}") from exc
            time.sleep(0.05)


def _header(name: str) -> bytes:
    return name.encode("ascii").ljust(HEADER_BYTES)


def _det3(m: Sequence[float]) -> float:
    return (
        m[0] * (m[4] * m[8] - m[5] * m[7])
        - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
    )


def _batch_size(init: bytes) -> int:
    text = init.decode("utf-8", "replace")
    for token in text.split(","):
        key, sep, value = token.partition(":")
        if sep and key.strip() == "batch_size":
            try:
                return int(value.strip())
            except ValueError:
                return 1
    return 1


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
        return struct.unpack("=i", raw)[0]

    def raw(self, nbytes: int, what: str) -> bytes:
        raw = self._read(nbytes, what)
        assert raw is not None
        return raw


class Driver:
    """Answer i-PI requests from one :class:`cpmdc.Session`.

    ``params`` is a ``CPMDParams`` path or its bytes, as for
    :class:`cpmdc.Session`. ``numbers`` are the atomic numbers in the order
    the server sends positions. ``cluster=True`` ignores the cell the server
    sends, so the deck sets the box; otherwise an all-zero cell (what ASE
    sends for a structure without a cell) is a cluster and any other cell is
    periodic. ``log`` receives one line per protocol message.
    """

    def __init__(
        self,
        params: bytes | str | os.PathLike[str],
        numbers: Sequence[int],
        *,
        library: str | os.PathLike[str] | None = None,
        cluster: bool = False,
        log: IO[str] | None = None,
    ) -> None:
        self.numbers = [int(number) for number in numbers]
        if not self.numbers:
            raise ValueError("numbers must name at least one atom")
        for number in self.numbers:
            if not 1 <= number <= len(_SYMBOLS):
                raise ValueError(f"atomic number {number} is outside 1..{len(_SYMBOLS)}")
        self.cluster = cluster
        self.steps = 0
        self._log = log
        self._warned_virial = False
        self._session = Session(params, library=library)

    def _say(self, text: str) -> None:
        if self._log is not None:
            print(f"cpmdc-ipi: {text}", file=self._log, flush=True)

    def _cell_rows(self, h: Sequence[float]) -> list[list[float]] | None:
        if self.cluster or not any(h):
            return None
        rows = [[h[3 * i + k] * BOHR_TO_ANGSTROM for i in range(3)] for k in range(3)]
        lengths = [sum(c * c for c in row) ** 0.5 for row in rows]
        if min(lengths) <= 0.0 or abs(_det3(h)) * BOHR_TO_ANGSTROM**3 <= 1.0e-8 * (
            lengths[0] * lengths[1] * lengths[2]
        ):
            raise ProtocolError(
                "the server sent a cell without three independent vectors; send an "
                "all-zero cell or start cpmdc-ipi with --cluster for an isolated system"
            )
        return rows

    def _check_count(self, n_atoms: int) -> None:
        if n_atoms != len(self.numbers):
            raise ProtocolError(
                f"the server sent {n_atoms} atoms, but this driver was started for "
                f"{len(self.numbers)}; the i-PI protocol carries no species, so a "
                "changed atom count or composition needs a new cpmdc-ipi process "
                "with matching --species"
            )

    def compute(
        self, cell: Sequence[float], positions: Sequence[float]
    ) -> tuple[float, list[float], list[float]]:
        """Energy, forces, and virial for one step, in atomic units.

        ``cell`` is the wire's nine reals: h, lattice vectors as columns,
        row-major, in Bohr. ``positions`` are 3N Bohr. The result is the
        energy in Hartree, 3N forces in Hartree/Bohr, and the row-major
        virial ``-dE/d(strain)`` in Hartree, zero for a cluster.
        """
        if len(positions) % 3:
            raise ProtocolError("positions must have three components per atom")
        n_atoms = len(positions) // 3
        self._check_count(n_atoms)
        rows = self._cell_rows(cell)
        coords = [
            [positions[3 * i + k] * BOHR_TO_ANGSTROM for k in range(3)] for i in range(n_atoms)
        ]
        done = self._session.energy_forces(coords, self.numbers, rows)
        forces = [component for row in done.forces_h_bohr for component in row]
        virial = [0.0] * 9
        if rows is not None and done.stress_h_bohr3 is not None:
            volume = abs(_det3(cell))
            virial = [volume * component for component in done.stress_h_bohr3]
        elif rows is not None and not self._warned_virial:
            self._warned_virial = True
            print(
                "cpmdc-ipi: warning: the library computed no stress for this periodic "
                "cell (isolated deck or CPMDC_STRESS=0); the virial sent is zero",
                file=sys.stderr,
                flush=True,
            )
        self.steps += 1
        return float(done.energy_h), forces, virial

    def serve(self, sock: socket.socket) -> int:
        """Answer the server on ``sock`` until ``EXIT`` and return the step count.

        A server that closes the connection between messages ends the run
        like ``EXIT``: ASE's ``SocketIOCalculator.close`` does that. A
        connection closed in the middle of a message, a message out of
        order, an unknown header, or a batched ``INIT`` raises
        :class:`ProtocolError`.
        """
        wire = _Wire(sock)
        initialised = False
        pending: tuple[float, list[float], list[float]] | None = None
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
                batch = _batch_size(text)
                if batch > 1:
                    raise ProtocolError(
                        f"the server asked for batch_size {batch}; cpmdc-ipi evaluates "
                        "one structure per request, so set batch_size to 1"
                    )
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
                pending = self.compute(cell, positions)
                initialised = True
            elif head == "GETFORCE":
                if pending is None:
                    raise ProtocolError("GETFORCE arrived before POSDATA")
                energy, forces, virial = pending
                n_atoms = len(forces) // 3
                transposed = [virial[3 * j + i] for i in range(3) for j in range(3)]
                sock.sendall(
                    _header("FORCEREADY")
                    + struct.pack("=d", energy)
                    + struct.pack("=i", n_atoms)
                    + struct.pack(f"={3 * n_atoms}d", *forces)
                    + struct.pack("=9d", *transposed)
                    + struct.pack("=i", 0)
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
        raise CpmdcError("--structure needs ASE: pip install 'cpmdc[ase]'") from exc
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
    except (CpmdcError, ValueError, OSError) as exc:
        print(f"cpmdc-ipi: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
