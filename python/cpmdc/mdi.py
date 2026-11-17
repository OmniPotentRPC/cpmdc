"""MDI engine for one cpmdc session.

``cpmdc-mdi`` (``python -m cpmdc.mdi``) is an engine for the MolSSI Driver
Interface, so a driver such as LAMMPS ``fix mdi/qm`` can request energies,
forces, and stress. The commands are the ``@DEFAULT`` set drivers send,
registered with the MDI library so ``MDI_Check_command_exists`` sees them.
``<NAME`` is answered by the library from ``-name``.

Every quantity on the wire is in atomic units, the MDI standard: Bohr,
Hartree, Hartree/Bohr, and the intensive virial in Hartree/Bohr^3. The
engine passes coordinates, atomic numbers, and the cell straight to
:meth:`cpmdc.Session.calculate_au`, which calls
``cpmdc_session_calculate_au()`` and applies the library's own Bohr
conversion. This module adds no conversion constant.

The cell is nine Bohr values, one lattice vector per row, the layout of
``>CELL`` and of ``calculate_au``. ``>CELL_DISPL`` is stored and does not
move the coordinates: MDI coordinates are Cartesian. An all-zero cell is a
cluster. ``--cluster`` ignores the cell and lets the deck set the box.

``>CHARGES`` is accepted when every value is zero, and when each value
equals the atomic number already sent with ``>ELEMENTS``. Any other charge
is an error: the total charge lives in the ``CPMDParams`` message, not in
per-atom MDI charges.

One :class:`cpmdc.Session` serves every step. The first successful
evaluation fixes the atom count and the elemental composition, the same
check the session makes. A later ``>NATOMS`` or ``>ELEMENTS`` that fails
that check stops the engine with a protocol error and does not replace the
system. A reordering of the same elements stays, as it does on the session.

``pymdi`` is imported when the engine starts. It is an optional dependency,
``pip install 'cpmdc[mdi]'``, and importing :mod:`cpmdc` does not import it.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import IO, Any, Sequence

from cpmdc import CPMDCError, Session

__all__ = [
    "COMMANDS",
    "NODE",
    "Engine",
    "ProtocolError",
    "MDI_Plugin_init_cpmdc",
    "main",
]

NODE = "@DEFAULT"

# Registered at @DEFAULT. <NAME is also answered inside pymdi from -name,
# before this engine sees the command. <KE_ELEC is absent on purpose:
# LAMMPS fix mdi/qm asks for <PE + <KE_ELEC only when both exist, and
# otherwise asks for <ENERGY.
COMMANDS = (
    "<NAME",
    "<@",
    ">NATOMS",
    "<NATOMS",
    ">ELEMENTS",
    ">COORDS",
    "<COORDS",
    ">CELL",
    "<CELL",
    ">CELL_DISPL",
    ">CHARGES",
    "<ENERGY",
    "<PE",
    "<FORCES",
    "<STRESS",
    "EXIT",
)

_RESULTS = frozenset({"<ENERGY", "<PE", "<FORCES", "<STRESS"})


class ProtocolError(CPMDCError):
    """The driver sent a command or a system this engine cannot accept."""


def _pymdi() -> Any:
    try:
        import mdi
    except ImportError as exc:
        raise CPMDCError(
            "pymdi is not installed; install it with pip install 'cpmdc[mdi]'"
        ) from exc
    return mdi


def _det3(m: Sequence[float]) -> float:
    return (
        m[0] * (m[4] * m[8] - m[5] * m[7])
        - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
    )


def _counts(numbers: Sequence[int]) -> dict[int, int]:
    found: dict[int, int] = {}
    for number in numbers:
        found[number] = found.get(number, 0) + 1
    return found


def _flag(options: str, flag: str) -> str | None:
    parts = options.split()
    for index, part in enumerate(parts):
        if part == flag and index + 1 < len(parts):
            return parts[index + 1]
    return None


def _method(options: str) -> str:
    return _flag(options, "-method") or ""


def _mdi_name(text: str) -> str:
    """A length-``MDI_NAME_LENGTH`` string, null padded, for ``MDI_Send``."""
    length = _pymdi().MDI_NAME_LENGTH
    if len(text) >= length:
        raise ProtocolError(f"MDI name {text!r} is longer than {length} bytes")
    return text + "\0" * (length - len(text))


class Engine:
    """Answer MDI ``@DEFAULT`` commands from one :class:`cpmdc.Session`.

    ``params`` is a ``CPMDParams`` path or its bytes, as for
    :class:`cpmdc.Session`. ``cluster=True`` ignores the cell the driver
    sends. ``log`` receives one line per command.
    """

    def __init__(
        self,
        params: bytes | str | os.PathLike[str],
        *,
        library: str | os.PathLike[str] | None = None,
        cluster: bool = False,
        log: IO[str] | None = None,
        name: str = "CPMDC",
    ) -> None:
        self.cluster = cluster
        self.steps = 0
        self._name = name
        self._log = log
        self._failed = ""
        self._reported = False
        self._session = Session(params, library=library)
        self._natoms: int | None = None
        self._elements: list[int] | None = None
        self._coords: list[float] | None = None
        self._cell: list[float] | None = None
        self._cell_displ = [0.0, 0.0, 0.0]
        self._topology: list[int] | None = None
        self._energy = 0.0
        self._forces: list[float] = []
        self._stress = [0.0] * 9
        self._fresh = False
        self._warned_stress = False

    def _say(self, text: str) -> None:
        if self._log is not None:
            print(f"cpmdc-mdi: {text}", file=self._log, flush=True)

    def _fail(self, message: str) -> None:
        self._failed = message
        if not self._reported:
            self._reported = True
            print(f"cpmdc-mdi: error: {message}", file=sys.stderr, flush=True)

    def _register(self) -> None:
        mdi = _pymdi()
        mdi.MDI_Register_Node(NODE)
        for command in COMMANDS:
            mdi.MDI_Register_Command(NODE, command)

    def _recv_doubles(self, count: int, comm: int, what: str) -> list[float]:
        mdi = _pymdi()
        raw = mdi.MDI_Recv(count, mdi.MDI_DOUBLE, comm)
        values = [raw] if count == 1 else [float(value) for value in raw]
        if len(values) != count:
            raise ProtocolError(f"{what} returned {len(values)} values, expected {count}")
        for value in values:
            if value != value or value in (float("inf"), float("-inf")):
                raise ProtocolError(f"{what} must be finite")
        return values

    def _check_topology(self, natoms: int, elements: Sequence[int] | None) -> None:
        if self._topology is None:
            return
        if natoms != len(self._topology):
            raise ProtocolError(
                "topology change requires a new session: "
                f"the driver sent {natoms} atoms after this engine evaluated "
                f"{len(self._topology)}; start another cpmdc-mdi process"
            )
        if elements is not None and _counts(elements) != _counts(self._topology):
            raise ProtocolError(
                "topology change requires a new session: the elemental composition "
                "changed; a reordering of the same elements is accepted, a different "
                "composition needs a new cpmdc-mdi process"
            )

    def _set_natoms(self, natoms: int) -> None:
        if natoms < 1:
            raise ProtocolError(f">NATOMS sent {natoms}; the engine needs at least one atom")
        self._check_topology(natoms, None)
        if self._natoms != natoms:
            self._elements = None
            self._coords = None
            self._fresh = False
        self._natoms = natoms

    def _set_elements(self, elements: Sequence[int]) -> None:
        if self._natoms is None:
            raise ProtocolError("send >NATOMS before >ELEMENTS")
        if len(elements) != self._natoms:
            raise ProtocolError(
                f">ELEMENTS sent {len(elements)} values for {self._natoms} atoms"
            )
        for number in elements:
            if not 1 <= int(number) <= 118:
                raise ProtocolError(f"atomic number {number} is outside 1..118")
        values = [int(number) for number in elements]
        self._check_topology(self._natoms, values)
        if values != self._elements:
            self._fresh = False
        self._elements = values

    def _cell_rows(self, cell: Sequence[float]) -> list[float] | None:
        if not any(cell):
            return None
        lengths = [
            sum(cell[3 * i + k] * cell[3 * i + k] for k in range(3)) ** 0.5 for i in range(3)
        ]
        if min(lengths) <= 0.0 or abs(_det3(cell)) <= 1.0e-8 * (
            lengths[0] * lengths[1] * lengths[2]
        ):
            raise ProtocolError(
                "the driver sent a cell without three independent vectors; send an "
                "all-zero cell, or start cpmdc-mdi with --cluster, for a system with no cell"
            )
        return [float(value) for value in cell]

    def compute(self) -> tuple[float, list[float], list[float]]:
        """Energy, forces, and stress for the system the driver last sent.

        Returns Hartree, ``3N`` Hartree/Bohr forces ``-dE/dR``, and the
        nine-component row-major intensive virial ``-(1/V) dE/d(strain)`` in
        Hartree/Bohr^3. The stress is zero when the step has no cell or the
        library computed no tensor. A repeated call with the same system
        returns the stored step.
        """
        if self._fresh:
            return self._energy, self._forces, self._stress
        if self._natoms is None or self._elements is None or self._coords is None:
            missing = [
                name
                for name, present in (
                    (">NATOMS", self._natoms is not None),
                    (">ELEMENTS", self._elements is not None),
                    (">COORDS", self._coords is not None),
                )
                if not present
            ]
            raise ProtocolError(
                "cannot evaluate yet; the driver has not sent " + ", ".join(missing)
            )
        cell = None if self.cluster else self._cell
        try:
            energy, forces, virial = self._session.calculate_au(
                self._coords, self._elements, cell
            )
        except CPMDCError as exc:
            text = str(exc)
            if "topology" in text:
                raise ProtocolError(text) from exc
            raise
        force_list = [float(value) for value in forces]
        if cell is None:
            stress = [0.0] * 9
        else:
            volume = abs(_det3(cell))
            stress = [float(value) / volume for value in virial] if volume else [0.0] * 9
            if volume and not any(virial) and not self._warned_stress:
                self._warned_stress = True
                print(
                    "cpmdc-mdi: warning: the library computed no stress for this cell "
                    "(isolated deck or CPMDC_STRESS=0); <STRESS is zero",
                    file=sys.stderr,
                    flush=True,
                )
        self._energy = float(energy)
        self._forces = force_list
        self._stress = stress
        self._fresh = True
        self._topology = list(self._elements)
        self.steps += 1
        return self._energy, self._forces, self._stress

    def _handle(self, command: str, comm: int) -> None:
        mdi = _pymdi()
        self._say(f"received {command}")
        if self._failed:
            raise ProtocolError(self._failed)
        if command == "<@":
            mdi.MDI_Send(_mdi_name(NODE), mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, comm)
        elif command == "<NAME":
            # pymdi answers <NAME from -name before this loop sees it. A library
            # that passes the command through gets the same name.
            mdi.MDI_Send(_mdi_name(self._name), mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, comm)
        elif command == ">NATOMS":
            natoms = int(mdi.MDI_Recv(1, mdi.MDI_INT, comm))
            self._set_natoms(natoms)
        elif command == "<NATOMS":
            if self._natoms is None:
                raise ProtocolError("the driver asked for <NATOMS before >NATOMS")
            mdi.MDI_Send(self._natoms, 1, mdi.MDI_INT, comm)
        elif command == ">ELEMENTS":
            if self._natoms is None:
                raise ProtocolError("send >NATOMS before >ELEMENTS")
            raw = mdi.MDI_Recv(self._natoms, mdi.MDI_INT, comm)
            values = [raw] if self._natoms == 1 else raw
            self._set_elements(values)
        elif command == ">COORDS":
            if self._natoms is None:
                raise ProtocolError("send >NATOMS before >COORDS")
            coords = self._recv_doubles(3 * self._natoms, comm, ">COORDS")
            if coords != self._coords:
                self._fresh = False
            self._coords = coords
        elif command == "<COORDS":
            if self._coords is None:
                raise ProtocolError("the driver asked for <COORDS before >COORDS")
            mdi.MDI_Send(self._coords, len(self._coords), mdi.MDI_DOUBLE, comm)
        elif command == ">CELL":
            cell = self._cell_rows(self._recv_doubles(9, comm, ">CELL"))
            if cell != self._cell:
                self._fresh = False
            self._cell = cell
        elif command == "<CELL":
            mdi.MDI_Send(self._cell or [0.0] * 9, 9, mdi.MDI_DOUBLE, comm)
        elif command == ">CELL_DISPL":
            # Cartesian coordinates already include the origin. The displacement
            # is recorded so a driver can send it; the step uses >COORDS as sent.
            self._cell_displ = self._recv_doubles(3, comm, ">CELL_DISPL")
        elif command == ">CHARGES":
            if self._natoms is None:
                raise ProtocolError("send >NATOMS before >CHARGES")
            charges = self._recv_doubles(self._natoms, comm, ">CHARGES")
            zeros = all(charge == 0.0 for charge in charges)
            nuclear = self._elements is not None and all(
                charge == float(number) for charge, number in zip(charges, self._elements)
            )
            if not zeros and not nuclear:
                raise ProtocolError(
                    ">CHARGES is ignored when every charge is 0, and when each charge "
                    "equals the atomic number already sent with >ELEMENTS; any other "
                    "charge is refused, because cpmdc takes the total charge from the "
                    "CPMDParams message, not from per-atom MDI charges"
                )
        elif command in _RESULTS:
            energy, forces, stress = self.compute()
            if command in ("<ENERGY", "<PE"):
                mdi.MDI_Send(energy, 1, mdi.MDI_DOUBLE, comm)
            elif command == "<FORCES":
                mdi.MDI_Send(forces, len(forces), mdi.MDI_DOUBLE, comm)
            else:
                mdi.MDI_Send(stress, 9, mdi.MDI_DOUBLE, comm)
        else:
            raise ProtocolError(
                f"unsupported MDI command {command!r} at {NODE}; "
                "this engine accepts " + ", ".join(COMMANDS)
            )
        self._say(f"answered {command}")

    def serve(self, comm: int) -> int:
        """Answer the driver on ``comm`` until ``EXIT`` and return the step count.

        ``<`` and ``>`` commands that a plugin driver delivers through the
        execute callback are handled there. This loop answers the commands
        ``MDI_Recv_Command`` returns, which for TCP is every command and for
        a plugin is ``EXIT`` after the driver's callback returns.
        """
        mdi = _pymdi()
        while True:
            command = mdi.MDI_Recv_Command(comm)
            if command == "EXIT":
                self._say("received EXIT")
                return self.steps
            code = _execute(command, comm, self)
            if code != 0:
                raise ProtocolError(self._failed or f"MDI command {command} failed")

    def close(self) -> None:
        """Close the session."""
        self._session.close()

    def __enter__(self) -> Engine:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def _execute(command: str, comm: int, engine: Engine) -> int:
    """pymdi execute callback. Zero on success, one on a protocol error."""
    try:
        engine._handle(command, comm)
    except (ProtocolError, CPMDCError, ValueError) as exc:
        engine._fail(str(exc))
        return 1
    return 0


def _params_path(explicit: str | None) -> str:
    if explicit:
        return explicit
    env = os.environ.get("CPMDC_PARAMS")
    if env:
        return env
    raise CPMDCError(
        "pass --params, or set CPMDC_PARAMS, to a serialized CPMDParams file"
    )


def _listen(options: str, engine: Engine) -> int:
    """Initialize pymdi, register ``@DEFAULT``, and serve one connection."""
    mdi = _pymdi()
    method = _method(options)
    if method == "MPI":
        raise ProtocolError(
            "cpmdc-mdi does not use the MDI MPI method; use -method TCP, "
            "or the plugin entry MDI_Plugin_init_cpmdc with -method LINK"
        )
    if method not in ("TCP", "LINK", "PLUG"):
        raise ProtocolError(
            f"cpmdc-mdi does not use -method {method or '(missing)'}; "
            "use TCP, or LINK for the plugin"
        )
    mdi.MDI_Init(options)
    engine._register()
    if method in ("LINK", "PLUG"):
        mdi.MDI_Set_Execute_Command_Func(_execute, engine)
    comm = mdi.MDI_Accept_Communicator()
    return engine.serve(comm)


def _plugin_flag(flag: str) -> str | None:
    """One option from the plugin arguments MDI kept after ``MDI_Init``."""
    mdi = _pymdi()
    argc = mdi.MDI_Plugin_get_argc()
    for index in range(argc):
        found = _flag(mdi.MDI_Plugin_get_arg(index), flag)
        if found:
            return found
    return None


def MDI_Plugin_init_cpmdc(plugin_state: object) -> int:
    """Plugin entry a driver calls after loading this engine.

    ``MDI_Set_plugin_state`` runs ``MDI_Init`` on the plugin's ``-mdi``
    argument, which must be ``-role ENGINE -name CPMDC -method LINK``.
    The ``CPMDParams`` file is ``CPMDC_PARAMS``: MDI rejects engine arguments
    that are not its own, so ``--params`` cannot travel in the plugin
    options. The same function is ``MDI_Plugin_init_mdi`` when the plugin
    name is ``mdi`` and ``-plugin_path`` is the directory of this file.
    """
    mdi = _pymdi()
    mdi.MDI_Set_plugin_state(plugin_state)
    try:
        name = _plugin_flag("-name") or "CPMDC"
    except Exception:
        # pymdi answers <NAME from -name. This name is only the fallback.
        name = "CPMDC"
    engine = Engine(_params_path(None), name=name)
    try:
        engine._register()
        mdi.MDI_Set_Execute_Command_Func(_execute, engine)
        comm = mdi.MDI_Accept_Communicator()
        engine.serve(comm)
    finally:
        engine.close()
    return 0


# The file loader looks up MDI_Plugin_init_<plugin name>. The package file is
# mdi.py, so a driver that sets -plugin_path to this directory uses the name mdi.
MDI_Plugin_init_mdi = MDI_Plugin_init_cpmdc


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cpmdc-mdi",
        description=(
            "Serve cpmdc energies, forces, and stress to an MDI driver "
            "(LAMMPS fix mdi/qm, or any driver that sends @DEFAULT commands) "
            "over TCP. Pass the MDI options as one string, "
            'for example -mdi "-role ENGINE -name CPMDC -method TCP '
            '-hostname localhost -port 8021".'
        ),
    )
    parser.add_argument(
        "--params",
        help="serialized CPMDParams file; default: the CPMDC_PARAMS environment variable",
    )
    parser.add_argument("--library", help="path of libcpmdc (default: the cpmdc search order)")
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="ignore the driver's cell; the CPMDParams deck sets the box",
    )
    parser.add_argument("--verbose", action="store_true", help="log every MDI command")
    return parser


def _take_mdi(argv: Sequence[str]) -> tuple[str | None, list[str]]:
    options = None
    rest: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token in ("-mdi", "--mdi"):
            if index + 1 >= len(argv):
                raise CPMDCError(f"{token} needs the MDI options as the next argument")
            options = argv[index + 1]
            index += 2
            continue
        rest.append(token)
        index += 1
    return options, rest


def main(argv: Sequence[str] | None = None) -> int:
    """Run ``cpmdc-mdi`` and return its exit status."""
    engine: Engine | None = None
    try:
        options, rest = _take_mdi(tuple(sys.argv[1:] if argv is None else argv))
        args = _parser().parse_args(rest)
        if options is None:
            raise CPMDCError(
                'pass -mdi "-role ENGINE -name CPMDC -method TCP -hostname HOST -port PORT"'
            )
        if _method(options) in ("LINK", "PLUG"):
            raise ProtocolError(
                "the LINK method is the plugin entry MDI_Plugin_init_cpmdc; "
                "python -m cpmdc.mdi speaks TCP"
            )
        log = sys.stderr if args.verbose else None
        engine = Engine(
            _params_path(args.params),
            library=args.library,
            cluster=args.cluster,
            log=log,
            name=_flag(options, "-name") or "CPMDC",
        )
        try:
            steps = _listen(options, engine)
        finally:
            engine.close()
        if args.verbose:
            print(f"cpmdc-mdi: {steps} steps", file=sys.stderr)
        return 0
    except (CPMDCError, ValueError, OSError) as exc:
        # The command callback already wrote this message.
        if engine is None or not engine._reported:
            print(f"cpmdc-mdi: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
