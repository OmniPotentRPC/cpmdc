#!/usr/bin/env python3
"""Rank-count invariance of the gamma atomic H and S overlaps.

The stock replicated overlap is a DGEMM whose K dimension is the
rank-local G count, followed by summat of that rounded upper triangle.
The partial changes with the plane split, so ranks 1, 2, 3, and 4
disagree. Under -D__EXACT_ACC each G term enters a superaccumulator,
the digits are merged, and the result is rounded once. ksmat posts one
skip, ovlap posts one, and the two summat calls consume them.

This script runs that assertion on the reducer selected by the exact
patch series. The same assertion on the stock partial-sum reducer
fails. Dropping the ovlap fill, either skip, the summat consume, or
the tdmal guard selects the stock reducer, and the identity check fails.
"""

from __future__ import annotations

import math
import struct
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXACT = ROOT / "bench" / "opencpmd_invariance" / "patches" / "exact"
RANKS = (1, 2, 3, 4)

# Per-G contributions of one matrix element. Contiguous slices are the
# rank split. The exact sum is 4 and 3.5. The stock partial sum is a
# different binary64 value at each of ranks 1, 2, 3, and 4.
H_TERMS = (1.0, 1.0, 1.0, 1e16, -1e16, 1.0)
S_TERMS = (1.0, 1.0, 1.0, 1e16, -1e16, 0.5)


@dataclass(frozen=True)
class Wiring:
    h_exact: bool
    s_exact: bool
    h_note: bool
    s_note: bool
    consume: bool


OLD = Wiring(False, False, False, False, False)


def hex_f64(value: float) -> str:
    return struct.pack(">d", value).hex()


def to_units(value: float) -> int:
    """Signed integer multiple of 2^-1074."""
    bits = struct.unpack("<Q", struct.pack("<d", float(value)))[0]
    sign = bits >> 63
    exponent = (bits >> 52) & 0x7FF
    fraction = bits & ((1 << 52) - 1)
    if exponent == 0x7FF:
        raise ValueError(f"non-finite overlap term {value!r}")
    if exponent == 0:
        if fraction == 0:
            return 0
        mantissa = fraction
        power = -1074
    else:
        mantissa = fraction | (1 << 52)
        power = exponent - 1023 - 52
    if sign:
        mantissa = -mantissa
    return mantissa << (power + 1074)


def round_units(units: int) -> float:
    """One round-ties-to-even of an integer multiple of 2^-1074."""
    if units == 0:
        return 0.0
    sign = 0
    if units < 0:
        sign = 1
        units = -units
    top = units.bit_length() - 1
    exponent = top - 1074

    def pack(exponent_bits: int, fraction: int) -> float:
        bits = (sign << 63) | (exponent_bits << 52) | fraction
        return struct.unpack("<d", struct.pack("<Q", bits))[0]

    if exponent > 1023:
        return pack(0x7FF, 0)
    if exponent >= -1022:
        shift = top - 52
        if shift > 0:
            mantissa = units >> shift
            rest = units & ((1 << shift) - 1)
            round_bit = (rest >> (shift - 1)) & 1
            sticky = rest & ((1 << (shift - 1)) - 1)
        elif shift == 0:
            mantissa, round_bit, sticky = units, 0, 0
        else:
            mantissa, round_bit, sticky = units << (-shift), 0, 0
        if round_bit and (sticky or (mantissa & 1)):
            mantissa += 1
        if mantissa >> 53:
            mantissa >>= 1
            exponent += 1
        if exponent > 1023:
            return pack(0x7FF, 0)
        return pack(exponent + 1023, mantissa & ((1 << 52) - 1))
    if units >= (1 << 52):
        raise RuntimeError(f"subnormal boundary units={units} exp={exponent}")
    return pack(0, units)


def chunks(terms: tuple[float, ...], nrank: int) -> list[tuple[float, ...]]:
    count = len(terms)
    base, remainder = divmod(count, nrank)
    out: list[tuple[float, ...]] = []
    start = 0
    for rank in range(nrank):
        width = base + (1 if rank < remainder else 0)
        out.append(tuple(terms[start : start + width]))
        start += width
    return out


def fill(terms: tuple[float, ...], nrank: int, exact: bool) -> tuple[float | None, list[float]]:
    parts = chunks(terms, nrank)
    if exact:
        total = 0
        for part in parts:
            total += sum(to_units(term) for term in part)
        return round_units(total), []
    partials: list[float] = []
    for part in parts:
        acc = 0.0
        for term in part:
            acc += term
        partials.append(acc)
    return None, partials


def summat(
    nrank: int,
    global_value: float | None,
    partials: list[float],
    skip: list[int],
    consume: bool,
) -> float:
    # summat consumes a pending skip before the one-rank return.
    if consume and skip[0] > 0:
        skip[0] -= 1
        if global_value is not None:
            return global_value
        return partials[0]
    if nrank <= 1:
        if global_value is not None:
            return global_value
        return partials[0]
    if global_value is not None:
        acc = 0.0
        for _ in range(nrank):
            acc += global_value
        return acc
    acc = 0.0
    for partial in partials:
        acc += partial
    return acc


def atomic_hs(nrank: int, wiring: Wiring) -> tuple[float, float]:
    """H then S, matching ksmat, ovlap, summat(xxmat), summat(xsmat)."""
    skip = [0]
    h_global, h_partials = fill(H_TERMS, nrank, wiring.h_exact)
    if wiring.h_note:
        skip[0] += 1
    s_global, s_partials = fill(S_TERMS, nrank, wiring.s_exact)
    if wiring.s_note:
        skip[0] += 1
    h_value = summat(nrank, h_global, h_partials, skip, wiring.consume)
    s_value = summat(nrank, s_global, s_partials, skip, wiring.consume)
    return h_value, s_value


def _section(patch: str, filename: str) -> str:
    text = (EXACT / patch).read_text(encoding="utf-8")
    marker = f"a/src/{filename}"
    start = text.find(marker)
    if start < 0:
        raise AssertionError(f"{patch}: missing {filename}")
    nxt = text.find("\ndiff -ruN ", start)
    return text[start:] if nxt < 0 else text[start:nxt]


def wiring_from_patches() -> Wiring:
    ovlap = _section("0003-exact-overlap.patch", "ovlap_utils.mod.F90")
    ksmat = _section("0003-exact-overlap.patch", "ksmat_utils.mod.F90")
    summat_text = _section("0002-exact-reduction-skip.patch", "summat_utils.mod.F90")
    atomwf = _section("0004-exact-apply-and-guess.patch", "atomwf_utils.mod.F90")
    tdmal = (
        "+#ifndef __EXACT_ACC\n"
        "     IF (.NOT.tkpts%tkpnt.AND..NOT.cntl%cdft) cntl%tdmal = .TRUE.\n"
        "+#endif"
    )
    if tdmal not in atomwf:
        return OLD
    h_exact = "CALL exact_ovlap_fill(n1,n2,a,c1,c2,ngw,ibeg_c0,NGW_local,&" in ovlap
    s_exact = "CALL exact_ovlap_fill(nstate,nstate,a,C1_local,C1_local,&" in ovlap
    h_note = "IF (.NOT.tkpts%tkpnt) CALL exact_note_skip(1)" in ksmat
    s_note = "CALL exact_note_skip(1)" in ovlap
    consume = "IF (exact_consume_skip()) RETURN" in summat_text
    if not (h_exact and s_exact and h_note and s_note and consume):
        return OLD
    return Wiring(True, True, True, True, True)


def assert_rank_identical(wiring: Wiring) -> list[tuple[str, str]]:
    rows = [atomic_hs(nrank, wiring) for nrank in RANKS]
    hexs = [(hex_f64(h_value), hex_f64(s_value)) for h_value, s_value in rows]
    if len(set(hexs)) != 1:
        detail = ", ".join(
            f"r{nrank} H={pair[0]} S={pair[1]}" for nrank, pair in zip(RANKS, hexs)
        )
        raise AssertionError(f"atomic overlaps depend on the rank count: {detail}")
    return hexs


def test_old_partial_sum_fails_rank_identity() -> None:
    try:
        assert_rank_identical(OLD)
    except AssertionError:
        return
    raise AssertionError("stock partial-sum overlap was bit-identical across ranks")


def test_exact_sum_matches_one_rounding() -> None:
    for terms in (H_TERMS, S_TERMS):
        oracle = math.fsum(terms)
        for nrank in RANKS:
            merged, _ = fill(terms, nrank, True)
            if merged != oracle or hex_f64(merged) != hex_f64(oracle):
                raise AssertionError(
                    f"merged sum {merged!r} at {nrank} ranks, fsum {oracle!r}"
                )


def test_shipped_atomic_overlaps_are_rank_identical() -> None:
    wiring = wiring_from_patches()
    if wiring == OLD:
        raise AssertionError(
            "exact patches do not wire the atomic H and S overlaps; "
            "the stock partial sum is rank-dependent"
        )
    hexs = assert_rank_identical(wiring)
    h_hex, s_hex = hexs[0]
    if h_hex != hex_f64(math.fsum(H_TERMS)) or s_hex != hex_f64(math.fsum(S_TERMS)):
        raise AssertionError(f"shipped overlaps {h_hex} {s_hex} are not the once-rounded sums")
    # A missed skip leaves the exact matrix on every rank and summat scales it.
    missed = Wiring(True, True, False, False, True)
    try:
        assert_rank_identical(missed)
    except AssertionError:
        return
    raise AssertionError("dropping both skips left the overlaps bit-identical")


def main() -> int:
    test_old_partial_sum_fails_rank_identity()
    test_exact_sum_matches_one_rounding()
    test_shipped_atomic_overlaps_are_rank_identical()
    print("OK: atomic H and S overlaps are bit-identical at ranks 1, 2, 3, 4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
