"""Internal element-symbol parsing shared by Python frontends."""

from __future__ import annotations

import re

__all__ = ["species_numbers"]

_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es "
    "Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
_NUMBER = {symbol: index + 1 for index, symbol in enumerate(_SYMBOLS)}
_FORMULA_PART = re.compile(r"([A-Z][a-z]?)(\d*)")


def species_numbers(text: str) -> list[int]:
    """Return atomic numbers from symbols, formulas, or atomic numbers."""
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
