"""Pin the species grammar shared by the i-PI driver and the LAMMPS front end.

``cpmdc._elements.species_numbers`` is documented in its docstring and in
``docs/orgmode/reference/python.org``. These tests hold that grammar fixed
against an independent element table, so a reordered or truncated symbol
list, or a loosened token rule, fails here and not in a trajectory.
"""

from __future__ import annotations

import pytest

from cpmdc import _elements
from cpmdc._elements import species_numbers
from cpmdc.lammps import type_numbers


def test_the_symbol_table_is_the_periodic_table_through_oganesson():
    data = pytest.importorskip("ase.data")
    assert list(_elements._SYMBOLS) == list(data.chemical_symbols[1:119])


def test_every_symbol_and_its_number_name_the_same_element():
    for number, symbol in enumerate(_elements._SYMBOLS, start=1):
        assert species_numbers(symbol) == [number]
        assert species_numbers(str(number)) == [number]


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("O", [8]),
        ("O H H", [8, 1, 1]),
        ("O,H,H", [8, 1, 1]),
        ("O, H ,H", [8, 1, 1]),
        ("  O\tH\nH  ", [8, 1, 1]),
        ("OH2", [8, 1, 1]),
        ("Si3N4", [14, 14, 14, 7, 7, 7, 7]),
        ("8 1 1", [8, 1, 1]),
        ("O 1 H", [8, 1, 1]),
        ("H10", [1] * 10),
        ("Si0 O", [8]),
    ],
)
def test_the_documented_forms_read_in_the_order_written(text, expected):
    assert species_numbers(text) == expected


def test_a_symbol_is_a_capital_with_an_optional_lowercase_letter():
    assert species_numbers("Co") == [27]
    assert species_numbers("CO") == [6, 8]
    assert species_numbers("Hf") == [72]
    assert species_numbers("HF") == [1, 9]
    assert species_numbers("H F") == [1, 9]


def test_a_formula_count_applies_to_the_symbol_before_it():
    assert species_numbers("Fe2O3") == [26, 26, 8, 8, 8]
    assert species_numbers("C2H5OH") == [6, 6, 1, 1, 1, 1, 1, 8, 1]


@pytest.mark.parametrize("text", ["0", "119", "1000", "00"])
def test_an_atomic_number_outside_1_to_118_is_refused(text):
    with pytest.raises(ValueError, match="outside 1..118"):
        species_numbers(text)


@pytest.mark.parametrize("text", ["Xx", "Qq", "Zz2", "D", "Uue"])
def test_an_unknown_element_is_refused_by_name(text):
    with pytest.raises(ValueError, match="unknown element|cannot read species"):
        species_numbers(text)


@pytest.mark.parametrize(
    "text",
    ["o", "h2o", "H-2", "H+", "2H", "0H", "H.5", "H_2", "(OH)2", "O=H", "H2 O!"],
)
def test_a_token_that_is_not_a_symbol_formula_or_number_is_refused(text):
    with pytest.raises(ValueError, match="cannot read species"):
        species_numbers(text)


@pytest.mark.parametrize("text", ["", "   ", ",", " , ,", "H0", "Si0 N0"])
def test_a_string_that_names_no_atoms_is_refused(text):
    with pytest.raises(ValueError, match="at least one atom"):
        species_numbers(text)


def test_a_refused_token_names_the_token_it_could_not_read():
    with pytest.raises(ValueError) as caught:
        species_numbers("O H h2o")
    assert "'h2o'" in str(caught.value)


def test_the_i_pi_driver_exports_the_shared_parser():
    from cpmdc import ipi

    assert ipi.species_numbers is species_numbers
    assert "species_numbers" in ipi.__all__


@pytest.mark.parametrize(
    ("elements", "numbers"),
    [("O H", [8, 1]), (["Co", "C"], [27, 6]), ([14, "N", "7"], [14, 7, 7])],
)
def test_a_lammps_type_list_uses_the_same_grammar(elements, numbers):
    assert type_numbers(elements, len(numbers)).tolist() == numbers


@pytest.mark.parametrize("item", ["OH", "H2", "O,H", "0"])
def test_a_lammps_type_names_exactly_one_element(item):
    with pytest.raises(ValueError):
        type_numbers([item], 1)
