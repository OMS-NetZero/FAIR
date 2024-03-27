"""Module for io tests."""

from fair.io import read_properties


def test_read_properties():
    read_properties(species=["CO2"])
