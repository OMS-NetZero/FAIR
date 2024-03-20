"""Module for testing fill_from functions."""

import os

import pytest
from fair import FAIR
from fair.io import read_properties

from fair.exceptions import MissingColumnError

HERE = os.path.dirname(os.path.realpath(__file__))


def test__check_csv_raises():
    f = FAIR()
    species = ["CO2"]
    species, properties = read_properties(species=species)
    f.define_species(species, properties)
    f.define_time(1750, 1753, 1)
    f.define_scenarios(["historical"])
    f.define_configs(["UKESM1-0-LL"])
    f.allocate()
    with pytest.raises(MissingColumnError):
        f.fill_from_csv(emissions_file=os.path.join(HERE, "..", "..", "test_data", "fill_from", "variable-not-defined.csv"))