"""Module for testing fill_from functions."""

import os

import pytest

from fair import FAIR
from fair.io import read_properties

EMISSIONS_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", "test_data", "fill_from"
)
PARAMS_FILE = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "..",
    "..",
    "examples",
    "data",
    "importing-data",
    "calibrated_constrained_parameters.csv",
)


def minimal_problem_def(configs=["one", "two", "three"]):
    fair_obj = FAIR()
    species, properties = read_properties(species=["CO2"])
    properties["CO2"]["input_mode"] = "emissions"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 1753, 1)
    fair_obj.define_scenarios(["test"])
    fair_obj.define_configs(configs)
    fair_obj.allocate()
    return fair_obj


def test_override_defaults():
    f = minimal_problem_def()
    f.fill_from_csv(
        emissions_file=os.path.join(EMISSIONS_PATH, "minimal-emissions.csv")
    )
    f.fill_species_configs()
    f.override_defaults(PARAMS_FILE)

    # "four" is not in the configs file, so should raise an error
    f = minimal_problem_def(configs=["four"])
    f.fill_from_csv(
        emissions_file=os.path.join(EMISSIONS_PATH, "minimal-emissions.csv")
    )
    f.fill_species_configs()
    with pytest.raises(KeyError):
        f.override_defaults(PARAMS_FILE)
