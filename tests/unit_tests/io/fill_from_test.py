"""Module for testing fill_from functions."""

import logging
import os

import pytest

from fair import FAIR
from fair.exceptions import (
    MetaAfterValueError,
    MissingColumnError,
    NonMonotonicError,
    UnitParseError,
)
from fair.io import read_properties

TEST_DATA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", "test_data", "fill_from"
)


def minimal_problem_def(species=["CO2"]):
    fair_obj = FAIR()
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = "emissions"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 1753, 1)
    fair_obj.define_scenarios(["test"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    return fair_obj


def test__check_csv():
    f = minimal_problem_def()

    with pytest.raises(MissingColumnError):
        f.fill_from_csv(
            emissions_file=os.path.join(TEST_DATA_PATH, "variable-not-defined.csv")
        )

    with pytest.raises(MetaAfterValueError):
        f.fill_from_csv(
            emissions_file=os.path.join(TEST_DATA_PATH, "meta-after-value.csv")
        )

    with pytest.raises(NonMonotonicError):
        f.fill_from_csv(
            emissions_file=os.path.join(TEST_DATA_PATH, "non-monotonic.csv")
        )


def test__parse_unit(caplog):
    caplog.set_level(logging.DEBUG)
    f = minimal_problem_def()

    with pytest.raises(UnitParseError):
        f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "bad-unit.csv"))

    with pytest.raises(UnitParseError):
        f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "bad-prefix.csv"))

    with pytest.raises(UnitParseError):
        f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "bad-time.csv"))


def test__non_default_specie(caplog):
    f = FAIR()
    species = ["HFC-152", "Hydrogen"]
    properties = {}
    properties["HFC-152"] = {}
    properties["HFC-152"]["input_mode"] = "emissions"
    properties["HFC-152"]["greenhouse_gas"] = True
    properties["HFC-152"]["type"] = "f-gas"
    properties["HFC-152"]["aerosol_chemistry_from_emissions"] = False
    properties["HFC-152"]["aerosol_chemistry_from_concentration"] = False
    properties["HFC-152"]["aerosol_chemistry_from_emissions"] = False
    properties["Hydrogen"] = {}
    properties["Hydrogen"]["input_mode"] = "emissions"
    properties["Hydrogen"]["greenhouse_gas"] = False
    properties["Hydrogen"]["type"] = "other slcf"
    properties["Hydrogen"]["aerosol_chemistry_from_emissions"] = False
    properties["Hydrogen"]["aerosol_chemistry_from_concentration"] = False
    properties["Hydrogen"]["aerosol_chemistry_from_emissions"] = False
    f.define_species(species, properties)
    f.define_time(1750, 1753, 1)
    f.define_scenarios(["test"])
    f.define_configs(["UKESM1-0-LL"])
    f.allocate()
    f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "new-specie.csv"))
    assert "HFC152 is not in fair's default list" in caplog.text
    assert "H2 is not in fair's default list" in caplog.text


def test_fill_from_csv(caplog):
    caplog.set_level(logging.DEBUG)
    f = minimal_problem_def()
    f.fill_from_csv(
        emissions_file=os.path.join(TEST_DATA_PATH, "minimal-emissions.csv")
    )
    assert (
        "The last time in the emissions file (1752) is earlier than the last time "
        "in the problem definition (1752.5)" in caplog.text
    )
