"""Module for io tests."""

import os

import pytest

from fair import FAIR
from fair.exceptions import FromCsvError
from fair.io import read_properties
from fair_test import minimal_ghg_run

HERE = os.path.dirname(os.path.realpath(__file__))

def minimal_init(startyear=2010, endyear=2030):
    fair_obj = FAIR()
    species, properties = read_properties()
    fair_obj.define_species(species, properties)
    fair_obj.define_time(startyear, endyear, 1)
    fair_obj.define_scenarios(["charlie", "delta"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    return fair_obj


def test_read_properties():
    read_properties(species=["CO2"])


def test_fill_from_csv_wrong_columns():
    ftest = minimal_init(endyear=2012)
    with pytest.raises(FromCsvError):
        ftest.fill_from_csv(
            os.path.join(HERE, "..", "test_data", "emissions_pyam_style_wrong_columns.csv"), 
            "emissions"
        )


def test_fill_from_csv_time_out_of_range():
    ftest = minimal_init(startyear=2000)
    with pytest.raises(FromCsvError):
        ftest.fill_from_csv(
            os.path.join(HERE, "..", "test_data", "emissions_pyam_style.csv"), 
            "emissions"
        )


def test_fill_from_csv():
    ftest = minimal_init()
    ftest.fill_from_csv(
        os.path.join(HERE, "..", "test_data", "emissions_pyam_style.csv"), 
        "emissions"
    )


def test_from_rcmip():
    ftest = minimal_ghg_run()
    ftest.fill_from_rcmip()


def test_fill_from_rcmip_missing_concentration_data():
    ftest = minimal_ghg_run()
    ftest.scenarios = ["ADVANCE"]
    with pytest.raises(ValueError):
        ftest.fill_from_rcmip()


def test_fill_from_rcmip_missing_emissions_data():
    fair_obj = FAIR()
    species = ["CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = "emissions"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 270)
    fair_obj.define_scenarios(["ADVANCE"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    with pytest.raises(ValueError):
        fair_obj.fill_from_rcmip()


def test_fill_from_rcmip_missing_forcing_data():
    fair_obj = FAIR()
    species = ["CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = "forcing"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 270)
    fair_obj.define_scenarios(["ADVANCE"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    with pytest.raises(ValueError):
        fair_obj.fill_from_rcmip()
