"""Module for io tests."""

import pytest

from fair import FAIR
from fair.io import read_properties
from fair_test import minimal_ghg_run

def test_read_properties():
    read_properties(species=["CO2"])


def test_fill_from_csv():
    assert 0==1


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
