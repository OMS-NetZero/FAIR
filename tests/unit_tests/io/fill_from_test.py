"""Module for testing fill_from functions."""

import os

import numpy as np
import pytest

from fair import FAIR
from fair.exceptions import (
    DuplicateScenarioError,
    MetaAfterValueError,
    MissingColumnError,
    MissingDataError,
    NonMonotonicError,
    UnitParseError,
)
from fair.io import read_properties

TEST_DATA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "..", "test_data", "fill_from"
)


def minimal_problem_def(input_mode="emissions", species=["CO2"]):
    fair_obj = FAIR()
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = input_mode
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 1753, 1)
    fair_obj.define_scenarios(["test"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    return fair_obj


def minimal_ghg_run(timestep=270, stochastic_run=False, seed=37):
    fair_obj = FAIR()
    species = ["CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = "concentration"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, timestep)
    fair_obj.define_scenarios(["historical"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    fair_obj.climate_configs["ocean_heat_capacity"][0, :] = np.array(
        [2.917300055, 11.28317472, 73.2487238]
    )
    fair_obj.climate_configs["ocean_heat_transfer"][0, :] = np.array(
        [0.65576633, 2.597877675, 0.612933889]
    )
    fair_obj.climate_configs["deep_ocean_efficacy"][0] = 1.133708775
    fair_obj.climate_configs["gamma_autocorrelation"][0] = 3.548407499
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(timestep)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(timestep)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = stochastic_run
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = seed
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [277, 731, 270]
    fair_obj.species_configs["forcing_reference_concentration"][0, :] = [277, 731, 270]
    fair_obj.concentration[0, 0, 0, :] = [277, 731, 270]
    fair_obj.concentration[1:, 0, 0, :] = [410, 1900, 325]
    fair_obj.forcing[0, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
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
    properties["Hydrogen"] = {}
    properties["Hydrogen"]["input_mode"] = "emissions"
    properties["Hydrogen"]["greenhouse_gas"] = False
    properties["Hydrogen"]["type"] = "other slcf"
    properties["Hydrogen"]["aerosol_chemistry_from_emissions"] = False
    properties["Hydrogen"]["aerosol_chemistry_from_concentration"] = False
    f.define_species(species, properties)
    f.define_time(1750, 1753, 1)
    f.define_scenarios(["test"])
    f.define_configs(["UKESM1-0-LL"])
    f.allocate()
    f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "new-specie.csv"))


def test__parse_unit():
    f = minimal_problem_def()

    with pytest.raises(UnitParseError):
        f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "bad-unit.csv"))

    with pytest.raises(UnitParseError):
        f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "bad-prefix.csv"))

    with pytest.raises(UnitParseError):
        f.fill_from_csv(emissions_file=os.path.join(TEST_DATA_PATH, "bad-time.csv"))


def test__concentration_unit_convert():
    f = FAIR()
    species = ["CO2"]
    properties = {}
    properties["CO2"] = {}
    properties["CO2"]["input_mode"] = "concentration"
    properties["CO2"]["greenhouse_gas"] = True
    properties["CO2"]["type"] = "co2"
    properties["CO2"]["aerosol_chemistry_from_emissions"] = False
    properties["CO2"]["aerosol_chemistry_from_concentration"] = False
    f.define_species(species, properties)
    f.define_time(1750, 1753, 1)
    f.define_scenarios(["test"])
    f.define_configs(["UKESM1-0-LL"])
    f.allocate()
    with pytest.raises(UnitParseError):
        f.fill_from_csv(
            concentration_file=os.path.join(TEST_DATA_PATH, "bad-mixing-ratio.csv")
        )

    f = FAIR()
    species = ["PF3"]
    properties = {}
    properties["PF3"] = {}
    properties["PF3"]["input_mode"] = "concentration"
    properties["PF3"]["greenhouse_gas"] = True
    properties["PF3"]["type"] = "f-gas"
    properties["PF3"]["aerosol_chemistry_from_emissions"] = False
    properties["PF3"]["aerosol_chemistry_from_concentration"] = False
    f.define_species(species, properties)
    f.define_time(1750, 1753, 1)
    f.define_scenarios(["test"])
    f.define_configs(["UKESM1-0-LL"])
    f.allocate()
    f.fill_from_csv(
        concentration_file=os.path.join(TEST_DATA_PATH, "new-concentration-specie.csv")
    )


# this one unfinished
def test_fill_from_csv():
    f = minimal_problem_def(input_mode="concentration")
    f.fill_from_csv(
        concentration_file=os.path.join(TEST_DATA_PATH, "minimal-concentration.csv")
    )

    f = minimal_problem_def(input_mode="forcing", species=["Solar", "Volcanic"])
    f.fill_from_csv(forcing_file=os.path.join(TEST_DATA_PATH, "minimal-forcing.csv"))

    f = minimal_problem_def()
    f.fill_from_csv(
        emissions_file=os.path.join(TEST_DATA_PATH, "minimal-emissions.csv")
    )
    with pytest.raises(DuplicateScenarioError):
        f.fill_from_csv(
            emissions_file=os.path.join(TEST_DATA_PATH, "duplicate-emissions.csv")
        )

    f = FAIR()
    species = ["CO2", "PF3", "Volcanic"]
    properties = {}
    properties["CO2"] = {}
    properties["CO2"]["input_mode"] = "emissions"
    properties["CO2"]["greenhouse_gas"] = True
    properties["CO2"]["type"] = "co2"
    properties["CO2"]["aerosol_chemistry_from_emissions"] = False
    properties["CO2"]["aerosol_chemistry_from_concentration"] = False
    properties["PF3"] = {}
    properties["PF3"]["input_mode"] = "concentration"
    properties["PF3"]["greenhouse_gas"] = True
    properties["PF3"]["type"] = "f-gas"
    properties["PF3"]["aerosol_chemistry_from_emissions"] = False
    properties["PF3"]["aerosol_chemistry_from_concentration"] = False
    properties["Volcanic"] = {}
    properties["Volcanic"]["input_mode"] = "forcing"
    properties["Volcanic"]["greenhouse_gas"] = False
    properties["Volcanic"]["type"] = "volcanic"
    properties["Volcanic"]["aerosol_chemistry_from_emissions"] = False
    properties["Volcanic"]["aerosol_chemistry_from_concentration"] = False
    f.define_species(species, properties)
    f.define_time(1750, 1753, 1)
    f.define_scenarios(["test"])
    f.define_configs(["UKESM1-0-LL"])
    f.allocate()
    f.fill_from_csv(
        emissions_file=os.path.join(TEST_DATA_PATH, "minimal-emissions.csv"),
        concentration_file=os.path.join(TEST_DATA_PATH, "new-concentration-specie.csv"),
        forcing_file=os.path.join(TEST_DATA_PATH, "minimal-forcing.csv"),
    )


def test_from_rcmip():
    ftest = minimal_ghg_run()
    ftest.fill_from_rcmip()


def test_fill_from_rcmip_missing_concentration_data():
    ftest = minimal_ghg_run()
    ftest.scenarios = ["ADVANCE"]
    with pytest.raises(MissingDataError):
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
    with pytest.raises(MissingDataError):
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
    with pytest.raises(MissingDataError):
        fair_obj.fill_from_rcmip()
