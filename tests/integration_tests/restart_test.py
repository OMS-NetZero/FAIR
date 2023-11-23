"""Module for integration test of fair."""

import copy

import numpy as np

from fair import FAIR
from fair.interface import initialise
from fair.io import read_properties


def minimal_ghg_run(start=1750, end=2020, timestep=27, stochastic_run=False, seed=37):
    fair_obj = FAIR()
    species = ["CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = "concentration"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(start, end, timestep)
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


def test_ohc_restart():
    # the whole restart mechanism is clunky right now and needs improving
    # I only care about OHC - the rest of this will give silly output
    fhist = minimal_ghg_run(stochastic_run=False, start=1750, end=2020, timestep=27)
    fhist.run()
    ohc_out = copy.deepcopy(fhist.ocean_heat_content_change)[-1, ...]
    ffut = minimal_ghg_run(stochastic_run=False, start=2020, end=2101, timestep=27)
    initialise(ffut.ocean_heat_content_change, ohc_out)
    ffut.run()
    # check if first OHC timebound in restart is same as last in historical
    np.testing.assert_allclose(ffut.ocean_heat_content_change[0, ...], ohc_out)
