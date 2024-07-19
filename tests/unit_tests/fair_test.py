"""Module for unit test of fair."""

import os
import tempfile
import warnings

import numpy as np
import pytest

from fair import FAIR
from fair.forcing.ghg import meinshausen2020
from fair.interface import fill
from fair.io import read_properties

f = FAIR()
HERE = os.path.dirname(os.path.realpath(__file__))


def minimal_empty_run(mode="concentration", timestep=270, temperature_prescribed=False):
    fair_obj = FAIR(temperature_prescribed=temperature_prescribed)
    species = ["CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    for specie in species:
        properties[specie]["input_mode"] = mode
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, timestep)
    fair_obj.define_scenarios(["historical"])
    fair_obj.define_configs(["UKESM1-0-LL"])
    fair_obj.allocate()
    return fair_obj


def minimal_fair_run(
    stochastic_run=False,
    seed=37,
    mode="concentration",
    timestep=270,
    temperature_prescribed=False,
):
    fair_obj = minimal_empty_run(
        mode=mode, timestep=timestep, temperature_prescribed=temperature_prescribed
    )
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


def test_ghg_method():
    f.ghg_method = "LEACH2021"
    assert f.ghg_method == "leach2021"


def test_invalid_ghg_method():
    with pytest.raises(ValueError):
        f.ghg_method = "Quaas2022"


def test_ch4_method():
    f.ch4_method = "LEACH2021"
    assert f.ch4_method == "leach2021"


def test_invalid_ch4_method():
    with pytest.raises(ValueError):
        f.ch4_method = "Quaas2022"


def test_specie_not_in_properties():
    species, properties = read_properties()
    with pytest.raises(ValueError):
        f.define_species(["Kryponite"], properties)


def test_species_invalid_type():
    species, properties = read_properties()
    properties["N2O"]["type"] = "laughing gas"
    with pytest.raises(ValueError):
        f.define_species(species, properties)


def test_species_invalid_input_mode():
    species, properties = read_properties()
    properties["N2O"]["input_mode"] = "inhaled"
    with pytest.raises(ValueError):
        f.define_species(species, properties)


def test_duplicate_species():
    species, properties = read_properties()
    properties["CO2 FFI"]["type"] = "co2"
    properties["CO2 AFOLU"]["type"] = "co2"
    with pytest.raises(ValueError):
        f.define_species(species, properties)


def test_allocate_before_definitions():
    with pytest.raises(AttributeError):
        f.allocate()


# If there is a numerics issue with scipy.sparse.linalg, this is often the first test
# to fail (for reasons unknown) and may be unrelated to the GHG formulae.
def test_ghg_routines():
    EXPECTED_RESULTS = {
        "myhre1998": np.array([2.2028445, 0.44916397, 0.19313647]),
        "etminan2016": np.array([2.21891894, 0.55302311, 0.18624564]),
        "meinshausen2020": np.array([2.1849852, 0.55574659, 0.18577101]),
        "leach2021": np.array([2.20722625, 0.54091863, 0.18102735]),
    }
    for method, results in EXPECTED_RESULTS.items():
        ftest = minimal_fair_run()
        ftest.ghg_method = method
        ftest.run(progress=False)
        forcing_out = ftest.forcing.squeeze()
        np.testing.assert_allclose(forcing_out[1, ...], results)


def test_ghg_forcing_offset():
    EXPECTED_RESULTS = {
        "no_offset": np.array([2.1849852, 0.55574659, 0.18577101]),
        "offset": np.array([2.1849852, 1.7776573, 2.09247788]),
    }
    # test that providing the offset gives the same results as not providing it.
    ftest = minimal_fair_run()
    ftest.ghg_forcing_offset = meinshausen2020(
        np.array([277, 731, 270]).reshape((1, 1, 1, 3)),
        np.array([277, 731, 270]).reshape((1, 1, 1, 3)),
        np.array([1, 1, 1]).reshape((1, 1, 1, 3)),
        np.ones((1, 1, 1, 3)),
        0,
        1,
        2,
        [],
    )
    ftest.run(progress=False)
    forcing_out = ftest.forcing.squeeze()
    np.testing.assert_allclose(forcing_out[1, ...], EXPECTED_RESULTS["no_offset"])

    # now check the results differ if the user-specified offset is different.
    ftest = minimal_fair_run()
    ftest.ghg_forcing_offset = meinshausen2020(
        np.array([277, 0, 0]).reshape((1, 1, 1, 3)),
        np.array([277, 731, 270]).reshape((1, 1, 1, 3)),
        np.array([1, 1, 1]).reshape((1, 1, 1, 3)),
        np.ones((1, 1, 1, 3)),
        0,
        1,
        2,
        [],
    )
    ftest.run(progress=False)
    forcing_out = ftest.forcing.squeeze()
    np.testing.assert_allclose(forcing_out[1, ...], EXPECTED_RESULTS["offset"])


def test_calculate_iirf0():
    EXPECTED_RESULTS = np.array([52.35538747, 8.2499551, 65.44969575])
    ftest = minimal_fair_run()
    ftest.calculate_iirf0()
    np.testing.assert_allclose(
        np.squeeze(ftest.species_configs["iirf_0"]), EXPECTED_RESULTS
    )


def test_calculate_g():
    EXPECTED_RESULTS = {
        "g0": np.array([0.01017828826538349, 0.36785516988915923, 0.07675558835522626]),
        "g1": np.array([11.412622431258765, 8.24941081407049, 25.495288175200994]),
    }
    ftest = minimal_fair_run()
    ftest.calculate_g()
    for variable, results in EXPECTED_RESULTS.items():
        np.testing.assert_allclose(np.squeeze(ftest.species_configs[variable]), results)


def test__make_ebms_climate_configs_nan():
    ftest = minimal_fair_run()
    ftest.climate_configs["ocean_heat_transfer"][0, :] = np.nan
    with pytest.raises(ValueError):
        ftest._make_ebms()


def test__make_ebms_stochastic_climate_configs_nan():
    ftest = minimal_fair_run(stochastic_run=True)
    ftest.climate_configs["sigma_eta"][0] = np.nan
    with pytest.raises(ValueError):
        ftest._make_ebms()


def test__check_properties_raise_if_nan():
    for mode in ["emissions", "concentration", "forcing"]:
        ftest = minimal_empty_run(mode=mode)
        with pytest.raises(ValueError):
            ftest._check_properties()
    ftest = minimal_fair_run(temperature_prescribed=True)
    ftest.temperature[:] = np.nan
    with pytest.raises(ValueError):
        ftest._check_properties()


def test_run_runtime_warning():
    # need to run stochastic (to trigger the problem in the first place) and at a
    # big enough time step to trigger the warning but not too big, else the matrix
    # really is wrong. Doesn't seem the most robust test to future scipy whims.
    ftest = minimal_fair_run(stochastic_run=True, timestep=27)
    # I want a warning (Idlewild; 2005)
    with pytest.warns(RuntimeWarning):
        ftest.run(suppress_warnings=False)
    # I don't want a warning
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ftest.run(suppress_warnings=True)
        ftest.run()


def test_to_netcdf():
    ftest = minimal_fair_run()
    with tempfile.TemporaryFile() as tf:
        ftest.to_netcdf(tf)


def test_fill_species_configs_aci_skips():
    f.define_time(2000, 2003, 1)
    f.define_configs(["high"])
    f.define_scenarios(["abrupt"])
    species_properties_filepath = os.path.join(
        HERE, "..", "test_data", "species_configs_properties_ch4.csv"
    )
    species, properties = read_properties(species_properties_filepath)
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(species_properties_filepath)


def test_fill_species_configs_ch4_skips():
    f.define_time(2000, 2003, 1)
    f.define_configs(["high"])
    f.define_scenarios(["abrupt"])
    species_properties_filepath = os.path.join(
        HERE, "..", "test_data", "species_configs_properties_aci.csv"
    )
    species, properties = read_properties(species_properties_filepath)
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(species_properties_filepath)


def test_ghg_method_raise():
    fair_obj = FAIR(ghg_method="Meinshausen2020")
    species = ["CO2 FFI", "CO2 AFOLU", "CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [
        np.nan,
        np.nan,
        277,
        731,
        270,
    ]
    fair_obj.species_configs["forcing_reference_concentration"][:] = np.nan
    fair_obj.emissions[:, 0, 0, :] = 0
    fair_obj.forcing[0, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_co2_run_mode_incompatible_raise():
    fair_obj = FAIR()
    species = ["CO2"]
    species, properties = read_properties(species=species)
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [277]
    fair_obj.species_configs["forcing_reference_concentration"][0, :] = [277]
    fair_obj.emissions[:, 0, 0, :] = 0
    fair_obj.forcing[0, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_landuse_run_mode_incompatible_raise():
    fair_obj = FAIR()
    species = ["Land use"]
    species, properties = read_properties(species=species)
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.forcing[0, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_thornhill_raises():
    fair_obj = FAIR(ch4_method="thornhill2021")
    species = ["Equivalent effective stratospheric chlorine", "Ozone"]
    species, properties = read_properties(species=species)
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [277]
    fair_obj.species_configs["forcing_reference_concentration"][0, :] = [277]
    fair_obj.forcing[0, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_leach_raises():
    fair_obj = FAIR(ch4_method="leach2021")
    species = ["Equivalent effective stratospheric chlorine", "Ozone"]
    species, properties = read_properties(species=species)
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [277]
    fair_obj.species_configs["forcing_reference_concentration"][0, :] = [277]
    fair_obj.forcing[0, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_meinshausen_raises():
    fair_obj = FAIR(ghg_method="Meinshausen2020")
    species = ["CO2", "CH4", "N2O"]
    species, properties = read_properties(species=species)
    properties["CO2"]["input_mode"] = "forcing"
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [277, 731, 270]
    fair_obj.species_configs["forcing_reference_concentration"][0, :] = [277, 731, 270]
    fair_obj.emissions[:, 0, 0, :] = 0
    fair_obj.forcing[:, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_myhre_raises():
    fair_obj = FAIR(ghg_method="Myhre1998")
    species = ["CH4"]
    species, properties = read_properties(species=species)
    fair_obj.define_species(species, properties)
    fair_obj.define_time(1750, 2020, 27)
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
    fair_obj.climate_configs["sigma_xi"][0] = 0.439126403 / np.sqrt(27)
    fair_obj.climate_configs["sigma_eta"][0] = 0.497441140 / np.sqrt(27)
    fair_obj.climate_configs["forcing_4co2"][0] = 7.378788155
    fair_obj.climate_configs["stochastic_run"][0] = False
    fair_obj.climate_configs["use_seed"][0] = True
    fair_obj.climate_configs["seed"][0] = 0
    fair_obj.fill_species_configs()
    fair_obj.species_configs["baseline_concentration"][0, :] = [731]
    fair_obj.species_configs["forcing_reference_concentration"][0, :] = [731]
    fair_obj.emissions[:, 0, 0, :] = 0
    fair_obj.forcing[:, 0, 0, :] = 0
    fair_obj.temperature[0, 0, 0, :] = 0
    fair_obj.cumulative_emissions[0, 0, 0, :] = 0
    fair_obj.airborne_emissions[0, 0, 0, :] = 0
    with pytest.raises(ValueError):
        fair_obj.run()


def test_run_forcing_driven():
    ftest = minimal_fair_run(mode="forcing")
    fill(ftest.forcing, 0)
    ftest.run()


def test_run_temperature_prescribed():
    ftest = minimal_fair_run(temperature_prescribed=True)
    fill(ftest.temperature, 0)
    ftest.run()
