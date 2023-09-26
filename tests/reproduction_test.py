"""Reproduction tests."""

import os

import pandas as pd
import pytest
import xarray as xr

from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties


@pytest.mark.filterwarnings("ignore:numpy.ndarray size changed")
def test_ssp_emissions_cmip6_ebm3_calibrations():
    f = FAIR(ch4_method="thornhill2021")
    f.define_time(1750, 2100, 1)

    scenarios = [
        "ssp119",
        "ssp126",
        "ssp245",
        "ssp370",
        "ssp434",
        "ssp460",
        "ssp534-over",
        "ssp585",
    ]
    f.define_scenarios(scenarios)

    HERE = os.path.dirname(os.path.realpath(__file__))

    df = pd.read_csv(os.path.join(HERE, "test_data", "4xCO2_cummins_ebm3.csv"))
    models = df["model"].unique()
    configs = []

    for imodel, model in enumerate(models):
        for run in df.loc[df["model"] == model, "run"]:
            configs.append(f"{model}_{run}")
    f.define_configs(configs)

    species, properties = read_properties()
    f.define_species(species, properties)

    f.allocate()

    f.fill_species_configs()
    fill(f.species_configs["unperturbed_lifetime"], 10.8537568, specie="CH4")
    fill(f.species_configs["baseline_emissions"], 19.01978312, specie="CH4")
    fill(f.species_configs["baseline_emissions"], 0.08602230754, specie="N2O")

    f.fill_from_rcmip()
    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)

    models = df["model"].unique()

    seed = 1355763

    for config in configs:
        model, run = config.split("_")
        condition = (df["model"] == model) & (df["run"] == run)
        fill(
            f.climate_configs["ocean_heat_capacity"],
            df.loc[condition, "C1":"C3"].values.squeeze(),
            config=config,
        )
        fill(
            f.climate_configs["ocean_heat_transfer"],
            df.loc[condition, "kappa1":"kappa3"].values.squeeze(),
            config=config,
        )
        fill(
            f.climate_configs["deep_ocean_efficacy"],
            df.loc[condition, "epsilon"].values[0],
            config=config,
        )
        fill(
            f.climate_configs["gamma_autocorrelation"],
            df.loc[condition, "gamma"].values[0],
            config=config,
        )
        fill(
            f.climate_configs["sigma_eta"],
            df.loc[condition, "sigma_eta"].values[0],
            config=config,
        )
        fill(
            f.climate_configs["sigma_xi"],
            df.loc[condition, "sigma_xi"].values[0],
            config=config,
        )
        fill(f.climate_configs["stochastic_run"], True, config=config)
        fill(f.climate_configs["use_seed"], True, config=config)
        fill(f.climate_configs["seed"], seed, config=config)

        seed = seed + 399

    f.run(progress=False)

    # These lines write expected test results. If the model is changed, and we mean it,
    # uncomment to make new test results.
    # f.temperature[-1,...,0].to_netcdf(os.path.join(HERE, "test_data",
    # "cmip6_ssp_emissions_run_temperature_2100.nc"))
    # f.concentration[-1,...,2].to_netcdf(os.path.join(HERE, "test_data",
    # "cmip6_ssp_emissions_run_co2_concentration_2100.nc"))
    # f.forcing_sum[-1,...].to_netcdf(os.path.join(HERE, "test_data",
    # "cmip6_ssp_emissions_run_forcing_sum_2100.nc"))

    expected_temperature = xr.open_dataarray(
        os.path.join(HERE, "test_data", "cmip6_ssp_emissions_run_temperature_2100.nc")
    )
    expected_co2_concentration = xr.open_dataarray(
        os.path.join(
            HERE, "test_data", "cmip6_ssp_emissions_run_co2_concentration_2100.nc"
        )
    )
    expected_forcing_sum = xr.open_dataarray(
        os.path.join(HERE, "test_data", "cmip6_ssp_emissions_run_forcing_sum_2100.nc")
    )

    xr.testing.assert_allclose(expected_temperature, f.temperature[-1, ..., 0])
    xr.testing.assert_allclose(expected_co2_concentration, f.concentration[-1, ..., 2])
    xr.testing.assert_allclose(expected_forcing_sum, f.forcing_sum[-1, ...])
