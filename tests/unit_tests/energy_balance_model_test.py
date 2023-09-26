"""Module for energy balance model tests."""

import os

import numpy as np
import pytest

from fair.energy_balance_model import EnergyBalanceModel

HERE = os.path.dirname(os.path.realpath(__file__))

# this model has the potential to cause problems for scipy's linear algebra
EBM_CAMS_STOCHASTIC = EnergyBalanceModel(
    ocean_heat_capacity=[2.632438882, 9.262194928, 52.92769715],
    ocean_heat_transfer=[1.876253552, 5.15359085, 0.643546006],
    deep_ocean_efficacy=1.285458434,
    gamma_autocorrelation=28.2398724,
    sigma_xi=0.439493317,
    sigma_eta=2.690512385,
    forcing_4co2=8.870602356,
    seed=23,
    stochastic_run=True,
    n_timesteps=5,
)

EBM_CAMS_DETERMINISTIC = EnergyBalanceModel(
    ocean_heat_capacity=[2.632438882, 9.262194928, 52.92769715],
    ocean_heat_transfer=[1.876253552, 5.15359085, 0.643546006],
    deep_ocean_efficacy=1.285458434,
    gamma_autocorrelation=28.2398724,
    forcing_4co2=8.870602356,
    seed=None,
    stochastic_run=False,
    n_timesteps=5,
)


def test_ebm_init_array_mismatch_error():
    with pytest.raises(ValueError):
        EnergyBalanceModel(
            ocean_heat_capacity=[2, 10, 75], ocean_heat_transfer=[1.0, 3.0]
        )


def test_ebm_init_ocean_layers_less_than_two_error():
    with pytest.raises(ValueError):
        EnergyBalanceModel(ocean_heat_capacity=[2], ocean_heat_transfer=[1.0])


def test_ebm_stochastic_d():
    np.testing.assert_array_equal(EBM_CAMS_DETERMINISTIC.stochastic_d, 0)


def test_ebm_emergent_parameters():
    EBM_CAMS_STOCHASTIC.impulse_response()
    EBM_CAMS_STOCHASTIC.emergent_parameters()
    # Generates the test data. Uncomment next lines if you want it.
    # np.savetxt(os.path.join(
    #     HERE, "..", "test_data", "ebm3_cams-csm1-0_timescales.txt"),
    #     EBM_CAMS_STOCHASTIC.timescales
    # )
    # np.savetxt(os.path.join(
    #     HERE, "..", "test_data", "ebm3_cams-csm1-0_response_coefficients.txt"),
    #     EBM_CAMS_STOCHASTIC.response_coefficients
    # )
    # np.savetxt(os.path.join(
    #     HERE, "..", "test_data", "ebm3_cams-csm1-0_ecs_tcr.txt"),
    #     np.array([EBM_CAMS_STOCHASTIC.ecs, EBM_CAMS_STOCHASTIC.tcr])
    # )
    timescales = np.loadtxt(
        os.path.join(HERE, "..", "test_data", "ebm3_cams-csm1-0_timescales.txt")
    )
    response_coefficients = np.loadtxt(
        os.path.join(
            HERE, "..", "test_data", "ebm3_cams-csm1-0_response_coefficients.txt"
        )
    )
    ecs_tcr = np.loadtxt(
        os.path.join(HERE, "..", "test_data", "ebm3_cams-csm1-0_ecs_tcr.txt")
    )
    np.testing.assert_allclose(EBM_CAMS_STOCHASTIC.timescales, timescales)
    np.testing.assert_allclose(
        EBM_CAMS_STOCHASTIC.response_coefficients, response_coefficients
    )
    np.testing.assert_allclose(
        np.array([EBM_CAMS_STOCHASTIC.ecs, EBM_CAMS_STOCHASTIC.tcr]), ecs_tcr
    )


@pytest.mark.filterwarnings("ignore:covariance")
def test_ebm_run():
    EBM_CAMS_STOCHASTIC.add_forcing(np.zeros(5), timestep=1)
    EBM_CAMS_STOCHASTIC.run()
    # Generates the test data. Uncomment next lines if you want it.
    # np.savetxt(os.path.join(
    #     HERE, "test_data", "ebm3_cams-csm1-0_temperature.txt"),
    #     EBM_CAMS_STOCHASTIC.temperature
    # )
    temperature = np.loadtxt(
        os.path.join(HERE, "..", "test_data", "ebm3_cams-csm1-0_temperature.txt")
    )
    # implement a fairly generous absolute tolerance on the temperature differences
    # because the scipy linalg routines seem to change with each version, and if they
    # are out by less than one microkelvin I am sure we can accept this.
    np.testing.assert_allclose(EBM_CAMS_STOCHASTIC.temperature, temperature, atol=1e-6)
