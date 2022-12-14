"""Module containing gas cycle functions."""

import warnings

import numpy as np


def calculate_alpha(
    airborne_emissions,
    cumulative_emissions,
    g0,
    g1,
    iirf_0,
    iirf_airborne,
    iirf_temperature,
    iirf_uptake,
    temperature,
    iirf_max,
):
    """
    Calculate greenhouse-gas time constant scaling factor.

    Parameters
    ----------
    airborne_emissions : ndarray
        Total cumulative emissions remaining in the atmosphere.
    cumulative_emissions : ndarray
        Total cumulative emissions.
    g0 : ndarray
        control parameter: see Leach et al. (2021)
    g1 : ndarray
        control parameter: see Leach et al. (2021)
    iirf_0 : ndarray
        baseline time-integrated airborne fraction.
    iirf_airborne : ndarray
        sensitivity of time-integrated airborne fraction with airborne
        emissions.
    iirf_temperature : ndarray
        sensitivity of time-integrated airborne fraction with temperature
        anomaly.
    iirf_uptake : ndarray
        sensitivity of time-integrated airborne fraction to the amount of gas
        taken up in the Earth system (cumulative minus airborne).
    temperature : ndarray or float
        K temperature anomaly.
    iirf_max : float
        maximum allowable value of time-integrated airborne fraction

    Notes
    -----
    Where array input is taken, the arrays always have the dimensions of
    (scenario, species, time, gas_box). Dimensionality can be 1, but we
    retain the singleton dimension in order to preserve clarity of
    calculation and speed.

    Returns
    -------
    alpha : float
        scaling factor for lifetimes
    """
    iirf = (
        iirf_0
        + iirf_uptake * (cumulative_emissions - airborne_emissions)
        + iirf_temperature * temperature
        + iirf_airborne * airborne_emissions
    )
    iirf = (iirf > iirf_max) * iirf_max + iirf * (iirf < iirf_max)
    # overflow and invalid value errors occur with very large and small values
    # in the exponential. This happens with very long lifetime GHGs, in practice
    # only CF4. Usually these GHGs don't have a temperature dependence on IIRF
    # but even if they did the lifetimes are so long that it is unlikely to have
    # a measurable effect, so we just set alpha to 1.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        alpha = g0 * np.exp(iirf / g1)
        alpha[np.isnan(alpha)] = 1

    return alpha
