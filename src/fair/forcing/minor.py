"""Module for a generic linear emissions or concentration to forcing calculation."""

import numpy as np

from ..constants import SPECIES_AXIS


def calculate_linear_forcing(
    driver,
    baseline_driver,
    forcing_scaling,
    radiative_efficiency,
):
    """
    Calculate effective radiative forcing from a linear relationship.

    Inputs
    ------
    driver : ndarray
        input emissions, concentration or forcing
    baseline_driver : ndarray
        baseline, possibly pre-industrial, emissions, concentration or forcing.
    forcing_scaling : ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : ndarray
        radiative efficiency (W m-2 (<driver unit>)-1) of each species.

    Returns
    -------
    erf_out : ndarray
        effective radiative forcing (W/m2)
    """
    erf_out = np.nansum(
        ((driver - baseline_driver) * radiative_efficiency) * forcing_scaling,
        axis=SPECIES_AXIS,
        keepdims=True,
    )
    return erf_out
