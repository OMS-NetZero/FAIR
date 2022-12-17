"""Module for a generic linear emissions or concentration to forcing calculation."""

import numpy as np

from ..constants import SPECIES_AXIS


def calculate_linear_forcing(
    driver,
    baseline_driver,
    forcing_scaling,
    radiative_efficiency,
):
    r"""
    Calculate effective radiative forcing from a linear relationship.

    Parameters
    ----------
    driver : np.ndarray
        input emissions, concentration or forcing
    baseline_driver : np.ndarray
        baseline, possibly pre-industrial, emissions, concentration or forcing.
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : np.ndarray
        radiative efficiency (W m\ :sup:`-2` (<driver unit>)\ :sup:`-1`) of each
        species.

    Returns
    -------
    erf_out : np.ndarray
        effective radiative forcing (W m\ :sup:`-2`)
    """
    erf_out = np.nansum(
        ((driver - baseline_driver) * radiative_efficiency) * forcing_scaling,
        axis=SPECIES_AXIS,
        keepdims=True,
    )
    return erf_out
