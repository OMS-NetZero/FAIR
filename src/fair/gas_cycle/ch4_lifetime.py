"""Methane lifetime definition that is based on multiple species."""

import numpy as np

from ..constants import SPECIES_AXIS


def calculate_alpha_ch4(
    emissions,
    concentration,
    temperature,
    baseline_emissions,
    baseline_concentration,
    ch4_lifetime_chemical_sensitivity,
    ch4_lifetime_temperature_sensitivity,
    emissions_indices,
    concentration_indices,
):
    """Calculate methane lifetime scaling factor.

    Parameters
    ----------
    emissions : np.ndarray
        Emitted species
    concentration : np.ndarray
        Species derived from concentration (including EESC)
    temperature : np.ndarray
        Global mean surface level temperature anomaly
    baseline_emissions : np.ndarray
        reference, possibly pre-industrial, emissions
    baseline_concentration : np.ndarray
        reference, possibly pre-industrial, concentration
    ch4_lifetime_chemical_sensitivity : np.ndarray
        per-species sensitivity of change in CH4 lifetime per unit emission.
    ch4_lifetime_temperature_sensitivity : np.ndarray
        per-species sensitivity of change in CH4 lifetime per unit concentration
        increase.
    emissions_indices : np.ndarray of bool
        Which species indices to use emissions from.
    concentration_indices : np.ndarray of bool
        Which species indices to use concentration from.

    Returns
    -------
    alpha_ch4 : np.ndarray
        Scaling factor to baseline methane lifetime.
    """
    log_lifetime_scaling = (
        np.sum(
            np.log(
                1
                + (
                    emissions[..., emissions_indices]
                    - baseline_emissions[..., emissions_indices]
                )
                * ch4_lifetime_chemical_sensitivity[..., emissions_indices]
            ),
            axis=SPECIES_AXIS,
            keepdims=True,
        )
        + np.sum(
            np.log(
                1
                + (
                    concentration[..., concentration_indices]
                    - baseline_concentration[..., concentration_indices]
                )
                * ch4_lifetime_chemical_sensitivity[..., concentration_indices],
            ),
            axis=SPECIES_AXIS,
            keepdims=True,
        )
        + np.log(1 + temperature * ch4_lifetime_temperature_sensitivity)
    )

    return np.exp(log_lifetime_scaling)
