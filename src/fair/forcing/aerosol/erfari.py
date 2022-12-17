"""Module for calculating forcing from aerosol-radiation interactions."""

import numpy as np

from ...constants import SPECIES_AXIS


def calculate_erfari_forcing(
    emissions,
    concentration,
    baseline_emissions,
    baseline_concentration,
    forcing_scaling,
    radiative_efficiency,
    emissions_indices,
    concentration_indices,
):
    r"""Calculate effective radiative forcing from aerosol-radiation interactions.

    Parameters
    ----------
    emissions : np.ndarray
        emissions in timestep
    concentration: np.ndarray
        concentrations in timestep
    baseline_emissions : np.ndarray
        baseline emissions to evaluate forcing against
    baseline_concentration : np.ndarray
        baseline concentrations to evaluate forcing against
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    radiative_efficiency : np.ndarray
        radiative efficiency (W m-2 (emission_unit yr-1)-1) of each species.
    emissions_indices : list of int
        provides a mapping of which aerosol species corresponds to which emitted
        species index along the SPECIES_AXIS.
    concentration_indices : list of int
        provides a mapping of which aerosol species corresponds to which
        atmospheric GHG concentration along the SPECIES_AXIS.

    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W m\ :sup:`-2`) from aerosol-radiation interactions
    """
    erf_out = np.ones_like(emissions) * np.nan

    # emissions-driven forcers
    erf_out[..., emissions_indices] = (
        (emissions[..., emissions_indices] - baseline_emissions[..., emissions_indices])
        * radiative_efficiency[..., emissions_indices]
    ) * forcing_scaling[..., emissions_indices]

    # concentration-driven forcers
    erf_out[..., concentration_indices] = (
        (
            concentration[..., concentration_indices]
            - baseline_concentration[..., concentration_indices]
        )
        * radiative_efficiency[..., concentration_indices]
    ) * forcing_scaling[..., concentration_indices]

    # in future we can retain contributions from each species. Will need one ERFari
    # array index for each species so we don't do this here yet.
    return np.nansum(erf_out, axis=SPECIES_AXIS, keepdims=True)
