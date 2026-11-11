"""Module for calculating forcing from aerosol-cloud interactions from aviation."""

import numpy as np

from ...constants import SPECIES_AXIS


def calculate_erfaciavi_forcing(
    emissions,
    forcing_scaling,
    linear_radiative_efficiency,
    linear_emissions_indices,
):
    r"""Calculate effective radiative forcing from aerosol-cloud
    interactions from aviation.

    Parameters
    ----------
    emissions : np.ndarray
        emissions in timestep
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    linear_radiative_efficiency : np.ndarray
        radiative efficiency (W m-2 (emission_unit yr-1)-1) of each species.
    linear_emissions_indices : list of int
        provides a mapping of which aerosol species corresponds to which emitted
        species index along the SPECIES_AXIS.
    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W m\ :sup:`-2`) from aerosol-cloud
        interactions from aviation
    """
    erf_out = np.ones_like(emissions) * np.nan

    # we have only the linear component here - the nonlinear interaction is
    # small based on current data
    erf_out[..., linear_emissions_indices] = (
        emissions[..., linear_emissions_indices]
        * linear_radiative_efficiency[..., linear_emissions_indices]
    ) * forcing_scaling


    return np.nansum(erf_out, axis=SPECIES_AXIS, keepdims=True)
