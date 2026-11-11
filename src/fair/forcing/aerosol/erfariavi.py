"""Module for calculating forcing from aerosol-radiation interactions from aviation."""

import numpy as np

from ...constants import SPECIES_AXIS


def calculate_erfariavi_forcing(
    emissions,
    forcing_scaling,
    linear_radiative_efficiency,
    nonlinear_radiative_efficiency,
    linear_emissions_indices,
    nonlinear_emissions_indices,
):
    r"""Calculate effective radiative forcing from aerosol-radiation
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
    nonlinear_radiative_efficiency : np.ndarray
        radiative efficiency species multiplication term.
    nonlinear_emissions_indices : list of int
        provides the set of aerosol species which contribute to the nonlinear forcing
    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W m\ :sup:`-2`) from aerosol-radiation
        interactions from aviation
    """
    erf_out = np.ones_like(emissions) * np.nan

    # we have the linear component
    erf_out[..., linear_emissions_indices] = (
        emissions[..., linear_emissions_indices]
        * linear_radiative_efficiency[..., linear_emissions_indices]
    ) * forcing_scaling

    # then we add the nonlinear term
    # we add this (arbitrarily) to the first nonlinear specie as we sum anyway,
    # and it's not directly attributable to any specie
    nonlinear_idx = np.where(nonlinear_emissions_indices)[0]

    if len(nonlinear_idx) > 0:
        # first the multiplication of emissions
        interaction = np.prod(
            emissions[..., nonlinear_emissions_indices],
            axis=SPECIES_AXIS,
            keepdims=True,
        )

        # then the sum of the coefficients - they are zero by default and then
        # filled in the call
        interaction_coeff = np.sum(
            nonlinear_radiative_efficiency[..., nonlinear_emissions_indices],
            axis=SPECIES_AXIS,
            keepdims=True,
        )

        # need to be sure if the linear term is nan for this specie we don't
        # lose the nonlinear term from "nan + x = nan"
        erf_out[..., nonlinear_idx[0]] = (
            np.nan_to_num(erf_out[..., nonlinear_idx[0]])
            + (interaction * interaction_coeff * forcing_scaling)[..., 0]
        )

    return np.nansum(erf_out, axis=SPECIES_AXIS, keepdims=True)
