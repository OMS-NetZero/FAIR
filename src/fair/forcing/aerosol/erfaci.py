"""Module for forcing from aerosol-cloud interactions."""

import numpy as np

from ...constants import SPECIES_AXIS


def logsum(
    emissions,
    concentration,
    baseline_emissions,
    baseline_concentration,
    forcing_scaling,
    scale,
    sensitivity,
    slcf_indices,
    ghg_indices,
):
    r"""Calculate effective radiative forcing from aerosol-cloud interactions.

    This uses the relationship to calculate ERFaci as follows

    .. math::

        F = \beta \log \left( 1 + \sum_{i} s_i A_i \right)

    where :math:`A_i` is the emissions or concentration of a specie, :math:`\beta` is
    the scale factor and :math:`s_i` describes how sensitive a specie is in contributing
    to ERFaci.

    The calculation is performed for the emissions/concentration of interest and then
    for the baseline. The difference between the two values is the forcing.

    This relationship is a generalisation of [Stevens2015]_. To recover [Stevens2015]_,
    set :math:`s_i` to zero for all except SO\ :sub:`2`, and :math:`s_{\text{SO2}} =
    1/Q_n` where :math:`Q_n` is the natural emissions source in [Stevens2015]_.

    Parameters
    ----------
    emissions : np.ndarray
        input emissions
    concentration : np.ndarray
        input emissions
    baseline_emissions : np.ndarray
        baseline emissions
    baseline_concentration : np.ndarray
        baseline concentration
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    scale : np.ndarray
        per-species scaling factor to apply to the logarithm
    sensitivity : np.ndarray
        per-species sensitivity factor for the logarithm
    slcf_indices : list of int
        provides a mapping of which aerosol species corresponds to which emitted
        species index along the SPECIES_AXIS.
    ghg_indices : list of int
        provides a mapping of which aerosol species corresponds to which
        atmospheric GHG concentration along the SPECIES_AXIS.

    Returns
    -------
    effective_radiative_forcing : np.ndarray
        effective radiative forcing (W m\ :sup:`-2`) from aerosol-cloud interactions
    """
    radiative_effect = scale * (
        np.log(
            1
            + np.nansum(
                slcf_indices[None, None, None, :] * sensitivity * emissions,
                axis=SPECIES_AXIS,
            )
            + np.nansum(
                ghg_indices[None, None, None, :] * sensitivity * concentration,
                axis=SPECIES_AXIS,
            )
        )
    )
    baseline_radiative_effect = scale * (
        np.log(
            1
            + np.nansum(
                slcf_indices[None, None, None, :] * sensitivity * baseline_emissions,
                axis=SPECIES_AXIS,
            )
            + np.nansum(
                ghg_indices[None, None, None, :] * sensitivity * baseline_concentration,
                axis=SPECIES_AXIS,
            )
        )
    )
    erf_out = (
        radiative_effect[..., None] - baseline_radiative_effect[..., None]
    ) * forcing_scaling
    return erf_out
