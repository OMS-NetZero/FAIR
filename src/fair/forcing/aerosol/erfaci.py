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

    F = \beta \log \left( 1 + \sum_{i} s_i A_i \right)

    where
        $A_i$ is the emissions or concentration of a specie
        $\beta$ is the scale factor
        $s_i$ describes how sensitive a specie is in contributing to ERFaci

    The calculation is performed for the emissions/concentration of interest and then
    for the baseline. The difference between the two values is the forcing.

    This relationship is a generalisation of Stevens (2015) [1]_. To recover Stevens
    (2015), set $s_i$ to zero for all except SO2, and $s_{\{textrm{SO2}}$ = $1/Q_n$
    where $Q_n$ is the natural emissions source in Stevens (2015).

    Inputs
    ------
    emissions : ndarray
        input emissions
    concentration : ndarray
        input emissions
    baseline_emissions : ndarray
        baseline emissions
    baseline_concentration : ndarray
        baseline concentration
    forcing_scaling : ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    scale : ndarray
        per-species scaling factor to apply to the logarithm
    sensitivity : ndarray
        per-species sensitivity factor for the logarithm
    slcf_indices : list of int
        provides a mapping of which aerosol species corresponds to which emitted
        species index along the SPECIES_AXIS.
    ghg_indices : list of int
        provides a mapping of which aerosol species corresponds to which
        atmospheric GHG concentration along the SPECIES_AXIS.

    Returns
    -------
    effective_radiative_forcing : ndarray
        effective radiative forcing (W/m2) from aerosol-cloud interactions

    References
    ----------
    .. [1] Stevens, B. (2015). Rethinking the Lower Bound on Aerosol Radiative
        Forcing, Journal of Climate, 28(12), 4794-4819.
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
