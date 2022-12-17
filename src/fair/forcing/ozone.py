"""Module for ozone forcing."""

import numpy as np

from ..constants import SPECIES_AXIS


def thornhill2021(
    emissions,
    concentration,
    baseline_emissions,
    baseline_concentration,
    forcing_scaling,
    ozone_radiative_efficiency,
    slcf_indices,
    ghg_indices,
):
    """Determine ozone effective radiative forcing.

    Calculates total ozone forcing from precursor emissions and
    concentrations based on AerChemMIP and CMIP6 Historical behaviour in
    [Skeie2020]_ and [Thornhill2021a]_.

    The treatment is identical to FaIR2.0 [Leach2021]_.

    Parameters
    ----------
    emissions : np.ndarray
        emissions in timestep
    concentration: np.ndarray
        concentrations in timestep
    baseline_emissions : np.ndarray
        reference, possibly pre-industrial, emissions
    baseline_concentration : np.ndarray
        reference, possibly pre-industrial, concentrations
    forcing_scaling : np.ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    ozone_radiative_efficiency : np.ndarray
        the radiative efficiency at which ozone precursor emissions or
        concentrations are converted to ozone radiative forcing. The unit is
        W m-2 (<native emissions or concentration unit>)-1, where the
        emissions unit is usually Mt/yr for short-lived forcers, ppb for CH4
        and N2O concentrations, and ppt for halogenated species. Note this is
        not the same radiative efficiency that is used in converting GHG concentrations
        to forcing.
    slcf_indices : list of int
        provides a mapping of which aerosol species corresponds to which emitted
        species index along the SPECIES_AXIS.
    ghg_indices : list of int
        provides a mapping of which aerosol species corresponds to which
        atmospheric GHG concentration along the SPECIES_AXIS.

    Returns
    -------
    _erf_out : np.ndarray
        ozone forcing
    """
    array_shape = emissions.shape
    n_timesteps, n_scenarios, n_configs, n_species = array_shape

    # revisit this if we ever want to dump out intermediate calculations like the
    # feedback strength.
    _erf = np.ones((n_timesteps, n_scenarios, n_configs, 2)) * np.nan

    # GHGs, with a concentration-given ozone radiative_efficiency, including EESC
    _erf[:, :, :, 0] = np.sum(
        (
            concentration[:, :, :, ghg_indices]
            - baseline_concentration[:, :, :, ghg_indices]
        )
        * ozone_radiative_efficiency[:, :, :, ghg_indices]
        * forcing_scaling,
        axis=SPECIES_AXIS,
    )

    # Emissions-based precursors
    _erf[:, :, :, 1] = np.sum(
        (emissions[:, :, :, slcf_indices] - baseline_emissions[:, :, :, slcf_indices])
        * ozone_radiative_efficiency[:, :, :, slcf_indices]
        * forcing_scaling,
        axis=SPECIES_AXIS,
    )

    # Temperature feedback has been moved later in the code

    erf_out = np.sum(_erf, axis=SPECIES_AXIS, keepdims=True)
    return erf_out
