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
    Skeie et al. (2020) [1]_ and Thornhill et al. (2021) [2]_.

    The treatment is identical to FaIR2.0 [3]_.

    Parameters
    ----------
    emissions : ndarray
        emissions in timestep
    concentration: ndarray
        concentrations in timestep
    baseline_emissions : ndarray
        reference, possibly pre-industrial, emissions
    baseline_concentration : ndarray
        reference, possibly pre-industrial, concentrations
    forcing_scaling : ndarray
        scaling of the calculated radiative forcing (e.g. for conversion to
        effective radiative forcing and forcing uncertainty).
    ozone_radiative_efficiency : ndarray
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
    _erf_out : ndarray
        ozone forcing

    References
    ----------
    .. [1] Skeie, R.B., Myhre, G., Hodnebrog, Ø., Cameron-Smith, P.J.,
        Deushi, M., Hegglin, M.I., Horowitz, L.W., Kramer, R.J., Michou, M.,
        Mills, M.J., Olivié, D.J., Connor, F.M., Paynter, D., Samset, B.H.,
        Sellar, A., Shindell, D., Takemura, T., Tilmes, S., Wu, T., 2020.
        Historical total ozone radiative forcing derived from CMIP6 simulations,
        npj Climate and Atmospheric Science, 3, 1–10.

    .. [2] Thornhill, G.D., Collins, W.J., Kramer, R.J., Olivié, D., Skeie,
        R.B., O'Connor, F.M., Abraham, N.L., Checa-Garcia, R., Bauer, S.E.,
        Deushi, M., Emmons, L.K., Forster, P.M., Horowitz, L.W., Johnson, B.,
        Keeble, J., Lamarque, J.-F., Michou, M., Mills, M.J., Mulcahy, J.P.,
        Myhre, G., Nabat, P., Naik, V., Oshima, N., Schulz, M., Smith, C.J.,
        Takemura, T., Tilmes, S., Wu, T., Zeng, G., Zhang, J. (2021). Effective
        radiative forcing from emissions of reactive gases and aerosols – a
        multi-model comparison, Atmospheric Chemistry and  Physics, 21, 853–874

    .. [3] Leach, N.J., Jenkins, S., Nicholls, Z., Smith, C.J., Lynch, J.,
        Cain, M., Walsh, T., Wu, B., Tsutsui, J., Allen, M.R. (2021). FaIRv2.0.0:
        a generalized impulse response model for climate uncertainty and future
        scenario exploration, Geoscientific Model Development, 14, 3007–3036.
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
