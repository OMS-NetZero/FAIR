from __future__ import division

import numpy as np
from ..constants import cl_atoms, br_atoms, fracrel


def thornhill_skeie(
        emissions,
        concentrations,
        temperature=0,
        feedback=-0.037,
        beta=np.array([2.33379720e-04,  1.27179106e-03, -6.69347820e-05,
                       1.14647701e-04,  5.14366051e-12,  3.78354423e-03]),
        emissions_pi=np.zeros(40),
        concentrations_pi=np.zeros(31),
    ):
    """Calculates total ozone forcing from precursor emissions and
    concentrations based on AerChemMIP and CMIP6 Historical behaviour

    Skeie et al. (2020)
    Thornhill et al. (2021)

    Unlike Stevenson CMIP5 no distinction is made for tropospheric and
    stratospheric.

    With this formulation, ozone forcing depends on concentrations of
    CH4, N2O, ozone-depleting halogens, and emissions of CO, NVMOC and NOx,
    but any combination of emissions and concentrations are allowed.

    Inputs:
        emissions: (nt x 40) numpy array in FaIR default units
        concentrations: (nt x 31) numpy array of GHGs in FaIR default units
        temperature: global mean surface temperature (for feedback)
        feedback: temperature feedback on ozone forcing (W/m2/K) - set to zero
            or False to turn off
        beta: 6-element array of radiative efficiency coefficients in order of
            CH4 concentrations,  W m-2 ppb-1 
            N2O concentrations,  W m-2 ppb-1
            ODS concentrations in EESC,  W m-2 ppt-1
            CO emissions, W m-2 (Mt yr-1)-1
            NMVOC emissions, W m-2 (Mt yr-1)-1
            NOx emissions, W m-2 (MtN yr-1)-1

        emissions_pi: pre-industrial/reference emissions
        concentrations_pi: pre-industrial/reference concentrations

    Outputs:
        ozone ERF time series.
    """

    # we allow 2D output for quick calculation if feedbacks turned off
    if emissions.ndim == 1:
        nspec = len(emissions)
        emissions = emissions.reshape((1, nspec))
    if concentrations.ndim == 1:
        nspec = len(concentrations)
        concentrations = concentrations.reshape((1, nspec))

    nt = emissions.shape[0]

    # calculate EESC for halogens
    cl = np.array(cl_atoms.aslist)
    br = np.array(br_atoms.aslist)
    fc = np.array(fracrel.aslist)

    def eesc(c_ods, c_ods_pi):
        return (
            np.sum(cl * (c_ods-c_ods_pi) * fc/fc[0]) + 
            45 * np.sum(br * (c_ods-c_ods_pi) * fc/fc[0])
        ) * fc[0]


    c_ch4, c_n2o = concentrations[:, [1, 2]].T
#    delta_c_ods = eesc(concentrations[:,15:].T, concentrations_pi[None, 15:])
    c_ods = concentrations[:,15:]
    e_co, e_nmvoc, e_nox = emissions[:,[6, 7, 8]].T
    c_ch4_pi, c_n2o_pi = concentrations_pi[[1, 2]]
    c_ods_pi = concentrations_pi[15:]
    e_co_pi, e_nmvoc_pi, e_nox_pi = emissions_pi[[6, 7, 8]]


    forcing = np.zeros(nt)
    if np.isscalar(temperature):
        temperature = np.ones(nt) * temperature

    for i in range(nt):
        f_ch4   = beta[0] * (c_ch4[i] - c_ch4_pi)
        f_n2o   = beta[1] * (c_n2o[i] - c_n2o_pi)
        f_ods   = beta[2] * eesc(c_ods[i], c_ods_pi)
        f_co    = beta[3] * (e_co[i] - e_co_pi)
        f_nmvoc = beta[4] * (e_nmvoc[i] - e_nmvoc_pi)
        f_nox   = beta[5] * (e_nox[i] - e_nox_pi)
        forcing[i] = f_ch4 + f_n2o + f_ods + f_co + f_nmvoc + f_nox + feedback * temperature[i]

    return forcing
