from __future__ import division

import numpy as np
from ..constants import molwt
from ..RCPs.rcp45 import Emissions as r45e


def Stevens(emissions, stevens_params=np.array([0.001875, 0.634, 60.]),
    ref_isSO2=True, E_pi=0):
    """Calculates aerosol forcing based on Stevens (2015) that relates sulphate
    aerosol forcing to SOx emissions in a logarithmic fashion.

    Input:
        emissions:   anthropogenic emissions database
    Keywords:
        stevens_params: 3 element array
            0. scaling parameter for ERFari (alpha)
            1. scaling parameter for ERFaci (beta)
            2. natural emissions of SOx in Mt/yr
        ref_isSO2:   True if E_SOx_nat is in units of SO2 rather than S.
        E_pi: pre-industrial/reference emissions of SO2 (or S).
    Output:
        ERFari, ERFaci:  aerosol effective radiative forcing due to 
            aerosol-radiation interactions and aerosol-cloud interactions.
    """

    alpha, beta, E_SOx_nat = stevens_params

    factor = 1
    if ref_isSO2:
        factor = molwt.SO2/molwt.S
    em_SOx = emissions[:,5] * factor
    em_pi  = E_pi * factor

    ERFari = -alpha * (em_SOx-em_pi)
#    ERFaci = (
#        (-beta * np.log(em_SOx/E_SOx_nat + 1)) - 
#        (-beta * np.log(em_pi/E_SOx_nat + 1)) )
    ERFaci = (-beta * np.log((em_SOx-em_pi)/E_SOx_nat + 1))
    return ERFari, ERFaci


def aerocom_direct(emissions,
        beta = np.array(
          [-6.2227e-3, 0.0, -3.8392e-4, -1.16551e-3, 1.601537e-2, -1.45339e-3,
           -1.55605e-3]), E_pi=np.zeros(40), diagnostics=None
    ):

    """Calculates direct aerosol forcing based on linear relationships between
    emissions and forcing in Aerocom models.

    Reference: Myhre et al., 2013: https://www.atmos-chem-phys.net/13/1853/2013

    If inputs from an RCPs SCEN file are used, the units will be correct.

    Inputs:
        emissions: (nt x 40) emissions array
    Keywords:
        beta: 7-element array of forcing efficiencies in W m-2 (Mt yr-1)-1 for
            SOx, CO, NMVOC, NOx, BC, OC, NH3 (in that order)
        E_pi: pre-industrial emissions (40 element array)
        diagnostics: if 'AR6', give split of direct aerosol effect by species
    Outputs:
        Forcing time series
    """

    em_SOx, em_CO, em_NMVOC, em_NOx, em_BC, em_OC, em_NH3 = \
        emissions[:,[5, 6, 7, 8, 9, 10, 11]].T

    F_SOx    = beta[0] * (em_SOx   - E_pi[5])
    F_CO     = beta[1] * (em_CO    - E_pi[6])
    F_NMVOC  = beta[2] * (em_NMVOC - E_pi[7])
    F_NOx    = beta[3] * (em_NOx   - E_pi[8])
    F_BC     = beta[4] * (em_BC    - E_pi[9])
    F_OC     = beta[5] * (em_OC    - E_pi[10])
    F_NH3    = beta[6] * (em_NH3   - E_pi[11])

    if diagnostics == 'AR6':
        ERFari = np.column_stack([F_SOx, F_CO, F_NMVOC, F_NOx, F_BC, F_OC, F_NH3])
    else:
        ERFari = F_SOx+F_CO+F_NMVOC+F_NOx+F_BC+F_OC+F_NH3
    
    return ERFari


def ghan_indirect(emissions, fix_pre1850_RCP=True, scale_AR5=False,
    ghan_params=np.array([-1.95011431, 0.01107147, 0.01387492]),
    E_pi=np.zeros(40)):
    """Estimates the aerosol indirect effect based on the simple model in
    Ghan et al., (2013), doi:10.1002/jgrd.50567.

    This function is just an emulator - a full implementation in Python of the
    Ghan routine (originally coded in Fortran) exists, but will require
    optimisation before it can be used in FaIR. I hope to make the full version
    available in a future version.

    A 500-member Latin Hypercube sample of emissions of SOx, NMVOC, BC and OC
    was prepared offline and run through the Ghan simple model and a functional
    relationship fit to the output. SOA aerosol (based on NMVOC emissions) is
    sometimes unreliable and does not exert a strong dependence on the ERF, and
    OC+BC is parameterised as primary organic matter, so the resulting output
    is a function of SOx and (BC+OC) emissions.

    Inputs:
        emissions: (nt x 40) numpy emissions array
    Keywords:
        fix_pre1850_RCP: Use different relationship for 1750/65 to 1850 based
                         on anthropogenic emissions from Skeie et al (2011)
                         for 1750 (atmos-chem-phys.net/11/11827/2011)
        scale_AR5:       If True, scale the forcing output so that the best
                         estimate forcing in 2011 is -0.45 W/m2 based on 2011
                         emissions from the RCPs. The Ghan emulator is built on
                         results from the CAM5 GCM. As reported in AR5 WG1 Ch7,
                         GCMs tend to overestimate forcing from aerosol-cloud
                         interactions.
        ghan_params:     3-element numpy array
                         0: scale factor
                         1: sensitivity to SOx emissions
                         2: sensitivity to BC+OC emissions
    Outputs:
        Forcing timeseries
    """

    year, em_SOx, em_BC, em_OC = emissions[:,[0, 5, 9, 10]].T

    def _ERFaci(em,
        ghan_params=np.array([-1.95011431, 0.01107147, 0.01387492])):
        scale = ghan_params[0]
        b_SOx = ghan_params[1]
        b_POM = ghan_params[2]
        return scale*np.log(1+b_SOx*em[0]+b_POM*em[1])

    # PI forcing was not zero as there were some emissions. Use estimates
    # from Skeie et al, 2011 for 1750 forcing.
    E_1765 = np.array([1.0, 11.2])
    nt = len(year)
    F_pd = np.zeros(nt)
    for i in range(nt):
        if year[i]>=1850 or fix_pre1850_RCP==False:
            F_pd[i] = _ERFaci([em_SOx[i], em_BC[i]+em_OC[i]],
                              ghan_params=ghan_params)
        else:
            # linearly interpolate between 1765 and 1850
            E_1850 = np.array([r45e.sox[85], r45e.bc[85]+r45e.oc[85]])
            F_pd[i] = _ERFaci((year[i]-1765)/85.*E_1850 +
                              (1850-year[i])/85.*E_1765,
                              ghan_params=ghan_params)

    # 1765 emissions = zero forcing
    F_1765 = -0.3002836449793625
    F_2011 = -1.5236182344467388

    # are we rescaling to AR5 best estimate with the default parameters?
    if scale_AR5:
        scale=-0.45/(F_2011-F_1765)
    else:
        scale=1.0

    ERFaci = (F_pd - F_1765) * scale
    return ERFaci


def ghan2(emissions, E_pi, ghan_params):
    """temphack for fair1.6"""
    beta, n_so2, n_pom = ghan_params
    pd_re = -beta * np.log(1 + emissions[:,5]/n_so2 + emissions[:,9:11].sum(axis=1)/n_pom)
    pi_re = -beta * np.log(1 + E_pi[5]/n_so2 + E_pi[9:11].sum()/n_pom)
    return pd_re - pi_re
