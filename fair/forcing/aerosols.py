from __future__ import division

import pickle
import os
import sys
import numpy as np
from scipy.interpolate import Rbf
from ..constants import molwt
from ..RCPs.rcp45 import Emissions as r45e


def Stevens(emissions, stevens_params=np.array([0.001875, 0.634, 60.]),
    ref_isSO2=True):
    """Calculates aerosol forcing based on Stevens (2015) that relates sulphate
    aerosol forcing to SOx emissions in a logarithmic fashion.

    Input:
        emissions:   anthropogenic emissions database
    Keywords:
        stevens_params: 3 element array
            0. natural emissions of SOx in Mt/yr
            1. scaling parameter for ERFari (alpha)
            2. scaling parameter for ERFaci (beta)
        ref_isSO2:   True if E_SOx_nat is in units of SO2 rather than S.
    Output:
        F:           aerosol effective radiative forcing
    """

    E_SOx_nat, alpha, beta = stevens_params

    factor = 1
    if ref_isSO2:
        factor = molwt.SO2/molwt.S
    em_SOx = emissions[:,5] * factor

    ERFari = -alpha * em_SOx
    ERFaci = -beta * np.log(em_SOx/E_SOx_nat + 1)
    F = ERFari + ERFaci
    return F


def aerocom_direct(emissions, beta = np.array(
    [-35.29e-4, 0.0, -5.034e-4, -5.763e-4, 453e-4, -37.83e-4, -10.35e-4]),
    scale_AR5=False):

    """Calculates direct aerosol forcing based on linear relationships between
    emissions and forcing in Aerocom models.

    Reference: Myhre et al., 2013: https://www.atmos-chem-phys.net/13/1853/2013

    If inputs from an RCPs SCEN file are used, the units will be correct.

    Inputs: 
        emissions: (nt x 40) emissions array
    Keywords:
        beta: 7-element array of forcing efficiencies in W m-2 (Mt yr-1)-1 for
            SOx, CO, NMVOC, NOx, BC, OC, NH3 (in that order)
        scale_AR5:       If True, scale the forcing output so that the best
                         estimate forcing in 2011 is -0.45 W/m2 based on 2011
                         emissions from the RCPs. The discrepancy between the
                         Aerocom and RCP values comes from two sources; (1)
                         there is no dust forcing considered in Aerocom/FAIR,
                         and (2) the aerosol semi-direct is accounted for by
                         applying this scaling factor.
    Outputs:
        Forcing time series
    """

    em_SOx, em_CO, em_NMVOC, em_NOx, em_BC, em_OC, em_NH3 = \
        emissions[:,[5, 6, 7, 8, 9, 10, 11]].T

    F_SOx    = beta[0] * em_SOx    * molwt.SO2 / molwt.S
    F_CO     = beta[1] * em_CO
    F_NMVOC  = beta[2] * em_NMVOC
    F_NOx    = beta[3] * em_NOx    * molwt.NO / molwt.N
    F_BC     = beta[4] * em_BC
    F_OC     = beta[5] * em_OC
    F_NH3    = beta[6] * em_NH3

    F = F_SOx+F_CO+F_NMVOC+F_NOx+F_BC+F_OC+F_NH3

    if scale_AR5:
        scale=1.3741
    else:
        scale=1.0

    return scale * F


def ghan_indirect(emissions, fix_pre1850_RCP=True, scale_AR5=False,
    ghan_params=np.array([-1.95011431, 0.01107147, 0.01387492])):
    """Estimates the aerosol indirect effect based on the simple model in
    Ghan et al., (2013), doi:10.1002/jgrd.50567.

    This function is just an emulator - a full implementation in Python of the
    Ghan routine (originally coded in Fortran) exists, but will require
    optimisation before it can be used in FAIR. I hope to make the full version
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

    return (F_pd - F_1765) * scale
