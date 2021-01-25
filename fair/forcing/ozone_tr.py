from __future__ import division

import numpy as np
from ..constants import molwt

def regress(emissions,
            beta=np.array([2.8249e-4, 1.0695e-4, -9.3604e-4, 99.7831e-4])):

    """Calculates tropospheric ozone forcing from precursor emissions.

    Inputs: (nt x 40) emissions array

    Keywords:
        beta: 4-element array of regression coefficients of precursor
              radiative efficiency, W m-2 (Mt yr-1)-1.
              order is [CH4, CO, NMVOC, NOx]

    Outputs:
        tropospheric ozone ERF time series.
    """

    if emissions.ndim==2:
        em_CH4, em_CO, em_NMVOC, em_NOx = emissions[:,[3, 6, 7, 8]].T
    else:
        em_CH4, em_CO, em_NMVOC, em_NOx = emissions[[3, 6, 7, 8]]

    F_CH4   = beta[0] * em_CH4
    F_CO    = beta[1] * em_CO
    F_NMVOC = beta[2] * em_NMVOC
    F_NOx   = beta[3] * em_NOx

    F = F_CH4 + F_CO + F_NMVOC + F_NOx
    return F


def cmip6_stevenson(emissions, C_CH4, T=0, feedback=False,
    PI=np.array([722, 170, 10, 4.29]), 
    beta=np.array([1.77871043e-04, 5.80173377e-05, 2.09151270e-03,
        1.94458719e-04])):
    """Calculates tropospheric ozone forcing from precursor emissions based on
    Stevenson et al, 2013 10.5194/acp-13-3063-2013

    Inputs:
        emissions: (nt x 40) numpy array     
        C_CH4    : (nt) numpy array of methane concentrations, ppb

    Keywords:
        T              : change in surface temperature since pre-industrial
        feedback       : True or False - include temperature feedback on ozone
                         forcing?
        PI:            : 4-element array of pre-industrial CH4 concentrations,
                         CO emissions, NMVOC emissions and NOx emissions
        beta:          : coefficients of how CH4 concentrations, CO emissions,
                         NMVOC emissions and NOx emissions affect forcing

    Outputs:
        tropospheric ozone ERF time series.
    """

    # expand to 2D/1D if not already
    if emissions.ndim == 1:
        nspec = len(emissions)
        emissions = emissions.reshape((1, nspec))
    if np.isscalar(C_CH4):
        C_CH4 = np.ones(1)*C_CH4

    year, em_CO, em_NMVOC, em_NOx = emissions[:,[0, 6, 7, 8]].T
    nt = len(year)
    F_CH4, F_CO, F_NMVOC, F_NOx = np.zeros((4,nt))

    for i in range(nt):
        F_CH4[i]   = beta[0] * (C_CH4[i]-PI[0])
        F_CO[i]    = beta[1] * (em_CO[i]-PI[1])
        F_NMVOC[i] = beta[2] * (em_NMVOC[i]-PI[2])
        F_NOx[i]   = beta[3] * (em_NOx[i]-PI[3])

    # Include the effect of climate feedback? We fit a curve to the 2000, 2030
    # and 2100 best estimates of feedback based on middle-of-the-road
    # temperature projections.
    def temperature_feedback(T, a=0.03189267, b=1.34966941, c=-0.03214807):
        if T<=0:
            return 0
        else:
            return a*np.exp(-b*T)+c

    if feedback:
        F = F_CH4 + F_CO + F_NMVOC + F_NOx + temperature_feedback(T)
    else:
        F = F_CH4 + F_CO + F_NMVOC + F_NOx
    return F


def stevenson(emissions, C_CH4, T=0, feedback=False, fix_pre1850_RCP=False,
    PI=np.array([722, 170, 10, 4.29])):
    """Calculates tropospheric ozone forcing from precursor emissions based on
    Stevenson et al, 2013 10.5194/acp-13-3063-2013

    Inputs:
        emissions: (nt x 40) numpy array     
        C_CH4    : (nt) numpy array of methane concentrations, ppb

    Keywords:
        T              : change in surface temperature since pre-industrial
        feedback       : True or False - include temperature feedback on ozone
                         forcing?
        fix_pre1850_RCP: Use different relationship for 1750/65 to 1850 based 
                         on anthropogenic emissions from Skeie et al (2011)
                         for 1750 (atmos-chem-phys.net/11/11827/2011)
        PI:            : 4-element array of pre-industrial CH4 concentrations,
                         CO emissions, NMVOC emissions and NOx emissions

    Outputs:
        tropospheric ozone ERF time series.
    """

    # expand to 2D/1D if not already
    if emissions.ndim == 1:
        nspec = len(emissions)
        emissions = emissions.reshape((1, nspec))
    if np.isscalar(C_CH4):
        C_CH4 = np.ones(1)*C_CH4

    # numbers in denominator are 2000-1750 concs or emissions used in 
    # Stevenson and traced back to Lamarque et al 2010 for 2000
    # https://www.atmos-chem-phys.net/10/7017/2010/
    year, em_CO, em_NMVOC, em_NOx = emissions[:,[0, 6, 7, 8]].T
    nt = len(year)
    F_CH4, F_CO, F_NMVOC, F_NOx = np.zeros((4,nt))

    for i in range(nt):
        if year[i]>=1850 or fix_pre1850_RCP==False:
            F_CH4[i]   = 0.166/960    * (C_CH4[i]-PI[0])
            F_CO[i]    = 0.058/681.8 * (em_CO[i]-PI[1])
            F_NMVOC[i] = 0.035/155.84 * (em_NMVOC[i]-PI[2])
            F_NOx[i]   = 0.119/61.16  * (em_NOx[i] * 
                molwt.NO / molwt.N - PI[3])
        # The RCP scenarios give a negative forcing prior to ~1780. This is 
        # because the anthropogenic emissions are given to be zero in RCPs but
        # not zero in the Skeie numbers which are used here. This can be fixed
        # to give a more linear behaviour.
        else:
            F_CH4[i]   = 0.166/960    * (C_CH4[i]-722)
            F_CO[i]    = 0.058/681.8 * 215.59  * em_CO[i] / 385.59
            F_NMVOC[i] = 0.035/155.84 * 51.97 * em_NMVOC[i] / 61.97
            F_NOx[i]   = 0.119/61.16  * 7.31 * (em_NOx[i]
                 * molwt.NO / molwt.N) / 11.6

    # Include the effect of climate feedback? We fit a curve to the 2000, 2030
    # and 2100 best estimates of feedback based on middle-of-the-road
    # temperature projections.
    def temperature_feedback(T, a=0.03189267, b=1.34966941, c=-0.03214807):
        if T<=0:
            return 0
        else:
            return a*np.exp(-b*T)+c

    if feedback:
        F = F_CH4 + F_CO + F_NMVOC + F_NOx + temperature_feedback(T)
    else:
        F = F_CH4 + F_CO + F_NMVOC + F_NOx

    return F
