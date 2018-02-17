from __future__ import division

import numpy as np
from ..constants import molwt

def regress(em_CH4, em_CO, em_NMVOC, em_NOx,
            beta=np.array([2.8249e-4, 1.0695e-4, -9.3604e-4, 99.7831e-4])):

    """Calculates tropospheric ozone forcing from precursor emissions.

    Inputs: emissions databases; can deal with 2D (year, species) or 1D (year)
        em_CH4
        em_CO
        em_NMVOC
        em_NOx

    Keywords:
        beta: 4-element array of regression coefficients of precursor
              radiative efficiency, W m-2 (Mt yr-1)-1.
              order is [CH4, CO, NMVOC, NOx]

    Outputs:
        tropospheric ozone ERF time series.
    """

    em       = np.array([em_CH4, em_CO, em_NMVOC, em_NOx]).T

    F = np.sum(beta * em, axis=-1)
    return F


def stevenson(year, 
              C_CH4, 
              em_CO, 
              em_NMVOC,
              em_NOx,
              T=0,
              feedback=False,
              fix_pre1850_RCP=False,
             ):
    """Calculates tropospheric ozone forcing from precursor emissions based on
    Stevenson et al, 2013 10.5194/acp-13-3063-2013

    Inputs (all 1D numpy arrays)
        year     
        C_CH4    : methane concentration (ppb)
        em_CO    : CO emissions (Mt/yr)
        em_NMVOC : NMVOC emissions (Mt/yr)
        em_NOx   : NOx emissions (Mt/yr)

    Keywords:
        T              : change in surface temperature since pre-industrial
        feedback       : True or False - include temperature feedback on ozone
                         forcing?
        fix_pre1850_RCP: Use different relationship for 1750/65 to 1850 based 
                         on anthropogenic emissions from Skeie et al (2011)
                         for 1750 (atmos-chem-phys.net/11/11827/2011)

    Outputs:
        tropospheric ozone ERF time series.
    """

    # numbers in denominator are 2000-1750 concs or emissions used in 
    # Stevenson and traced back to Lamarque et al 2010 for 2000
    # https://www.atmos-chem-phys.net/10/7017/2010/
    if year>=1850 or fix_pre1850_RCP==False:
        F_CH4   = 0.166/960    * (C_CH4-722)
        F_CO    = 0.058/893.39 * (em_CO-170)
        F_NMVOC = 0.035/155.84 * (em_NMVOC-10)
        F_NOx   = 0.119/61.16  * (em_NOx * molwt.NO / molwt.N - 4.29)
    # The RCP scenarios give a negative forcing prior to ~1780. This is 
    # because the anthropogenic emissions are given to be zero in RCPs but
    # not zero in the Skeie numbers which are used here. This can be fixed
    # to give a more linear behaviour.
    else:
        F_CH4   = 0.166/960    * (C_CH4-722)
        F_CO    = 0.058/893.39 * 215.59  * em_CO / 385.59
        F_NMVOC = 0.035/155.84 * 51.97 * em_NMVOC / 61.97
        F_NOx   = 0.119/61.16  * 7.31 * (em_NOx * molwt.NO / molwt.N) / 11.6

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
