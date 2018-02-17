import numpy as np

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


def stevenson(C_CH4, em_CO, em_NMVOC, em_NOx, C_CH4pi=722.0, T=0, feedback=-0.02):
    """Calculates tropospheric ozone forcing from precursor emissions based on
    Stevenson et al, 2013 10.5194/acp-13-3063-2013

    Inputs:
        C_CH4    : methane concentration (ppb)
        em_CO    : CO emissions (Mt/yr)
        em_NMVOC : NMVOC emissions (Mt/yr)
        em_NOx   : NOx emissions (Mt/yr)

    Keywords:
        C_CH4pi  : pre-industrial CH4 concentration (ppb)
        T        : change in surface temperature since pre-industrial
        feedback : change in tropospheric ozone forcing per Kelvin temperature
                   increase (table 12 Stevenson, value for year 2000 used, and
                   assumed 1850-2000 temperature change is about 0.83K)

    Outputs:
        tropospheric ozone ERF time series.
    """

    F_CH4   = 0.166/960 * (C_CH4-C_CH4pi)
    F_CO    = 0.058/681.8 * em_CO
    F_NMVOC = 0.035/145.84 * em_NMVOC
    F_NOx   = 0.119/56.87 * em_NOx

    F = F_CH4 + F_CO + F_NMVOC + F_NOx + feedback*T
    return F
