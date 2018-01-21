from __future__ import division

import numpy as np
from ..constants import molwt

def regress(emissions,
           aNOx =  -2.62e-3,
           aSOx =  -5.66e-3,
           aBC  =  14.28e-3,
           aOC  =  -8.90e-3,
           aNH3 =  -5.89e-3):

    """Calculates aerosol forcing based on a multi-linear regression of
    aerosol emissions to the aerosol ERF time series in AR5.

    Input:
        emissions:   anthropogenic emissions database
    Keywords:
        aSOx, aNOx, aBC, aOC, aNH3:  radiative efficiency of each aerosol
                     precursor, W/m2/(Mt/yr)
    Output:
        F:           aerosol effective radiative forcing
    """

    em_SOx   = emissions[:,5]
    em_NOx   = emissions[:,8]
    em_BC    = emissions[:,9]
    em_OC    = emissions[:,10]
    em_NH3   = emissions[:,11]

    em       = np.array([em_SOx, em_NOx, em_BC, em_OC, em_NH3]).T
    beta     = np.array([aSOx, aNOx, aBC, aOC, aNH3])

    F = np.sum(beta * em, axis=-1)
    return F


def Stevens(emissions, E_SOx_nat=60, alpha=0.001875, beta=0.634, ref_isSO2=True):
    """Calculates aerosol forcing based on Stevens (2015) that relates sulphate
    aerosol forcing to SOx emissions in a logarithmic fashion.

    Input:
        emissions:   anthropogenic emissions database
    Keywords:
        E_SOx_nat:   natural emissions of SOx in Mt/yr
        alpha:       scaling parameter for ERFari
        beta:        scaling parameter for ERFaci
        ref_isSO2:   True if E_SOx_nat is in units of SO2 rather than S.
    Output:
        F:           aerosol effective radiative forcing
    """

    factor = 1
    if ref_isSO2:
        factor = molwt.SO2/molwt.S
    em_SOx = emissions[:,5] * factor

    ERFari = -alpha * em_SOx
    ERFaci = -beta * np.log(em_SOx/E_SOx_nat + 1)
    F = ERFari + ERFaci
    return F

