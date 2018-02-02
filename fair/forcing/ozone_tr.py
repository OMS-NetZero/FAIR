import numpy as np

def regress(emissions,
            beta=np.array([2.8249e-4, 1.0695e-4, -9.3604e-4, 99.7831e-4])):

    """Calculates tropospheric ozone forcing from precursor emissions.

    Inputs:
        emissions: emissions database
    Keywords:
        beta: 4-element array of regression coefficients of precursor
              radiative efficiency, W m-2 (Mt yr-1)-1.
              order is [CH4, CO, NMVOC, NOx]
    Outputs:
        tropospheric ozone ERF time series.
    """

    em_CH4   = emissions[:,3]
    em_CO    = emissions[:,6]
    em_NMVOC = emissions[:,7]
    em_NOx   = emissions[:,8]

    em       = np.array([em_CH4, em_CO, em_NMVOC, em_NOx]).T

    F = np.sum(beta * em, axis=-1)
    return F

