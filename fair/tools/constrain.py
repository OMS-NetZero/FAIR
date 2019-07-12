from __future__ import division

import numpy as np
from scipy import stats

def hist_temp(Tobs, Tmodel, years, inflate=True, CI=0.9):
    """Checks to see whether model-derived temperatures fall in observational
    uncertainty.

    Uses the method of Thompson et al., 2015, also used in IPCC AR5 to derive
    temperature trends which includes autocorrelation. The regression slope of
    the observations is compared to the regression slope of the model. If the
    trend of the model is within observational uncertainty, the test passes.

    Reference: J. Climate, 28, 6443-6456 10.1175/JCLI-D-14-00830.1

    inputs:
        Tobs: observed temperature time series. Numpy array
        Tmodel: modelled temperature time series. Numpy array
        years: Numpy array of years covered by Tobs and Tmodel

    keywords:
        inflate: True (default) if the uncertainty bounds should be inflated
            for lag-1 autocorrelation - as used in Thompson.
        CI: confidence interval around the mean regression slope to count as
            constrained. Default 0.9.

    returns:
        accept: True if ensemble member agrees with observations else False.
        slope_m: regression slope of modelled temperature
        intercept_m: intercept of modelled temperature
        slope_o: regression slope of observed temperature
        intercept_o: intercept of observed temperature
    """

    n = float(len(years))

    # detrend the 1880-2016 observations and apply the internal variability
    # estimate as detailed in eq. 8 of Thompson et al, J. Climate 2015

    slope_o, intercept_o, _, _, stderr_o = stats.linregress(years,Tobs)
    resid = slope_o * years + intercept_o - Tobs
    s = np.std(resid, ddof=2)
    if inflate:
        lag1ac = np.corrcoef(resid[:-1],resid[1:])[1,0]
    else:
        lag1ac = 0
    gamma = np.sqrt( n / ((n-2) * (1-lag1ac)/(1+lag1ac) - 2 ))
    g = np.sqrt(12./(n**3 - n))
    tcrit = stats.t.ppf(1-(1-CI)/2.0, df=n)
    CI_o = s * gamma * tcrit * g

    # check modelled T against observed
    slope_m,intercept_m,_,_,stderr_m = stats.linregress(years, Tmodel)

    if slope_o-CI_o <= slope_m <= slope_o+CI_o:
        accept = True
    else:
        accept = False
    return accept, slope_m, intercept_m, slope_o, intercept_o
