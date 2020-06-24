from __future__ import division

import numpy as np

def forcing_to_temperature(t0, q, d, f, e=1.0):
    """Calculate temperature from a given radiative forcing.

    This follows the forcing to temperature function described in Millar et
    al. (2015; 2017).

    Inputs:
        t0: Temperature in timestep t-1
        q: The matrix contributions to slow and fast temperature change
           calculated from ECS and TCR (2 element array)
        d: The slow and fast thermal response time constants (2 element array)
        f: radiative forcing (can be scalar or 1D array representing multiple
           species)

    Keywords:
        e: efficacy factor (default 1); if f is an array, e should be an array
           of the same length.

    Outputs:
        t1: slow and fast contributions to total temperature (2 element array)
        in timestep t
    """
    t1 = t0*np.exp(-1.0/d) + q*(1.0-np.exp((-1.0)/d))*np.sum(f*e)
    return t1


def calculate_q(tcrecs, d, f2x, tcr_dbl, nt):
    """If TCR and ECS are supplied, calculate the q model coefficients.
    See Eqs. (4) and (5) of Millar et al ACP (2017).

    Inputs:
        tcrecs  : 2-element array of transient climate response (TCR) and
                  equilibrium climate sensitivity (ECS).
        d       : The slow and fast thermal response time constants
        f2x     : Effective radiative forcing from a doubling of CO2
        tcr_dbl : time to a doubling of CO2 under 1% per year CO2 increase, yr
        nt      : number of timesteps

    Outputs:
        q       : coefficients of slow and fast temperature change in each
                  timestep ((nt, 2) array).
    """

    # TODO:
    # error checking before call
    # benchmark one call per timestep and if not slower do not convert to 2D
    #  - will make code cleaner

    k = 1.0 - (d/tcr_dbl)*(1.0 - np.exp(-tcr_dbl/d))
    # if ECS and TCR are not time-varying, expand them to 2D array anyway
    if tcrecs.ndim==1:
        if len(tcrecs)!=2:
            raise ValueError(
              "Constant TCR and ECS should be a 2-element array")
        tcrecs = np.ones((nt, 2)) * tcrecs
    elif tcrecs.ndim==2:
        if tcrecs.shape!=(nt, 2):
            raise ValueError(
              "Transient TCR and ECS should be a nt x 2 array")
    q  = (1.0 / f2x) * (1.0/(k[0]-k[1])) * np.array([
        tcrecs[:,0]-tcrecs[:,1]*k[1],tcrecs[:,1]*k[0]-tcrecs[:,0]]).T
    return q

