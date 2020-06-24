from __future__ import division

import numpy as np
from scipy.optimize import root

from ..constants.general import ppm_gtc

"""Carbon cycle function from FaIR v1.0.0."""

def carbon_cycle(e0, c_acc0, temp, r0, rc, rt, iirf_max, time_scale_sf0, a, tau,
    iirf_h, carbon_boxes0, c_pi, c0, e1):
    """Calculates CO2 concentrations from emissions.

    Inputs:
        e0            : emissions of CO2 (GtC) in timestep t-1
        c_acc0        : cumulative airborne carbon anomaly (GtC) since
                        pre-industrial, timestep t-1
        temp          : temperature anomaly above pre-industrial (K)
        r0            : pre-industrial time-integrated airborne fraction (yr)
        rc            : sensitivity of time-integrated airborne fraction to
                        airborne carbon (yr/GtC)
        rt            : sensitivity of time-integrated airborne fraction to
                        temperature (yr/K)
        iirf_max      : maximum value of time-integrated airborne fraction (yr)
        time_scale_sf0: initial guess of alpha scaling factor
        a             : partition coefficient of carbon boxes
        tau           : present-day decay time constants of CO2 (yr)
        iirf_h        : time horizon for time-integrated airborne fraction (yr)
        carbon_boxes0 : carbon stored in each atmospheric reservoir at timestep
                        t-1 (GtC)
        c_pi          : pre-industrial concentration of CO2, ppmv
        c0            : concentration of CO2 in timestep t-1, ppmv
        e1            : emissions of CO2 in timestep t, GtC

    Outputs:
        c1            : concentrations of CO2 in timestep t, ppmv
        c_acc1        : cumulative airborne carbon anomaly (GtC) since
                        pre-industrial, timestep t
        carbon_boxes1 : carbon stored in each atmospheric reservoir at timestep
                        t (GtC)
        time_scale_sf : scale factor for CO2 decay constants
    """
    iirf = _iirf_simple(c_acc0, temp, r0, rc, rt, iirf_max)
    time_scale_sf = root(_iirf_interp, time_scale_sf0,
      args=(a, tau, iirf_h, iirf))['x']
    tau_new = tau * time_scale_sf
    carbon_boxes1 = carbon_boxes0*np.exp(-1.0/tau_new) + a*e1 / ppm_gtc
    c1 = np.sum(carbon_boxes1) + c_pi
    c_acc1 = c_acc0 + 0.5*(e1 + e0) - (c1 - c0)*ppm_gtc
    return c1, c_acc1, carbon_boxes1, time_scale_sf


def _iirf_interp(alp_b,a,tau,iirf_h,targ_iirf):
    """Interpolation function for finding alpha, the CO2 decay time constant
    scaling factor, in iirf_h equation. See Eq. (7) of Millar et al ACP (2017).

    Inputs:
        alp_b    : Guess for alpha, the scale factor, for tau
        a        : partition fractions for CO2 boxes
        tau      : time constants for CO2 boxes
        iirf_h   : time horizon for time-integrated airborne fraction
        targ_iirf: iirf_h calculated using simple parameterisation (Eq. (8),
                   Millar et al (2017)).
    """

    iirf_arr = alp_b*(np.sum(a*tau*(1.0 - np.exp(-iirf_h/(tau*alp_b)))))
    return iirf_arr - targ_iirf


def _iirf_simple(c_acc, temp, r0, rc, rt, iirf_max):
    """Simple linear iIRF relationship. Eq. (8) of Millar et al ACP (2017).

    Inputs:
        c_acc    : cumulative airborne carbon anomaly (GtC) since
                   pre-industrial
        temp     : temperature anomaly since pre-industrial
        r0       : pre-industrial time-integrated airborne fraction (yr)
        rc       : sensitivity of time-integrated airborne fraction to airborne
                   carbon (yr/GtC)
        rt       : sensitivity of time-integrated airborne fraction to
                   temperature (yr/K)
        iirf_max : maximum value of time-integrated airborne fraction (yr)

    Outputs:
        iirf     : time-integrated airborne fraction of carbon (yr)
    """

    return np.min([r0 + rc * c_acc + rt * temp, iirf_max])
