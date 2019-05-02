from __future__ import division

import inspect
import numpy as np
import warnings
from scipy.optimize import root
from .ancil import natural, cmip6_volcanic, cmip6_solar, historical_scaling
from .constants import molwt, lifetime, radeff
from .constants.general import M_ATMOS
from .forcing import ozone_tr, ozone_st, h2o_st, contrails, aerosols, bc_snow,\
                                         landuse
from .forcing.ghg import co2_log


def carbon_cycle(e0, c_acc0, temp, r0, rc, rt, iirf_max, time_scale_sf0, a, tau,
    iirf_h, carbon_boxes0, ppm_gtc, c_pi, c0, e1):
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
        ppm_gtc       : conversion fraction from ppmv to GtC
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
    iirf = iirf_simple(c_acc0, temp, r0, rc, rt, iirf_max)
    time_scale_sf = root(iirf_interp, time_scale_sf0,
      args=(a, tau, iirf_h, iirf))['x']
    tau_new = tau * time_scale_sf
    carbon_boxes1 = carbon_boxes0*np.exp(-1.0/tau_new) + a*e1 / ppm_gtc
    c1 = np.sum(carbon_boxes1) + c_pi
    c_acc1 = c_acc0 + 0.5*(e1 + e0) - (c1 - c0)*ppm_gtc
    return c1, c_acc1, carbon_boxes1, time_scale_sf
