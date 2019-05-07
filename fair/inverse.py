from __future__ import division

import numpy as np
from scipy.optimize import root
from .forward import forc_to_temp, calculate_q, iirf_simple, iirf_interp
from .forcing.ghg import co2_log
from .defaults import carbon, thermal
from .constants import molwt
from .constants.general import ppm_gtc


def infer_emissions(e1, c1_prescribed, carbon_boxes0, tau_new, a, c_pi):
    """Matches prescribed concentrations to forward-calculated concentrations.

    Inputs:
        e1            : emissions in timestep t, GtC
        c1_prescribed : CO2 concentrations in timestep t, ppmv
        carbon_boxes0 : carbon stored in each atmospheric reservoir at timestep
                        t-1 (GtC)
        tau_new       : decay time constants of CO2 (yr)
        a             : partition coefficient of carbon boxes
        c_pi          : pre-industrial concentration of CO2, ppmv
    """
    
    c1_calculated = np.sum(carbon_boxes0*np.exp(-1.0/tau_new) + a*e1 / ppm_gtc) + c_pi
    return c1_calculated-c1_prescribed


def inverse_carbon_cycle(c1, c_acc0, temp, r0, rc, rt, iirf_max, time_scale_sf,
                         a, tau, iirf_h, carbon_boxes0, c_pi, c0, e0):
    """Calculates CO2 emissions from concentrations.
    
    Inputs:
        c1            : concentration of CO2 in timestep t, ppmv
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
        e0            : emissions of CO2 in timestep t, GtC

    Outputs:
        e1            : emissions of CO2 in timestep t, GtC
        c_acc1        : cumulative airborne carbon anomaly (GtC) since
                        pre-industrial, timestep t
        carbon_boxes1 : carbon stored in each atmospheric reservoir at timestep
                        t (GtC)
        time_scale_sf : scale factor for CO2 decay constants
    """
     
    iirf = iirf_simple(c_acc0, temp, r0, rc, rt, iirf_max)
    time_scale_sf = root(iirf_interp, time_scale_sf,
      args=(a, tau, iirf_h, iirf))['x']
    tau_new = tau * time_scale_sf
    e1 = root(infer_emissions, e0, args=(c1, carbon_boxes0, tau_new, a, c_pi))['x']
    c_acc1 = c_acc0 + 0.5*(e1 + e0) - (c1 - c0)*ppm_gtc
    carbon_boxes1 = carbon_boxes0*np.exp(-1.0/tau_new) + a*e1 / ppm_gtc
    return e1, c_acc1, carbon_boxes1, time_scale_sf

    
def inverse_fair_scm(
# Longer term goal: one calling interface only - appropriate function
# determined from call
    C             = None,
    other_rf      = 0.0,
    q             = thermal.q,
    tcrecs        = thermal.tcrecs,
    d             = thermal.d,
    F2x           = thermal.f2x,
    tcr_dbl       = thermal.tcr_dbl,
    a             = carbon.a,
    tau           = carbon.tau,
    r0            = carbon.r0,
    rc            = carbon.rc,
    rt            = carbon.rt,
    iirf_max      = carbon.iirf_max,
    iirf_h        = carbon.iirf_h,
    C_pi          = 278.,
    time_scale_sf = 0.16
    ):
    
    # Error checking and validation goes here...

    # Dimensions
    nt = len(C)
    carbon_boxes_shape = (nt, a.shape[0])
    thermal_boxes_shape = (nt, d.shape[0])
    
    # Thermal response
    q         = calculate_q(tcrecs, d, F2x, tcr_dbl, nt)
    
    # Allocate intermediate and output arrays
    C_acc     = np.zeros(nt)
    R_i       = np.zeros(carbon_boxes_shape)
    emissions = np.zeros(nt)
    T_j       = np.zeros(thermal_boxes_shape)
    F         = np.zeros(nt)
    
    if np.isscalar(other_rf):
        other_rf = other_rf * np.ones(nt)
    
    # First timestep
    emissions[0] = root(infer_emissions, 0., args=(C[0], R_i[0],
        tau, a, C_pi))['x']
    F[0]         = co2_log(C[0], C_pi, F2x=F2x) + other_rf[0]
    T_j[0,:]     = forc_to_temp(T_j[0,:], q[0,:], d, F[0])

    # Second timestep onwards
    for t in range(1,nt):
        emissions[t], C_acc[t], R_i[t,:], time_scale_sf = (
            inverse_carbon_cycle(
                C[t], C_acc[t-1], np.sum(T_j[t-1,:]), r0, rc, rt, iirf_max, 
                time_scale_sf, a, tau, iirf_h, R_i[t-1,:],
                C_pi, C[t-1], emissions[t-1]
            )
        )
        F[t]     = co2_log(C[t], C_pi, F2x=F2x) + other_rf[t]
        T_j[t,:] = forc_to_temp(T_j[t-1,:], q[t,:], d, F[t])

    # Output temperatures
    T = np.sum(T_j, axis=-1)
    
    return emissions, F, T
