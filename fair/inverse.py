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
    time_scale_sf = 0.16,
    restart_in    = False,
    restart_out   = False,
    ):

    """Diagnoses emissions from prescribed concentrations.

    Inputs:
        C             : concentrations of CO2, ppmv
        other_rf      : non-CO2 radiative forcing (scalar or numpy array, W/m2)
        q             : coefficients of slow and fast temperature change.
                        Overridden if tcrecs is specified.
        tcrecs        : transient climate response and equilibrium climate
                        sensitivity array (2-element or (nt, 2))
        d             : timescales of slow and fast contribution to temperature
                        change
        F2x           : radiative forcing from a doubling of CO2 concentrations
                        (W/m2)
        tcr_dbl       : timescale over which a 1% compound increase of CO2 acts
                        (yr)
        a             : partition fractions for CO2 boxes
        tau           : time constants for CO2 boxes
        r0            : pre-industrial time-integrated airborne fraction (yr)
        rc            : sensitivity of time-integrated airborne fraction to
                        airborne carbon (yr/GtC)
        rt            : sensitivity of time-integrated airborne fraction to
                        temperature (yr/K)
        iirf_max      : maximum value of time-integrated airborne fraction (yr)
        iirf_h        : time horizon for time-integrated airborne fraction
        C_pi          : pre-industrial concentration of CO2, ppmv
        time_scale_sf : initial guess for scaling factor of CO2 time constants.
                        Overridden if using a restart.
        restart_in    : Allows a restart of the carbon cycle from a non-initial
                        state. A 6-tuple of:
                          array of accumulated carbon in each atmospheric box,
                          array of slow and fast temperature contributions,
                          total accumulated carbon,
                          emissions in the timestep before restart
                          time constant scale factor
                          CO2 concentrations in the timestep before restart
        restart_out   : if True, return the restart state as an extra output.
                        See restart_in.
    Outputs:
        E             : Timeseries of diagnosed CO2 emissions in GtC
        F             : Timeseries of total radiative forcing, W/m2
        T             : Timeseries of temperature anomaly since pre-industrial
        restart       : if restart_out=True, 6-tuple of carbon cycle state
                        parameters. See restart_in.
    """
    
    # Error checking and validation goes here...

    # Dimensions
    nt = len(C)
    carbon_boxes_shape = (nt, a.shape[0])
    thermal_boxes_shape = (nt, d.shape[0])
    
    # Thermal response
    if type(tcrecs) is np.ndarray:
        q     = calculate_q(tcrecs, d, F2x, tcr_dbl, nt)
    
    # Allocate intermediate and output arrays
    C_acc     = np.zeros(nt)
    R_i       = np.zeros(carbon_boxes_shape)
    emissions = np.zeros(nt)
    T_j       = np.zeros(thermal_boxes_shape)
    F         = np.zeros(nt)
    
    if np.isscalar(other_rf):
        other_rf = other_rf * np.ones(nt)
    
    # First timestep
    if restart_in:
        R_i_minus1    = restart_in[0]
        T_j_minus1    = restart_in[1]
        C_acc_minus1  = restart_in[2]
        E_minus1      = restart_in[3]
        time_scale_sf = restart_in[4]
        C_minus1      = restart_in[5]
        emissions[0], C_acc[0], R_i[0,:], time_scale_sf = (
            inverse_carbon_cycle(
                C[0], C_acc_minus1, np.sum(T_j_minus1), r0, rc, rt, iirf_max,
                time_scale_sf, a, tau, iirf_h, R_i_minus1,
                C_pi, C_minus1, E_minus1
            )
        )
        F[0]          = co2_log(C[0], C_pi, F2x=F2x) + other_rf[0]
        T_j[0,:]      = forc_to_temp(T_j_minus1, q[0,:], d, F[0])
    else:
        emissions[0]  = root(infer_emissions, 0., args=(C[0], R_i[0,:],
            tau, a, C_pi))['x']
        F[0]          = co2_log(C[0], C_pi, F2x=F2x) + other_rf[0]
        T_j[0,:]      = forc_to_temp(T_j[0,:], q[0,:], d, F[0])

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
    if restart_out:
        restart_out_val = (R_i[-1], T_j[-1], C_acc[-1], emissions[-1],
            time_scale_sf, C[-1])
        return emissions, F, T, restart_out_val
    else:
        return emissions, F, T
