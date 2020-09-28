from __future__ import division

import numpy as np
from scipy.optimize import root
from .gas_cycle.fair1 import _iirf_simple, _iirf_interp
from .forcing.ghg import co2_log
from .defaults import carbon, thermal
from .constants import molwt
from .constants.general import ppm_gtc
from .temperature.millar import forcing_to_temperature, calculate_q


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
     
    iirf = _iirf_simple(c_acc0, temp, r0, rc, rt, iirf_max)
    time_scale_sf = root(_iirf_interp, time_scale_sf,
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
    F_in          = None,
    temperature_function = 'Millar',
    lambda_global=1.18,  # this and the below only used in two-layer model
    ocean_heat_capacity=np.array([8.2, 109.0]),
    ocean_heat_exchange=0.67,
    deep_ocean_efficacy=1.28,
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
                        - array of accumulated carbon in each atmospheric box,
                        - array of slow and fast temperature contributions,
                        - total accumulated carbon,
                        - emissions in the timestep before restart
                        - time constant scale factor
                        - CO2 concentrations in the timestep before restart
        restart_out   : if True, return the restart state as an extra output.
                        See restart_in.
        F_in          : either None, in which case calculate forcing from
                        simple Myhre logarithmic relationship, or numpy array
                        of prescribed total forcing
        temperature_function : 'Millar' (default) for impulse-response
                        relationship, 'Geoffroy' for two-layer energy balance
                        model
        ### the below parameters only affect 'Geoffroy' temperature function
        lambda_global : climate feedback parameter, W m-2 K-1
        ocean_heat_capacity : 2-element array of mixed layer and deep ocean
                        heat capacities (W yr m-2 K-1)
        ocean_heat_exchange : heat exchange coefficient between the two ocean
                        layers
        deep_ocean_efficacy : efficacy factor for deep ocean

    Outputs:
        E             : Timeseries of diagnosed CO2 emissions in GtC
        F             : Timeseries of total radiative forcing, W/m2
        T             : Timeseries of temperature anomaly since pre-industrial
        restart       : if restart_out=True, 6-tuple of carbon cycle state
                        parameters. See restart_in.
        ### the below outputs are included if 'Geoffroy' temperature function
        selected:
        lambda_eff    : effective climate feedback parameter (W m-2 K-1)
        ohc           : integrated ocean heat content (J)
        heatflux      : top of atmosphere energy imbalance, W m-2
    """
    
    # Error checking and validation goes here...

    # Dimensions
    nt = len(C)
    carbon_boxes_shape = (nt, a.shape[0])    

    # import correct conversion
    if temperature_function=='Millar':
        from .temperature.millar import forcing_to_temperature
        thermal_boxes_shape = (nt, d.shape[0])
    elif temperature_function=='Geoffroy':
        from .temperature.geoffroy import forcing_to_temperature
        thermal_boxes_shape = (nt, d.shape[0], 2)
        heatflux = np.zeros(nt)
        ohc = np.zeros(nt)
        lambda_eff = np.zeros(nt)
    else:
        raise ValueError('temperature_function must be "Millar" or "Geoffroy"')

    # Thermal response
    if type(tcrecs) is np.ndarray and temperature_function=='Millar':
        q     = calculate_q(tcrecs, d, F2x, tcr_dbl, nt)
    
    # Allocate intermediate and output arrays
    C_acc     = np.zeros(nt)
    R_i       = np.zeros(carbon_boxes_shape)
    emissions = np.zeros(nt)
    cumulative_emissions = np.zeros(nt)
    if F_in is None:
        F = np.zeros(nt)
        prescribed_forcing = False
    else:
        if np.isscalar(F_in):
            F_in = np.ones(nt) * F_in
        else:
            if len(F_in) != nt:
                raise ValueError('F_in must be same size as C, which is '+nt)
            F = F_in
        prescribed_forcing = True

    T_j       = np.zeros(thermal_boxes_shape)
    T         = np.zeros(nt)

    cumulative_emissions[0] = 0
    airborne_emissions = np.zeros_like(cumulative_emissions)

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
                C_pi, C_minus1, E_minus1,
            )
        )
        cumulative_emissions[0] = emissions[0]
        if not prescribed_forcing:
            F[0]          = co2_log(C[0], C_pi, F2x=F2x) + other_rf[0]
        if temperature_function=='Millar':
            T_j[0,:]      = forcing_to_temperature(T_j_minus1, q[0,:], d, F[0])
            T[0]=np.sum(T_j[0,:])
        else:
            # leave unimplemented unless somebody invents a use case
            raise(NotImplementedError('Restarts not implemented with Geoffroy '+
                'temperature function'))

    else:
        emissions[0]  = root(infer_emissions, 0., args=(C[0], R_i[0,:],
            tau, a, C_pi))['x']
        cumulative_emissions[0] = emissions[0]
        if not prescribed_forcing:
            F[0]          = co2_log(C[0], C_pi, F2x=F2x) + other_rf[0]
        if temperature_function=='Millar':
            T_j[0,:]      = forcing_to_temperature(T_j[0,:], q[0,:], d, F[0])
            T[0] = np.sum(T_j[0,:])
        else:
            T_j[0,:,:], heatflux[0], del_ohc, lambda_eff[0] = forcing_to_temperature(
                T_j[0,:,:],
                F[0],
                F[0],
                lambda_global=lambda_global,
                ocean_heat_capacity=ocean_heat_capacity,
                ocean_heat_exchange=ocean_heat_exchange,
                deep_ocean_efficacy=deep_ocean_efficacy,
                dt=1
            )
            T[0] = np.sum(T_j[0,:,:], axis=1)[0]
            ohc[0] = ohc[0] + del_ohc

    # Second timestep onwards
    for t in range(1,nt):
        emissions[t], C_acc[t], R_i[t,:], time_scale_sf = (
            inverse_carbon_cycle(
                C[t], C_acc[t-1], T[t-1], r0, rc, rt, iirf_max, 
                time_scale_sf, a, tau, iirf_h, R_i[t-1,:],
                C_pi, C[t-1], emissions[t-1]
            )
        )
        if not prescribed_forcing:
            F[t]     = co2_log(C[t], C_pi, F2x=F2x) + other_rf[t]
        if temperature_function=='Millar':
            T_j[t,:] = forcing_to_temperature(T_j[t-1,:], q[t,:], d, F[t])
            T[t] = np.sum(T_j[t,:])
        else:
            T_j[t,:,:], heatflux[t], del_ohc, lambda_eff[t] = forcing_to_temperature(
                T_j[t-1,:,:],
                F[t-1],
                F[t],
                lambda_global=lambda_global,
                ocean_heat_capacity=ocean_heat_capacity,
                ocean_heat_exchange=ocean_heat_exchange,
                deep_ocean_efficacy=deep_ocean_efficacy,
                dt=1
            )
            T[t] = np.sum(T_j[t,:,:], axis=1)[0]
            ohc[t] = ohc[t-1] + del_ohc


    cumulative_emissions = np.cumsum(emissions)
    airborne_emissions = np.sum(R_i, axis=1) * ppm_gtc

#    # Output temperatures
#    T = np.sum(T_j, axis=-1)
    if restart_out:
        restart_out_val = (R_i[-1], T_j[-1], C_acc[-1], emissions[-1],
            time_scale_sf, C[-1])
        return emissions, F, T, restart_out_val
    else:
        if temperature_function=='Millar':
            return emissions, F, T
        else:
            return emissions, F, T, lambda_eff, ohc, heatflux, airborne_emissions/cumulative_emissions
