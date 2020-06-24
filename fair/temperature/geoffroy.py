from __future__ import division

import numpy as np

from ..constants.general import EARTH_RADIUS, SECONDS_PER_YEAR

# This is a stop-gap until I figure out how to couple the openscm-twolayermodel

def forcing_to_temperature(
    temp,
    f0,
    f1,
    lambda_global=1.18,
    ocean_heat_capacity=np.array([8.2, 109.0]),
    ocean_heat_exchange=0.67,
    deep_ocean_efficacy=1.28,
    dt=1):
    """Calculate temperature from a given radiative forcing.

    This follows the forcing to temperature function described in Geoffroy
    et al. (2013a, 2013b).

    Thanks to Glen Harris, UK Met Office, for doing the grunt work.

    Inputs:
        temp: (2,2) numpy array (layer, component) of ocean temperatures
            in timestep t-1 (layer = mixed, deep); (component=fast; slow)
        f0: effective radiative forcing in timestep t-1
        f1: effective radiative forcing in timestep t

    Keywords:
        lambda_global: Climate feedback parameter (convention positive
            stable), W m-2 K-1
        ocean_heat_capacity: 2-element np.ndarray of ocean heat capacities
            (mixed layer, deep ocean) in W yr m-2 K-1
        ocean_heat_exchange: heat exchange coefficient between the two
            ocean layers in W m-2 K-1
        deep_ocean_efficacy: deep ocean efficacy parameter, non-dimensional           
        dt: timestep, year

    Outputs:
        t1: slow and fast contributions to total temperature (2 element array)
        in timestep t
    """

    # care with unit handling. This will be a large number - divide 1e21 for ZJ
    ntoa_joule = 4 * np.pi * EARTH_RADIUS**2 * SECONDS_PER_YEAR
             
    # Define derived constants                
    cdeep_p = ocean_heat_capacity[1] * deep_ocean_efficacy
    gamma_p = ocean_heat_exchange * deep_ocean_efficacy
    g1 = (lambda_global+gamma_p)/ocean_heat_capacity[0]
    g2      = gamma_p/cdeep_p
    g       = g1+g2
    gstar   = g1-g2
    delsqrt = np.sqrt(g*g - 4*g2*lambda_global/ocean_heat_capacity[0])
    afast   = (g + delsqrt)/2
    aslow   = (g - delsqrt)/2
    cc      = 0.5/(ocean_heat_capacity[0]*delsqrt)        
    amix_f  = cc*(gstar+delsqrt)
    amix_s  = -cc*(gstar-delsqrt)        
    adeep_f = -gamma_p/(ocean_heat_capacity[0]*cdeep_p*delsqrt)
    adeep_s = -adeep_f                                            

    temp_mix0 = np.copy(temp[0,:])
    temp_deep0 = np.copy(temp[1,:])

    adf = 1/(afast*dt)
    ads = 1/(aslow*dt)
    exp_f = np.exp(-1./adf)
    exp_s = np.exp(-1./ads)        
    int_f = (f0*adf + f1*(1-adf) - exp_f*(f0*(1+adf)-f1*adf))/afast
    int_s = (f0*ads + f1*(1-ads) - exp_s*(f0*(1+ads)-f1*ads))/aslow

    temp_mix1  = np.array([exp_f*temp_mix0[0] + amix_f*int_f, 
                           exp_s*temp_mix0[1] + amix_s*int_s]) 
                                         
    temp_deep1 = np.array([exp_f*temp_deep0[0] + adeep_f*int_f, 
                           exp_s*temp_deep0[1] + adeep_s*int_s])

    c_dtemp = (
        ocean_heat_capacity[0]*(temp_mix1.sum()-temp_mix0.sum()) + 
        ocean_heat_capacity[1]*(temp_deep1.sum()-temp_deep0.sum())
    )

    heatflux = c_dtemp/dt
    del_ohc  = ntoa_joule * c_dtemp
        
    factor_lambda_eff = (deep_ocean_efficacy-1.0)*ocean_heat_exchange
    if abs(np.sum(temp_mix1)) > 1e-6:
        ratio = (np.sum(temp_mix1) - np.sum(temp_deep1))/np.sum(temp_mix1)
        lambda_eff = lambda_global + factor_lambda_eff*ratio
    else:
        lambda_eff = lambda_global + factor_lambda_eff
        
    return np.array([temp_mix1, temp_deep1]), heatflux, del_ohc, lambda_eff

