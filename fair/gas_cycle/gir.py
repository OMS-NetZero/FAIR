from __future__ import division

import numpy as np

from ..constants.general import ppm_gtc

"""Gas cycle functions from Generalised Impulse Response Model v1.0.0.

Much of this has been adapted from:

Leach et al., 2020, Geoscientific Model Development
https://www.geosci-model-dev-discuss.net/gmd-2019-379/
"""

def calculate_alpha(cumulative_emissions,airborne_emissions,temperature,r0,rC,rT,g0,g1,iirf_max = 97.0):
    """
    Calculate CO2 time constant scaling factor.

    Inputs:
        cumulative_emissions: GtC cumulative emissions since pre-industrial.
        airborne_emissions: GtC total emissions remaining in the atmosphere.
        temperature: K temperature anomaly since pre-industrial.
        r0: pre-industrial 100-year time-integrated airborne fraction.
        rC: sensitivity of 100-year time-integrated airborne fraction with 
            atmospheric carbon stock.
        rT: sensitivity of 100-year time-integrated airborne fraction with
            temperature anomaly.
        g0: parameter for alpha
        g1: parameter for alpha
    Keywords:
        iirf_max: maximum allowable value to 100-year time-integrated airborne
            fraction
    Outputs:
        alpha: scaling factor.
    """
    iirf = r0 + rC * (cumulative_emissions-airborne_emissions) + rT * temperature
    iirf = (iirf>iirf_max) * iirf_max + iirf * (iirf<iirf_max)
    alpha = g0 * np.sinh(iirf / g1)
    return alpha

def step_concentration(carbon_boxes0,emissions,alpha,a,tau,Cpi,dt=1):
    """
    Calculate concentrations from emissions.

    Inputs:
        carbon_boxes0: CO2 boxes at the end of the previous timestep.
        emissions: GtC CO2 emissions this timestep.
        alpha: CO2 time constant scaling factor.
        a: CO2 partitioning coefficient
        tau: CO2 atmospheric time constants (unscaled).
        Cpi: pre-industrial CO2 concentrations (ppm).
    Keywords:
        dt: timestep in years.
    Outputs:
        C: CO2 concentration (ppm)
        carbon_boxes1: CO2 boxes at the end of this timestep.
        airbone_emissions: GtC total emissions remaining in atmosphere.
    """
    carbon_boxes1 = emissions / ppm_gtc * a * alpha * (tau/dt) * (
        1. - np.exp(-dt/(alpha*tau))) + carbon_boxes0 * np.exp(-dt/(alpha*tau))
    C = Cpi + np.sum(carbon_boxes1 + carbon_boxes0) / 2
    airborne_emissions = np.sum(carbon_boxes1) * ppm_gtc
    return C, carbon_boxes1, airborne_emissions
