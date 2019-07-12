from __future__ import division

import numpy as np
from ..constants import molwt, radeff
from ..constants.general import M_ATMOS

def gwp(H, tau, re, wt, f=0., a=np.array([0.2173, 0.2240, 0.2824, 0.2763]),
    tau_co2=np.array([np.inf, 394.4, 36.54, 4.304])):

    """Calculates Global Warming Potentials

    Inputs:
        H: Time horizon (years)
        tau: lifetime of gas
        re: radiative efficiency of gas
        wt: molecular weight of gas

    Keywords:
        f: feedback factor. f=0.65 for methane and -0.071874 for N2O
        a: impulse response partition for CO2
        tau_co2: lifetimes for impulse response partition for CO2

    Outputs:
        global warming potential, scalar

    TODO: recognise gas from list and select tau, re, wt, f accordingly
    """

    agwp_co2 = radeff.CO2 * (
      a[0]*H + np.sum(a[1:]*tau_co2[1:]*(1-np.exp(-H/tau_co2[1:])))) * (
      molwt.AIR/molwt.CO2 * 1e9/M_ATMOS)

    agwp_gas = re * (1.0+f) * tau * (1.0 - np.exp(-H/tau)) * (
      molwt.AIR/wt * 1e9/M_ATMOS)

    gwp = agwp_gas/agwp_co2
    return gwp
