import numpy as np

"""Default values of thermal response and radiative forcing parameters.

q       : slow and fast coefficients of contribution to temperature
tcr     : transient climate response (K)
ecs     : equilibrium climate sensitivity (K)
d       : slow and fast contributions to global mean temperature change (yr)
f2x     : radiative forcing from a doubling of CO2 (W/m2)
tcr_dbl : timescale for reaching a doubling of CO2 for TCR (yr)
"""

tcr     = 1.6
ecs     = 2.75

q       = np.array([0.33, 0.41])
tcrecs  = np.array([tcr, ecs])
d       = np.array([239.0,4.1])
f2x     = 3.71
tcr_dbl = 69.661 #TODO: v2.0 make explicitly log(2)/log(1.01)