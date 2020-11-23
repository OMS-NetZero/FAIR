import numpy as np

"""Default values of carbon cycle parameters.

r0      : pre-industrial time-integrated airborne fraction of CO2 (yr)
rc      : change in time-integrated airborne fraction with CO2 (yr/GtC)
rt      : change in time-integrated airborne fraction with temperature (yr/K)
iirf_h  : time horizon for time integrated airborne fraction (yr)
iirf_max: maximum allowed value of iirf (yr)
a       : partition coefficient of carbon boxes
tau     : present-day decay time constants of CO2 (yr)
"""

r0       = 35.0
rc       = 0.019
rt       = 4.165
iirf_h   = 100.0
iirf_max = 97.0
a        = np.array([0.2173,0.2240,0.2824,0.2763])
tau      = np.array([1000000,394.4,36.54,4.304])