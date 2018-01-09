import numpy as np

def cumulative(emissions, aCO2land=-0.00113789):
    E_CO2land = emissions[:,2]
    return np.cumsum(E_CO2land) * aCO2land
