# usage
# from fair.ancil import historical_scaling
# s = historical_scaling.all

import numpy as np
import os

filename = os.path.join(os.path.dirname(__file__), 'historical_scaling.csv')
all = np.loadtxt(filename, skiprows=0, delimiter=',')
co2 = all[:,0]
