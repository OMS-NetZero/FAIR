# usage
# from fair.ancil import historical_scaling
# s = historical_scaling.all

import os

import numpy as np

filename = os.path.join(os.path.dirname(__file__), 'historical_scaling.csv')
all = np.loadtxt(filename, skiprows=0, delimiter=',')
co2 = all[:,0]
