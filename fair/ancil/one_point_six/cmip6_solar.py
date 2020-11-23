import numpy as np
import os

forcing_filename = os.path.join(os.path.dirname(__file__), 'cmip6_solar.csv')

class Forcing:
    forcing  = np.loadtxt(forcing_filename, skiprows=7, delimiter=',')
    year     = forcing[:,0]
    solar    = forcing[:,1]
