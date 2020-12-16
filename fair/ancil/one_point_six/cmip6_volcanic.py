import numpy as np
import os

forcing_filename = os.path.join(os.path.dirname(__file__), 'cmip6_volcanic.csv')

class Forcing:
    forcing  = np.loadtxt(forcing_filename, skiprows=9, delimiter=',')
    year     = forcing[:,0]
    volcanic = forcing[:,1]
