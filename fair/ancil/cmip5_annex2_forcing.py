import numpy as np
import os

forcing_filename = os.path.join(os.path.dirname(__file__),
    'cmip5_annex2_forcing.csv')

class Forcing:
    forcing  = np.loadtxt(forcing_filename, skiprows=1, delimiter=',')
    year     = forcing[:,0]
    co2      = forcing[:,1]
    ghg_other= forcing[:,2]
    tropo3   = forcing[:,3]
    strato3  = forcing[:,4]
    aero     = forcing[:,5]
    landuse  = forcing[:,6]
    stwv     = forcing[:,7]
    bcsnow   = forcing[:,8]
    contrails= forcing[:,9]
    solar    = forcing[:,10]
    volcanic = forcing[:,11]
    total    = np.sum(forcing[:,1:], axis=1)
