import numpy as np
import os

emissions_filename = os.path.join(os.path.dirname(__file__), 'natural.csv')

class Emissions:
    alldata   = np.loadtxt(emissions_filename, skiprows=4)
    year      = alldata[:,0]
    ch4       = alldata[:,1]
    n2o       = alldata[:,2]
    emissions = alldata[:,1:]
