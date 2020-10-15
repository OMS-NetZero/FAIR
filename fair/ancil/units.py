import numpy as np
import os

units_filename = os.path.join(os.path.dirname(__file__), 'units.csv')

class Units(dict):
    def __init__(self):
        alldata   = np.genfromtxt(units_filename, dtype = "str", skip_header=2, delimiter = "\t")
        self.update(zip([alldata[i,0].replace('"','') for i in range(len(alldata[:,0]))],\
                        [alldata[i,1].replace('"','') for i in range(len(alldata[:,1]))]))