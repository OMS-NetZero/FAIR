# test that the carbon cycle provides expected results using default params
# and a 10Gt pulse of C

import fair
import numpy as np
import os

emissions = np.zeros(250)
emissions[125:] = 10.0
other_rf = np.zeros(emissions.size)
for x in range(0,emissions.size):
    other_rf[x] = 0.5*np.sin(2*np.pi*(x)/14.0)
    
C,F,T = fair.forward.fair_scm(emissions=emissions, other_rf=other_rf)

# here's one I prepared earlier
datadir = os.path.join(os.path.dirname(__file__), 'ten_GtC_pulse/')
C_expected = np.load(datadir + 'C.npy')
F_expected = np.load(datadir + 'F.npy')
T_expected = np.load(datadir + 'T.npy')

def test_ten_GtC_pulse():
    # use np.allclose rather than == to allow for small floating point errors
    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)
