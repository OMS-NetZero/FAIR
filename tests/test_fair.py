import pytest

import fair
import os
import numpy as np
from fair.RCPs import rcp3pd


def test_import():
    version = fair.__version__
    assert version[:3] == '1.1'


def test_no_arguments():
    with pytest.raises(ValueError):
        fair.forward.fair_scm()


def test_zero_emissions():
    nt = 250
    emissions = np.zeros(nt)
    C,F,T = fair.forward.fair_scm(emissions=emissions, other_rf=0.)
    assert np.allclose(C,np.ones(nt)*278.)
    assert np.allclose(F,np.zeros(nt))
    assert np.allclose(T,np.zeros(nt))


def test_ten_GtC_pulse():
    emissions = np.zeros(250)
    emissions[125:] = 10.0
    other_rf = np.zeros(emissions.size)
    for x in range(0,emissions.size):
        other_rf[x] = 0.5*np.sin(2*np.pi*(x)/14.0)
    
    C,F,T = fair.forward.fair_scm(emissions=emissions, other_rf=other_rf)
    
    datadir = os.path.join(os.path.dirname(__file__), 'ten_GtC_pulse/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp3pd():
    C,F,T = fair.forward.fair_scm(emissions=rcp3pd.Emissions.emissions)
    # Currently fails as useMultigas = False
    # Build in an exception in code to handle this.
