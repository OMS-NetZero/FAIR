import pytest

import fair
from fair.RCPs import rcp3pd, rcp45, rcp6, rcp85, rcp26, rcp60
import numpy as np
import os
from fair.tools import ensemble
from fair.constants import molwt

def test_ten_GtC_pulse():
    emissions = np.zeros(250)
    emissions[125:] = 10.0
    other_rf = np.zeros(emissions.size)
    for x in range(0,emissions.size):
        other_rf[x] = 0.5*np.sin(2*np.pi*(x)/14.0)

    C,F,T = fair.forward.fair_scm(
        emissions=emissions, other_rf=other_rf, useMultigas=False,
        r0=32.4, tcr_dbl=70)

    datadir = os.path.join(os.path.dirname(__file__), 'ten_GtC_pulse/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_multigas_fullemissions_error():
    with pytest.raises(ValueError):
        fair.forward.fair_scm(emissions=rcp3pd.Emissions.emissions,
            useMultigas=False)


# There must be a good way to avoid duplication here
def test_rcp3pd():
    C,F,T = fair.forward.fair_scm(
        emissions=rcp3pd.Emissions.emissions,
        b_aero = np.array([-35.29e-4*1.3741*molwt.SO2/molwt.S, 0.0, -5.034e-4*1.3741, -5.763e-4*1.3741*molwt.NO/molwt.N, 453e-4*1.3741,-37.83e-4*1.3741, -10.35e-4*1.3741]),
        efficacy=np.ones(13)
    )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp3pd/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp45():
    C,F,T = fair.forward.fair_scm(
        emissions=rcp45.Emissions.emissions,
        b_aero = np.array([-35.29e-4*1.3741*molwt.SO2/molwt.S, 0.0, -5.034e-4*1.3741, -5.763e-4*1.3741*molwt.NO/molwt.N, 453e-4*1.3741,-37.83e-4*1.3741, -10.35e-4*1.3741]),
        efficacy=np.ones(13)
    )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp45/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp6():
    C,F,T = fair.forward.fair_scm(
        emissions=rcp6.Emissions.emissions,
        b_aero = np.array([-35.29e-4*1.3741*molwt.SO2/molwt.S, 0.0, -5.034e-4*1.3741, -5.763e-4*1.3741*molwt.NO/molwt.N, 453e-4*1.3741,-37.83e-4*1.3741, -10.35e-4*1.3741]),
        efficacy=np.ones(13)
    )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp6/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp85():
    C,F,T = fair.forward.fair_scm(
        emissions=rcp85.Emissions.emissions,
        b_aero = np.array([-35.29e-4*1.3741*molwt.SO2/molwt.S, 0.0, -5.034e-4*1.3741, -5.763e-4*1.3741*molwt.NO/molwt.N, 453e-4*1.3741,-37.83e-4*1.3741, -10.35e-4*1.3741]),
        efficacy=np.ones(13)
    )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp85/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    print(np.sum(F, axis=1))

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


# rcp3pd and rcp6 have been renamed. The modules should still work otherwise
# the tests would not have got to this point. But we import directly here to
# ensure compatibility.
def test_rcp_aliases():

    # 1. rcp26
    C,F,T = fair.forward.fair_scm(
        emissions=rcp26.Emissions.emissions,
        b_aero = np.array([-35.29e-4*1.3741*molwt.SO2/molwt.S, 0.0, -5.034e-4*1.3741, -5.763e-4*1.3741*molwt.NO/molwt.N, 453e-4*1.3741,-37.83e-4*1.3741, -10.35e-4*1.3741]),
        efficacy=np.ones(13)
    )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp3pd/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected) 

    # 2. rcp60
    C,F,T = fair.forward.fair_scm(
        emissions=rcp60.Emissions.emissions,
        b_aero = np.array([-35.29e-4*1.3741*molwt.SO2/molwt.S, 0.0, -5.034e-4*1.3741, -5.763e-4*1.3741*molwt.NO/molwt.N, 453e-4*1.3741,-37.83e-4*1.3741, -10.35e-4*1.3741]),
        efficacy=np.ones(13)
    )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp6/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_ensemble():
    tcrecs = ensemble.tcrecs_generate('cmip5', dist='lognorm', n=100,
      seed=40000)
    # check first 100 lines
    filepath = os.path.join(os.path.dirname(__file__), 
      'tcrecs/tcrecs_output.txt')
    tcrecs_compare = np.loadtxt(filepath)
    assert np.allclose(tcrecs,tcrecs_compare)


def test_co2_concentration_driven():
    C, F, T = fair.forward.fair_scm(
        emissions_driven=False,
        C=rcp45.Concentrations.co2,
        useMultigas=False
        )
    assert (C==rcp45.Concentrations.co2).all()
    datadir = os.path.join(os.path.dirname(__file__), 'rcp45/')
    T_expected = np.load(datadir + 'T_concdriven.npy')
    assert np.allclose(T, T_expected)


def test_multigas_concentration_driven():
    C, F, T = fair.forward.fair_scm(
        emissions_driven=False,
        C=rcp45.Concentrations.gases,
        F_tropO3 = rcp45.Forcing.tropo3,
        F_aerosol = rcp45.Forcing.aero+rcp45.Forcing.cloud,
        F_bcsnow = rcp45.Forcing.bcsnow,
        useMultigas=True
        )
    datadir = os.path.join(os.path.dirname(__file__), 'rcp45/')
    T_expected = np.load(datadir + 'T_concdriven_multi.npy')
    assert np.allclose(T, T_expected)
