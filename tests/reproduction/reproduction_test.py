import pytest

import fair
from fair.RCPs import rcp3pd, rcp45, rcp6, rcp85, rcp26, rcp60
import numpy as np
import os
from fair.constants import molwt, radeff, lifetime
from fair.tools.constrain import hist_temp
from fair.tools.gwp import gwp

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


def test_inverse_fair():
    """Tests reproducibility of concentrations-to-emissions FaIR."""

    # initialise a 1% run
    nt = 140
    C = 1.01**np.arange(nt)*278.

    E,F,T = fair.inverse.inverse_fair_scm(C=C, tcrecs=np.array([1.7, 3.0]))

    datadir = os.path.join(os.path.dirname(__file__), '1pctCO2/')
    E_expected = np.load(datadir + 'E.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(E, E_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_forward_versus_reverse():
    """Does inverse FaIR recover the same emissions as forward FaIR?

    Both methods require numerical root finding methods so exact correspondence
    is quite unlikely, so accept a small tolerance"""

    E_forward = rcp85.Emissions.co2
    other_rf = np.sin(np.arange(736)) * 0.2
    C_forward, F_forward, T_forward = fair.forward.fair_scm(emissions=E_forward, other_rf=other_rf, useMultigas=False)
    E_inverse, F_inverse, T_inverse = fair.inverse.inverse_fair_scm(C=C_forward, other_rf=other_rf)

    assert np.allclose(E_forward, E_inverse, atol=0.01, rtol=0.01)
    assert np.allclose(F_forward, F_inverse, atol=0.01, rtol=0.01)
    assert np.allclose(T_forward, T_inverse, atol=0.01, rtol=0.01)


def test_restart_co2_continuous():
    """Tests to check that a CO2-only run with a restart produces the same
    results as a CO2-only run without a restart."""

    C, F, T = fair.forward.fair_scm(
        emissions   = rcp45.Emissions.co2[:20],
        useMultigas = False
        )

    C1, F1, T1, restart = fair.forward.fair_scm(
        emissions   = rcp45.Emissions.co2[:10],
        useMultigas = False,
        restart_out = True
        )

    C2, F2, T2 = fair.forward.fair_scm(
        emissions   = rcp45.Emissions.co2[10:20],
        useMultigas = False,
        restart_in  = restart
        )

    assert np.all(C == np.concatenate((C1, C2)))
    assert np.all(F == np.concatenate((F1, F2)))
    assert np.all(T == np.concatenate((T1, T2)))


def test_inverse_restart():
    """Tests restarts for inverse FaIR."""

    E, F, T = fair.inverse.inverse_fair_scm(
        C = rcp85.Concentrations.co2[:20])

    E1, F1, T1, restart = fair.inverse.inverse_fair_scm(
        C = rcp85.Concentrations.co2[:10], restart_out=True)

    E2, F2, T2 = fair.inverse.inverse_fair_scm(
        C = rcp85.Concentrations.co2[10:20], restart_in=restart)

    assert np.all(E == np.concatenate((E1, E2)))
    assert np.all(F == np.concatenate((F1, F2)))
    assert np.all(T == np.concatenate((T1, T2)))


def test_constrain():
    """Checks that the historical temperature constraining function works"""

    datadir = os.path.join(os.path.dirname(__file__),
        '../../fair/tools/tempobs/')
    tempobsdata = np.loadtxt(datadir+'had4_krig_annual_v2_0_0.csv')
    years   = tempobsdata[:,0]
    tempobs = tempobsdata[:,1]

    C,F,T = fair.forward.fair_scm(emissions=rcp45.Emissions.emissions)
    accept,_,_,_,_ = hist_temp(tempobs, T[85:252], years)
    assert accept==True

    accept,_,_,_,_ = hist_temp(tempobs, np.zeros(167), years)
    assert accept==False


def test_gwp():
    """Checks that GWP calculator produces correct GWPs."""

    # methane uses "perturbation lifetime" for GWP calculations and feedback
    # factor
    assert np.round(gwp(100, 12.4, radeff.CH4, molwt.CH4, f=0.65))==28

    # for N2O, I think the IPCC AR5 value is out by one year. Most likely
    # explanation is that they have rounded off the feedback somewhere.
    # This is calculated as 1-(1-0.36*(1.65)*radeff.CH4/radeff.N2O). See
    # eq. 8.SM.20 in the supplement to Chapter 8, AR5
    assert np.round(
        gwp(20, lifetime.N2O, radeff.N2O, molwt.N2O, f=-0.071874))==263
    assert np.round(
        gwp(100, lifetime.N2O, radeff.N2O, molwt.N2O, f=-0.071874))==264

    # Now check a nice straightforward example
    assert np.round(
        gwp(100, lifetime.CFC11, radeff.CFC11, molwt.CFC11), decimals=-1)==4660
