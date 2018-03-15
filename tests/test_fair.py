import pytest

import fair
import os
import numpy as np
from fair.RCPs import rcp3pd, rcp45, rcp6, rcp85
from fair.tools import magicc
from fair.ancil import natural
from fair.constants import molwt
from fair.forcing.ghg import myhre

def test_no_arguments():
    with pytest.raises(ValueError):
        fair.forward.fair_scm()


def test_zero_emissions():
    nt = 250
    emissions = np.zeros(nt)
    C,F,T = fair.forward.fair_scm(
        emissions=emissions, other_rf=0., useMultigas=False)
    assert np.allclose(C,np.ones(nt)*278.)
    assert np.allclose(F,np.zeros(nt))
    assert np.allclose(T,np.zeros(nt))


# Rather than change the testing values, I am going to check the old results
# can still be recreated with non-default options.


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
    C,F,T = fair.forward.fair_scm(emissions=rcp3pd.Emissions.emissions)
    datadir = os.path.join(os.path.dirname(__file__), 'rcp3pd/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp45():
    C,F,T = fair.forward.fair_scm(emissions=rcp45.Emissions.emissions)
    datadir = os.path.join(os.path.dirname(__file__), 'rcp45/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp6():
    C,F,T = fair.forward.fair_scm(emissions=rcp6.Emissions.emissions)
    datadir = os.path.join(os.path.dirname(__file__), 'rcp6/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_rcp85():
    C,F,T = fair.forward.fair_scm(emissions=rcp85.Emissions.emissions)
    datadir = os.path.join(os.path.dirname(__file__), 'rcp85/')
    C_expected = np.load(datadir + 'C.npy')
    F_expected = np.load(datadir + 'F.npy')
    T_expected = np.load(datadir + 'T.npy')

    assert np.allclose(C, C_expected)
    assert np.allclose(F, F_expected)
    assert np.allclose(T, T_expected)


def test_division():
    # Ensure parameters given as integers are treated as floats when dividing
    # (Python2 compatibility).
    _, _, T = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp6.Emissions.emissions,
        useMultigas=True,
        d=np.array([239.0, 4.0]),
        tcr_dbl=70.0
    )
    _, _, T_int_params = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp6.Emissions.emissions,
        useMultigas=True,
        d=np.array([239, 4]),
        tcr_dbl=70
    )
    assert (T == T_int_params).all()


def test_scenfile():
    datadir = os.path.join(os.path.dirname(__file__), 'rcp45/')
    # Purpose of this test is to determine whether the SCEN file for RCP4.5
    # which does not include CFCs, years before 2000, or emissions from every
    # year from 2000 to 2500, equals the emissions file from RCP4.5
    # after reconstruction.
    # The .SCEN and .XLS files at http://www.pik-potsdam.de/~mmalte/rcps
    # sometimes differ in the 4th decimal place. Thus we allow a tolerance of
    # 0.0002 in addition to machine error in this instance.

    E1 = magicc.scen_open(datadir + 'RCP45.SCEN')
    assert np.allclose(E1, rcp45.Emissions.emissions, rtol=1e-8, atol=2e-4)

    E2 = magicc.scen_open(datadir + 'RCP45.SCEN', include_cfcs=False)
    assert np.allclose(E2, rcp45.Emissions.emissions[:,:24], rtol=1e-8,
        atol=2e-4)

    E3 = magicc.scen_open(datadir + 'RCP45.SCEN', startyear=1950)
    assert np.allclose(E3, rcp45.Emissions.emissions[185:,:], rtol=1e-8,
        atol=2e-4)

    # Test `_import_emis_file`
    assert magicc._import_emis_file('rcp26') == rcp3pd.Emissions
    assert magicc._import_emis_file('rcp6') == rcp6.Emissions
    assert magicc._import_emis_file('rcp85') == rcp85.Emissions
    with pytest.raises(ValueError):
        magicc._import_emis_file('rcp19')

    test_files = os.path.join(os.path.dirname(__file__), "scenfiles")
    scenfile_2000 = os.path.join(test_files, "WORLD_ONLY.SCEN")
    scenfile_2010 = os.path.join(test_files, "WORLD_ONLY_2010.SCEN")

    # Test CFCs inclusion.
    with pytest.raises(ValueError):
        magicc.scen_open(datadir + 'RCP45.SCEN', include_cfcs=np.zeros(0))
    E4 = magicc.scen_open(
            scenfile_2000, startyear=2000, include_cfcs=np.ones((51, 16)))
    assert E4[0, -1] == 1
    with pytest.raises(ValueError):
        magicc.scen_open(datadir + 'RCP45.SCEN', include_cfcs="foo")

    # Test filling of history and harmonisation.
    with pytest.raises(ValueError):
        magicc.scen_open(scenfile_2010)
    with pytest.raises(ValueError):
        magicc.scen_open(scenfile_2000, harmonise=1950)
    with pytest.raises(ValueError):
        magicc.scen_open(scenfile_2000, harmonise=2060)
    E5 = magicc.scen_open(scenfile_2000, harmonise=2010)
    assert E5[0, 1] == rcp45.Emissions.co2_fossil[0]


def test_strat_h2o_scale_factor():
    # Default scale factor changed to 0.12 for Etminan but can be overridden
    _, F1, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        useMultigas=True,
        ghg_forcing='Etminan'
    )

    _, F2, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        useMultigas=True,
        ghg_forcing='Etminan',
        stwv_from_ch4=0.15
    )
    # Index 6 is stratospheric water vapour. Check latter is 80% of former
    assert np.allclose(F1[:,6],F2[:,6]*0.8)


def test_aerosol_regression_zeros_fix():
    _, F, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        b_aero = np.zeros(7)
    )
    # Index 7 is aerosol forcing
    assert (F[:,7]==np.zeros(736)).all()


def test_aerosol_regression_zeros_nofix():
    _, F, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        b_aero = np.zeros(7),
        fixPre1850RCP=False
    )
    # Index 7 is aerosol forcing
    assert (F[:,7]==np.zeros(736)).all()


def test_ozone_regression_zero():
    _, F, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        useStevenson=False,
        b_tro3 = np.zeros(4)
    )
    # Index 4 is ozone forcing
    assert (F[:,4]==np.zeros(736)).all()


def test_timevarying_ecs():
    emissions = np.zeros(250)
    emissions[125:] = 10.0
    tcr = np.ones(250)*1.6
    ecs = np.linspace(1.5,4.5,num=250)
    tcrecs = np.vstack((tcr,ecs)).T
    C,F,T = fair.forward.fair_scm(
        emissions=emissions,
        tcrecs=tcrecs,
        useMultigas=False
    )


def test_ozone_stevenson_zero_nofix():
    emissions = fair.RCPs.rcp85.Emissions.emissions
    # zero all emissions except methane which we fix to pi steady state and
    # CO, NMVOC and NOx which we fix to 1750 anthro values from Skeie et al
    emissions[:,1:] = 0
    emissions[:,3] = 209.279889
    emissions[:,6] = 170.
    emissions[:,7] = 10.
    emissions[:,8] = 4.29 * molwt.N / molwt.NO
    _, F, _ = fair.forward.fair_scm(
        emissions=emissions,
        natural=0,
        fixPre1850RCP=False,
        useTropO3TFeedback=False)

    # won't be exactly zero due to numerical precision in eyeballing natural
    # CH4 to balance
    assert np.allclose(F[:,4],np.zeros(736))


# Test if changing the scale factor for CO2 forcing has no effect on concs,
# temperature and forcing (this is desired: change F2x to change CO2 forcing)
def test_co2_scale():
    emissions = fair.RCPs.rcp85.Emissions.emissions
    scale = np.ones(13)
    scale[0] = 1.15
    C1, F1, T1 = fair.forward.fair_scm(emissions)
    C2, F2, T2 = fair.forward.fair_scm(emissions, scale=scale)
    assert (C2 == C1).all()
    assert (F2 == F1).all()
    assert (T2 == T1).all()


def test_myhre():
    C = [350, 1000, 500]
    Cpi = np.array([278., 722., 273.])
    rf = myhre(C, Cpi)
    assert np.allclose(rf, np.array([1.232721, 0.150752, 0.659443]))
