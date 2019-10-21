import pytest

import fair
from fair.RCPs import rcp26, rcp3pd, rcp45, rcp6, rcp60, rcp85
import numpy as np
import warnings
from fair.forcing.ghg import myhre, etminan
from fair.constants import molwt
from copy import deepcopy

# Rather than change the testing values, I am going to check the old results
# can still be recreated with non-default options.

# This means that all multi-gas scenarios will have aerosol scalings that
# differ from the defaults.

def test_zero_emissions():
    nt = 250
    emissions = np.zeros(nt)
    C,F,T = fair.forward.fair_scm(
        emissions=emissions, other_rf=0., useMultigas=False)
    assert np.allclose(C,np.ones(nt)*278.)
    assert np.allclose(F,np.zeros(nt))
    assert np.allclose(T,np.zeros(nt))


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

    
def test_contrails_ext():
    # impose contrail forcing of +0.2 W/m2 and scale factor of 1.5; verify that
    # contrail forcing of 0.3 W/m2 returned
    C,F,T = fair.forward.fair_scm(emissions=rcp45.Emissions.emissions,
        contrail_forcing='ext', F_contrails=0.2, scale=np.ones(13)*1.5)
    assert np.allclose(F[:,7],np.ones(736)*0.3)


def test_contrails_invalid():
    with pytest.raises(ValueError):
        fair.forward.fair_scm(emissions=rcp45.Emissions.emissions,
            contrail_forcing='other')


def test_aerosol_regression_zeros_fix():
    _, F, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        b_aero = np.zeros(7),
        aerosol_forcing='aerocom'
    )
    # Index 8 is aerosol forcing
    assert (F[:,8]==np.zeros(736)).all()


def test_aerosol_regression_zeros_nofix():
    _, F, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        b_aero = np.zeros(7),
        fixPre1850RCP=False,
        aerosol_forcing='aerocom'
    )
    # Index 7 is aerosol forcing
    assert (F[:,8]==np.zeros(736)).all()


def test_aerosol_noscale():
    _, F1, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        scaleAerosolAR5=False
    )
    _, F2, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        scaleAerosolAR5=True
    )
    assert (F1[:,8]!=F2[:,8]).any()


def test_stevens():
    # The Stevens aerosol causes an overflow in the root finder routine at
    # present. This isn't a problem, but switching it off gives nice clean
    # tests
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, F1, _ = fair.forward.fair_scm(
            emissions=fair.RCPs.rcp85.Emissions.emissions,
            aerosol_forcing='Stevens'
        )
    _, F2, _ = fair.forward.fair_scm(
        emissions=fair.RCPs.rcp85.Emissions.emissions,
        aerosol_forcing='aerocom+ghan'
    )
    assert (F1[:,8]!=F2[:,8]).any()


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
    emissions = deepcopy(fair.RCPs.rcp85.Emissions.emissions)
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
@pytest.mark.xfail(msg='I do not yet know why this is failing')
def test_co2_scale():
    emissions = deepcopy(fair.RCPs.rcp85.Emissions.emissions)
    scale = np.ones(13)
    scale[0] = 1.15
    C1, F1, T1 = fair.forward.fair_scm(emissions)
    C2, F2, T2 = fair.forward.fair_scm(emissions, scale=scale)
    assert (C2 == C1).all()
    assert (F2 == F1).all()
    assert (T2 == T1).all()


def test_etminan_noscale():
    C = [350, 1000, 500]
    Cpi = np.array([278., 722., 273.])
    F1 = etminan(C, Cpi, scale_F2x=False)
    F2 = etminan(C, Cpi, scale_F2x=True)
    assert (F1[0] != F2[0])


def test_myhre():
    C = [350, 1000, 500]
    Cpi = np.array([278., 722., 273.])
    rf = myhre(C, Cpi)
    assert np.allclose(rf, np.array([1.232721, 0.150752, 0.659443]))


def test_landuse_ext():
    # we impose a land-use ERF of -0.1 W/m2 and a forcing scale factor if 1.2
    # land use forcing should be -0.12 W/m2
    C,F,T = fair.forward.fair_scm(emissions=rcp45.Emissions.emissions,
        landuse_forcing='external', scale=np.ones(13)*1.2, F_landuse=-0.1)
    assert np.allclose(F[:,10],np.ones(736)*(-0.12))

    
def test_contrails_nox():
    # default option, should not break anything
    C,F,T = fair.forward.fair_scm(emissions=rcp45.Emissions.emissions,
        contrail_forcing='NOx')


def test_contrails_fuel():
    kerosene_supply = np.array([0.]*175 + [8.9534883721,
        9.6511627907, 10.4651162791, 11.2790697674, 12.2093023256,
        13.1395348837, 14.3023255814, 15.3488372093, 16.6279069767,
        17.9069767442, 19.4186046512, 20.9302325581, 22.5581395349,
        24.4186046512, 26.3953488372, 28.488372093, 30.8139534884,
        33.2558139535, 35.9302325581, 38.7209302326, 41.8604651163,
        45.9302325581, 50.5813953488, 53.3720930233, 55.8139534884,
        59.6511627907, 64.6511627907, 76.2790697674, 86.3953488372,
        90.4651162791, 90.6976744186, 104.6511627907, 111.6279069767,
        115.5813953488, 111.6279069767, 111.7441860465, 112.0930232558,
        118.7209302326, 122.9069767442, 128.023255814, 128.9534883721,
        127.0930232558, 128.488372093, 130.2325581395, 138.9534883721,
        143.488372093, 151.0465116279, 157.6744186047, 164.4186046512,
        170.3488372093, 170.8139534884, 166.7441860465, 165.1162790698,
        167.5581395349, 174.4186046512, 179.4186046512, 197.349, 204.095,
        205.029, 211.022, 218.996, 213.021, 214.353, 215.012, 228.23, 236.063,
        238.329, 244.104, 247.726, 239.072, 246.544, 254.792])
    C,F,T = fair.forward.fair_scm(emissions=rcp45.Emissions.emissions[:247,:],
        contrail_forcing='fuel', kerosene_supply=kerosene_supply, natural=0.,
        F_volcanic=0., F_solar=0.)
    F_expected = np.array([
        0.        , 0.00169919, 0.0018316 , 0.00198607, 0.00214054, 0.00231708,
        0.00249362, 0.00271429, 0.0029129 , 0.00315564, 0.00339838, 0.00368526,
        0.00397214, 0.00428108, 0.00463416, 0.00500931, 0.00540652, 0.00584787,
        0.00631128, 0.00681883, 0.00734845, 0.00794427, 0.00871663, 0.00959933,
        0.01012895, 0.01059236, 0.01132059, 0.01226949, 0.01447623, 0.0163961 ,
        0.01716846, 0.01721259, 0.01986068, 0.02118473, 0.02193502, 0.02118473,
        0.02120679, 0.021273  , 0.02253084, 0.02332527, 0.02429623, 0.02447277,
        0.02411969, 0.0243845 , 0.02471551, 0.02637057, 0.0272312 , 0.02866558,
        0.02992343, 0.03120334, 0.03232878, 0.03241705, 0.03164469, 0.03133574,
        0.03179916, 0.03310114, 0.03405004, 0.03745286, 0.03873312, 0.03891037,
        0.04004772, 0.04156103, 0.04042709, 0.04067988, 0.04080494, 0.04331345,
        0.0448    , 0.04523004, 0.04632602, 0.0470134 , 0.04537105, 0.04678908,
        0.04835439])
    assert np.allclose(F[174:,7], F_expected)


def test_co2only_scale_error():
    """Test to fix #42"""
    C, F, T = fair.forward.fair_scm(
        emissions = rcp85.Emissions.co2,
        useMultigas=False,
        scale=1.0)
