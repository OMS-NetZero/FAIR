import pytest

import numpy as np

import concentration_driven

def test_single_species():
    #emissions are based on example CO2 concentrations
    out_dict = concentration_driven._run_numpy( inp_ar = np.array([[399.6173423,\
                                                                    401.8213271,\
                                                                    406.3410747,\
                                                                    415.7639971,\
                                                                    435.6164218]]),
                                            a1 = np.array([0.2173]),
                                            a2 = np.array([0.224]),
                                            a3 = np.array([0.2824]),
                                            a4 = np.array([0.2763]),
                                            tau1 = np.array([1000000]),
                                            tau2 = np.array([394.4]),
                                            tau3 = np.array([36.54]),
                                            tau4 = np.array([4.304]),
                                            r0 = np.array([28.627296]),
                                            rC = np.array([0.019773]),
                                            rT = np.array([4.334433]),
                                            rA = np.array([0]),
                                            PI_conc = np.array([278]),
                                            emis2conc = np.array([0.468952343952344]),
                                            f1 = np.array([5.754389]),
                                            f2 = np.array([0.001215]),
                                            f3 = np.array([-0.069598]),
                                            d = np.array([  283,\
                                                            9.88,\
                                                            0.85]),
                                            q = np.array([  0.311333,\
                                                            0.165417,\
                                                            0.242]),
                                            ext_forcing = np.array([-0.256119925,\
                                                                    -0.304324144,\
                                                                    -0.501633962,\
                                                                    -0.904262779,\
                                                                    -0.12278275]),
                                            timestep = np.array([1,\
                                                                2,\
                                                                4,\
                                                                8,\
                                                                1])
                                          )
    emissions_out = out_dict["emissions"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]

    emissions_compare = np.array([[ 319.828087,\
                            46.772580,\
                            23.792089,\
                            18.293743,\
                            18.238578]]) 
    T_compare = np.array([  0.161252,\
                            0.405198,\
                            0.512490,\
                            0.531526,\
                            0.719456])
    #This should just be RF from the gas itself, not including external forcing
    RF_compare = np.array([[2.005091,\
                            2.035587,\
                            2.097619,\
                            2.224813,\
                            2.483857]])
    alpha_compare = np.array([[ 0.125078,\
                                0.156356,\
                                0.178314,\
                                0.184508,\
                                0.179664]])
    assert np.allclose(emissions_out, emissions_compare)
    assert np.allclose(T_out, T_compare)
    assert np.allclose(RF_out, RF_compare)
    assert np.allclose(alpha_out, alpha_compare)

def test_dual_species():
    #emissions are based on exemplar CO2 & CH4 emissions
    out_dict = concentration_driven._run_numpy( inp_ar = np.array([[399.6173423,\
                                                                    401.8213271,\
                                                                    406.3410747,\
                                                                    415.7639971,\
                                                                    435.6164218],\
                                                                [   733.822081,\
                                                                    738.822081,\
                                                                    750.822081,\
                                                                    780.822081,\
                                                                    850.822081,\
                                                                ]]),\
                                            a1 = np.array([0.2173, 1]),\
                                            a2 = np.array([0.224, 0]),\
                                            a3 = np.array([0.2824, 0]),\
                                            a4 = np.array([0.2763, 0]),\
                                            tau1 = np.array([1000000, 9.150000]),\
                                            tau2 = np.array([394.4, 1]),\
                                            tau3 = np.array([36.54, 1]),\
                                            tau4 = np.array([4.304, 1]),\
                                            r0 = np.array([28.627296, 9.078874]),\
                                            rC = np.array([0.019773, 0]),\
                                            rT = np.array([4.334433, 0.000000]),\
                                            rA = np.array([0, -0.287247]),\
                                            PI_conc = np.array([278, 733.822081]),\
                                            emis2conc = np.array([0.468952343952344, 0.351714]),\
                                            f1 = np.array([5.754389, 0.061736]),\
                                            f2 = np.array([0.001215, -0.000049]),\
                                            f3 = np.array([-0.069598, 0.038416]),\
                                            d = np.array([  283,\
                                                            9.88,\
                                                            0.85]),\
                                            q = np.array([  0.311333,\
                                                            0.165417,\
                                                            0.242]),\
                                            ext_forcing = np.array([-0.256119925,\
                                                                    -0.304324144,\
                                                                    -0.501633962,\
                                                                    -0.904262779,\
                                                                    -0.12278275]),\
                                            timestep = np.array([1,\
                                                                2,\
                                                                4,\
                                                                8,\
                                                                1])
                                          )
    emissions_out = out_dict["emissions"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]
    #Given in the same order as inputs (i.e. CO2, CH4)
    emissions_compare = np.array([[ 319.828087,\
                            46.772580,\
                            23.786699,\
                            18.280039,\
                            18.201683],\
                          [ 7.506675,\
                            14.802304,\
                            34.585501,\
                            25.434119,\
                            25.054660]])
    T_compare = np.array([  0.161252,\
                            0.405665,\
                            0.514882,\
                            0.539532,\
                            0.740753])
    #This should just be RF from the gasses themselves, not including external forcing
    RF_compare = np.array([[2.005091,\
                            2.035587,\
                            2.097619,\
                            2.224813,\
                            2.483857],\
                           [0.000000,\
                            0.003714,\
                            0.012566,\
                            0.034339,\
                            0.083294]])
    alpha_compare = np.array([[ 0.125078,\
                                0.156356,\
                                0.178378,\
                                0.184777,\
                                0.180490],\
                              [ 0.992273,\
                                0.793780,\
                                0.371651,\
                                2.373146,\
                                206.035092]])
    assert np.allclose(emissions_out, emissions_compare)
    assert np.allclose(T_out, T_compare)
    assert np.allclose(RF_out, RF_compare)
    assert np.allclose(alpha_out, alpha_compare)

def test_constant_concentrations():
    
    out_dict = concentration_driven._run_numpy( inp_ar = np.array([[278,\
                                                                    278,\
                                                                    278,\
                                                                    278,\
                                                                    278],\
                                                                [ 733.822081,\
                                                                    733.822081,\
                                                                    733.822081,\
                                                                    733.822081,\
                                                                    733.822081,\
                                                                ]]),\
                                            a1 = np.array([0.2173, 1]),\
                                            a2 = np.array([0.224, 0]),\
                                            a3 = np.array([0.2824, 0]),\
                                            a4 = np.array([0.2763, 0]),\
                                            tau1 = np.array([1000000, 9.150000]),\
                                            tau2 = np.array([394.4, 1]),\
                                            tau3 = np.array([36.54, 1]),\
                                            tau4 = np.array([4.304, 1]),\
                                            r0 = np.array([28.627296, 9.078874]),\
                                            rC = np.array([0.019773, 0]),\
                                            rT = np.array([4.334433, 0.000000]),\
                                            rA = np.array([0, -0.287247]),\
                                            PI_conc = np.array([278, 733.822081]),\
                                            emis2conc = np.array([0.468952343952344, 0.351714]),\
                                            f1 = np.array([5.754389, 0.061736]),\
                                            f2 = np.array([0.001215, -0.000049]),\
                                            f3 = np.array([-0.069598, 0.038416]),\
                                            d = np.array([  283,\
                                                            9.88,\
                                                            0.85]),\
                                            q = np.array([  0.311333,\
                                                            0.165417,\
                                                            0.242]),\
                                            ext_forcing = np.array([0,\
                                                                    0,\
                                                                    0,\
                                                                    0,\
                                                                    0]),\
                                            timestep = np.array([1,\
                                                                2,\
                                                                4,\
                                                                8,\
                                                                1])
                                          )
    emissions_out = out_dict["emissions"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]
    #Given in the same order as inputs (i.e. CO2, CH4)
    emissions_compare = np.array([[ 0,\
                                    0,\
                                    0,\
                                    0,\
                                    0],\
                          [ 0.000000e+00,\
                            -1.800721e-13,\
                            -3.560150e-14,\
                            -3.560150e-14,\
                            -3.560150e-14]])
    T_compare = np.array([  0.0,\
                            0.0,\
                            0.0,\
                            0.0,\
                            0.0])
    #This should just be RF from the gasses themselves, not including external forcing
    RF_compare = np.array([[0.0,\
                            0.0,\
                            0.0,\
                            0.0,\
                            0.0],\
                           [0.0,\
                            0.0,\
                            0.0,\
                            0.0,\
                            0.0]])
    alpha_compare = np.array([[ 0.125078,\
                                0.125078,\
                                0.125078,\
                                0.125078,\
                                0.125078],\
                              [ 0.992273,\
                                0.992273,\
                                0.992273,\
                                0.992273,\
                                0.992273]])
    assert np.allclose(emissions_out, emissions_compare)
    assert np.allclose(T_out, T_compare)
    assert np.allclose(RF_out, RF_compare)
    assert np.allclose(alpha_out, alpha_compare)
