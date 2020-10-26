import numpy as np
import pandas as pd
import pyam as pyam

import fair.concentration_driven as concentration_driven
from fair.tools import unifiedtools


def test_single_species():
    # emissions are based on example CO2 concentrations
    out_dict = concentration_driven._run_numpy(
        inp_ar=np.array(
            [[399.6173423, 401.8213271, 406.3410747, 415.7639971, 435.6164218]]
        ),
        a1=np.array([0.2173]),
        a2=np.array([0.224]),
        a3=np.array([0.2824]),
        a4=np.array([0.2763]),
        tau1=np.array([1000000]),
        tau2=np.array([394.4]),
        tau3=np.array([36.54]),
        tau4=np.array([4.304]),
        r0=np.array([28.627296]),
        rC=np.array([0.019773]),
        rT=np.array([4.334433]),
        rA=np.array([0]),
        PI_conc=np.array([278]),
        emis2conc=np.array([0.468952343952344]),
        f1=np.array([5.754389]),
        f2=np.array([0.001215]),
        f3=np.array([-0.069598]),
        d=np.array([283, 9.88, 0.85]),
        q=np.array([0.311333, 0.165417, 0.242]),
        ext_forcing=np.array(
            [-0.256119925, -0.304324144, -0.501633962,
             -0.904262779, -0.12278275]
        ),
        timestep=np.array([1, 2, 4, 8, 8]),
    )
    emissions_out = out_dict["emissions"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]

    emissions_compare = np.array(
        [[319.828087, 46.772580, 23.792089, 18.293743, 18.238578]]
    )
    T_compare = np.array([0.161252, 0.405198, 0.512490,
                          0.531526, 0.719456])
    # This is RF from the gas itself, not including external forcing
    RF_compare = np.array([[2.005091, 2.035587, 2.097619, 2.224813, 2.483857]])
    alpha_compare = np.array([[0.125078, 0.156356, 0.178314, 0.184508, 0.179664]])
    np.testing.assert_allclose(alpha_out, alpha_compare, atol=0.000001)
    np.testing.assert_allclose(emissions_out, emissions_compare, atol=0.000001)
    np.testing.assert_allclose(T_out, T_compare, atol=0.000001)
    np.testing.assert_allclose(RF_out, RF_compare, atol=0.000001)


def test_dual_species():
    # emissions are based on exemplar CO2 & CH4 emissions
    out_dict = concentration_driven._run_numpy(
        inp_ar=np.array(
            [
                [399.6173423, 401.8213271, 406.3410747, 415.7639971, 435.6164218],
                [733.822081, 738.822081, 750.822081, 780.822081, 850.822081,],
            ]
        ),
        a1=np.array([0.2173, 1]),
        a2=np.array([0.224, 0]),
        a3=np.array([0.2824, 0]),
        a4=np.array([0.2763, 0]),
        tau1=np.array([1000000, 9.150000]),
        tau2=np.array([394.4, 1]),
        tau3=np.array([36.54, 1]),
        tau4=np.array([4.304, 1]),
        r0=np.array([28.627296, 9.078874]),
        rC=np.array([0.019773, 0]),
        rT=np.array([4.334433, 0.000000]),
        rA=np.array([0, -0.287247]),
        PI_conc=np.array([278, 733.822081]),
        emis2conc=np.array([0.468952343952344, 0.351714]),
        f1=np.array([5.754389, 0.061736]),
        f2=np.array([0.001215, -0.000049]),
        f3=np.array([-0.069598, 0.038416]),
        d=np.array([283, 9.88, 0.85]),
        q=np.array([0.311333, 0.165417, 0.242]),
        ext_forcing=np.array(
            [-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]
        ),
        timestep=np.array([1, 2, 4, 8, 8]),
    )
    emissions_out = out_dict["emissions"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]
    # Given in the same order as inputs (i.e. CO2, CH4)
    emissions_compare = np.array(
        [
            [319.828087, 46.772580, 23.786699, 18.280039, 18.201683],
            [7.506675, 14.802304, 34.585501, 25.434119, 25.054660],
        ]
    )
    T_compare = np.array([0.161252, 0.405665, 0.514882, 0.539532, 0.740753])
    # This should just be RF from the gasses themselves, not including external forcing
    RF_compare = np.array(
        [
            [2.005091, 2.035587, 2.097619, 2.224813, 2.483857],
            [0.000000, 0.003714, 0.012566, 0.034339, 0.083294],
        ]
    )
    alpha_compare = np.array(
        [
            [0.125078, 0.156356, 0.178378, 0.184777, 0.180490],
            [0.992273, 0.793780, 0.371651, 2.373146, 206.035092],
        ]
    )
    np.testing.assert_allclose(alpha_out, alpha_compare, atol=0.000001)
    np.testing.assert_allclose(emissions_out, emissions_compare, atol=0.000001)
    np.testing.assert_allclose(T_out, T_compare, atol=0.000001)
    np.testing.assert_allclose(RF_out, RF_compare, atol=0.000001)


def test_constant_concentrations():

    out_dict = concentration_driven._run_numpy(
        inp_ar=np.array(
            [
                [278, 278, 278, 278, 278],
                [733.822081, 733.822081, 733.822081, 733.822081, 733.822081,],
            ]
        ),
        a1=np.array([0.2173, 1]),
        a2=np.array([0.224, 0]),
        a3=np.array([0.2824, 0]),
        a4=np.array([0.2763, 0]),
        tau1=np.array([1000000, 9.150000]),
        tau2=np.array([394.4, 1]),
        tau3=np.array([36.54, 1]),
        tau4=np.array([4.304, 1]),
        r0=np.array([28.627296, 9.078874]),
        rC=np.array([0.019773, 0]),
        rT=np.array([4.334433, 0.000000]),
        rA=np.array([0, -0.287247]),
        PI_conc=np.array([278, 733.822081]),
        emis2conc=np.array([0.468952343952344, 0.351714]),
        f1=np.array([5.754389, 0.061736]),
        f2=np.array([0.001215, -0.000049]),
        f3=np.array([-0.069598, 0.038416]),
        d=np.array([283, 9.88, 0.85]),
        q=np.array([0.311333, 0.165417, 0.242]),
        ext_forcing=np.array([0, 0, 0, 0, 0]),
        timestep=np.array([1, 2, 4, 8, 8]),
    )
    emissions_out = out_dict["emissions"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]
    # Given in the same order as inputs (i.e. CO2, CH4)
    emissions_compare = np.array(
        [
            [0, 0, 0, 0, 0],
            [0.000000e00, -1.800721e-13, -3.560150e-14, -3.560150e-14, -3.560150e-14],
        ]
    )
    T_compare = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    # This should just be RF from the gasses themselves, not including external forcing
    RF_compare = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]])
    alpha_compare = np.array(
        [
            [0.125078, 0.125078, 0.125078, 0.125078, 0.125078],
            [0.992273, 0.992273, 0.992273, 0.992273, 0.992273],
        ]
    )
    np.testing.assert_allclose(alpha_out, alpha_compare, atol=0.000001)
    np.testing.assert_allclose(emissions_out, emissions_compare, atol=0.000001)
    np.testing.assert_allclose(T_out, T_compare, atol=0.000001)
    np.testing.assert_allclose(RF_out, RF_compare, atol=0.000001)


def test_run_df():
    year_index_np = np.array([2020, 2021, 2023, 2027, 2035])

    SIMPLE_DF = pd.DataFrame(   [
                                ['model_a', 'scen_a', 'World', 'Atmospheric Concentrations|CO2', 'ppm', 399.6173423, 401.8213271, 406.3410747, 415.7639971, 435.6164218],
                                ['model_a', 'scen_a', 'World', 'Atmospheric Concentrations|CH4', 'ppb', 733.822081 , 738.822081 , 750.822081 , 780.822081 , 850.822081],
                                ],
                                columns=pyam.IAMC_IDX + [2020, 2021, 2023, 2027, 2035],
                            )

    inp_df = pyam.IamDataFrame(SIMPLE_DF)

    gas_parameter_value_np = np.array(
        [
            [0.2173, 1],
            [0.224, 0],
            [0.2824, 0],
            [0.2763, 0],
            [1000000, 9.150000],
            [394.4, 1],
            [36.54, 1],
            [4.304, 1],
            [28.627296, 9.078874],
            [0.019773, 0],
            [4.334433, 0.000000],
            [0, -0.287247],
            [278, 733.822081],
            [0.468952343952344, 0.351714],
            [5.754389, 0.061736],
            [0.001215, -0.000049],
            [-0.069598, 0.038416],
            [False, False],
        ]
    )
    gas_parameter_name_np = np.array(
        [
            "a1",
            "a2",
            "a3",
            "a4",
            "tau1",
            "tau2",
            "tau3",
            "tau4",
            "r0",
            "rC",
            "rT",
            "rA",
            "PI_conc",
            "emis2conc",
            "f1",
            "f2",
            "f3",
            "aer_conc",
        ]
    )

    gas_params_df = pd.DataFrame(
        data=gas_parameter_value_np, index=gas_parameter_name_np, columns=['CO2','CH4']
    )

    thermal_parameter_value_np = np.array(
        [[283, 9.88, 0.85], [0.311333, 0.165417, 0.242]]
    )

    thermal_parameter_name_np = ["d", "q"]

    thermal_params_df = pd.DataFrame(
        data=thermal_parameter_value_np,
        index=thermal_parameter_name_np,
        columns=[1, 2, 3],
    )

    ext_forcing_value_np = np.array(
        [-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]
    )

    ext_forcing_df = pd.DataFrame(
        data=ext_forcing_value_np, index=year_index_np, columns=["External Forcing"]
    )

    cfg = {
        "gas_params": gas_params_df,
        "thermal_params": thermal_params_df,
        "ext_forcing": ext_forcing_df,
    }

    res_df = concentration_driven.run(inp_df, cfg)

    SIMPLE_DF = pd.DataFrame([['model_a', 'scen_a', 'World', 'Alpha|CO2', 'None', 1.25078000e-01, 1.56356000e-01, 1.78378000e-01, 1.84777000e-01,
        1.80490000e-01],
       ['model_a', 'scen_a', 'World', 'Alpha|CH4', 'None', 9.92273000e-01, 7.93780000e-01, 3.71651000e-01, 2.37314600e+00,
        2.06035092e+02],
       ['model_a', 'scen_a', 'World', 'Atmospheric Concentrations|CO2', 'ppm', 399.6173423, 401.8213271, 406.3410747, 415.7639971, 435.6164218],
    ['model_a', 'scen_a', 'World', 'Atmospheric Concentrations|CH4', 'ppb', 733.822081 , 738.822081 , 750.822081 , 780.822081 , 850.822081],
       ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing', 'W/m**2', 1.74897108, 1.73497686, 1.60855104, 1.35488922, 2.44436825],
       ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing|CH4', 'W/m**2', 0.000000, 0.003714, 0.012566, 0.034339, 0.083294],
       ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing|CO2', 'W/m**2', 2.005091, 2.035587, 2.097619, 2.224813, 2.483857],
       ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing|External Forcing', 'W/m**2', -0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275],
       ['model_a', 'scen_a', 'World', 'Emissions|CO2', 'GtC/yr', 319.828087, 46.772580, 23.786699, 18.280039, 18.201683],
       ['model_a', 'scen_a', 'World', 'Emissions|CH4', 'MtCH4/yr', 7.506675, 14.802304, 34.585501, 25.434119, 25.054660],
       ['model_a', 'scen_a', 'World', 'Surface Temperature', 'K',0.161252, 0.405665, 0.514882, 0.539532, 0.740753]],
    columns=pyam.IAMC_IDX + [2020, 2021, 2023, 2027, 2035],
)   
    compare_df = pyam.IamDataFrame(SIMPLE_DF)

    pd.testing.assert_frame_equal(res_df.timeseries(), compare_df.timeseries(), check_exact = False, atol = 1e-04)
