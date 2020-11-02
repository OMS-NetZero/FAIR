import datetime as dt

import numpy as np
import pandas as pd
import pyam as pyam

from fair.tools import unifiedtools


def test_calculate_alpha():
    # Based on CO2
    G_1 = 619
    G_A_1 = 259
    T_1 = 1.24068
    r0_1 = 28.63
    rC_1 = 0.019773
    rT_1 = 4.334433
    rA_1 = 0
    g0_1 = 0.020369508
    g1_1 = 11.4137078
    iirf100_max = 15

    # Based on CH4
    G_2 = 34000
    G_A_2 = 3100
    T_2 = 1.24068
    r0_2 = 8.444641
    rC_2 = 0
    rT_2 = -0.287247
    rA_2 = 0.000343
    g0_2 = 0.850699166
    g1_2 = 9.148042797

    alpha_no_max_out = unifiedtools.calculate_alpha(
        np.array([G_1, G_2]),
        np.array([G_A_1, G_A_2]),
        np.array([T_1, T_2]),
        np.array([r0_1, r0_2]),
        np.array([rC_1, rC_2]),
        np.array([rT_1, rT_2]),
        np.array([rA_1, rA_2]),
        np.array([g0_1, g0_2]),
        np.array([g1_1, g1_2]),
    )
    alpha_with_max_out = unifiedtools.calculate_alpha(
        np.array([G_1, G_2]),
        np.array([G_A_1, G_A_2]),
        np.array([T_1, T_2]),
        np.array([r0_1, r0_2]),
        np.array([rC_1, rC_2]),
        np.array([rT_1, rT_2]),
        np.array([rA_1, rA_2]),
        np.array([g0_1, g0_2]),
        np.array([g1_1, g1_2]),
        iirf100_max,
    )

    alpha_no_max_compare = np.array([0.74788084, 2.31332918])
    alpha_with_max_compare = np.array([0.07581137, 2.31332918])

    np.testing.assert_allclose(alpha_no_max_out, alpha_no_max_compare)
    np.testing.assert_allclose(alpha_with_max_out, alpha_with_max_compare)


def test_calculate_g():
    a1 = np.array([0.2173, 1])
    a2 = np.array([0.224, 0])
    a3 = np.array([0.2824, 0])
    a4 = np.array([0.2763, 0])
    tau1 = np.array([1000000, 9.15])
    tau2 = np.array([394.4, 1])
    tau3 = np.array([36.54, 1])
    tau4 = np.array([4.304, 1])
    a = np.array([a1, a2, a3, a4])
    tau = np.array([tau1, tau2, tau3, tau4])

    g0_out, g1_out = unifiedtools.calculate_g(a, tau)

    g0_compare = np.array([0.010183697764797717, 0.3678073392252985])
    g1_compare = np.array([11.4137078, 9.1480428])
    np.testing.assert_allclose(g1_out, g1_compare)
    np.testing.assert_allclose(g0_out, g0_compare)


def test_step_concentration():

    emissions = np.array([10.2123, 361])
    a1 = np.array([0.2173, 1])
    a2 = np.array([0.224, 0])
    a3 = np.array([0.2824, 0])
    a4 = np.array([0.2763, 0])
    tau1 = np.array([1000000, 9.15])
    tau2 = np.array([394.4, 1])
    tau3 = np.array([36.54, 1])
    tau4 = np.array([4.304, 1])
    a = np.array([a1, a2, a3, a4])
    tau = np.array([tau1, tau2, tau3, tau4])
    alpha = np.array([0.74788084, 2.31332918])
    dt = 10
    R_1 = np.array([63.10423371, 1091.54903])
    R_2 = np.array([42.23272521, 0])
    R_3 = np.array([14.23992555, 0])
    R_4 = np.array([2.040457847, 0])
    R_old = np.array([R_1, R_2, R_3, R_4])
    G_A_old = np.array([259.3, 3104])
    PI_conc = np.array([278, 720])
    emis2conc = np.array([0.468952344, 0.351714258])

    C_out, R_out, G_A_out = unifiedtools.step_concentration(
        emissions, a, dt, alpha, tau, R_old, G_A_old, PI_conc, emis2conc
    )

    C_compare = np.array([383.68002648, 1891.49581172])
    R_compare = np.array(
        [
            [85.29456948, 3557.63389782],
            [63.31706123, 0.0],
            [34.02781502, 0.0],
            [8.767445190, 0.0],
        ]
    )
    G_A_compare = np.array([191.40689091, 3557.63389782])

    np.testing.assert_allclose(C_out, C_compare)
    np.testing.assert_allclose(R_out, R_compare)
    np.testing.assert_allclose(G_A_out, G_A_compare)


def test_step_forcing():
    C = np.array([399.6, 1811])
    PI_conc = np.array([278, 720])
    f1 = np.array([5.754389, 0.061736])
    f2 = np.array([0.001215, -0.000049])
    f3 = np.array([-0.069598, 0.038416])

    RF_out = unifiedtools.step_forcing(C, PI_conc, f1, f2, f3)
    RF_compare = np.array([2.0048501, 0.60750117])

    np.testing.assert_allclose(RF_out, RF_compare)


def test_step_temperature():
    d = np.array([283, 9.88, 0.85])
    q = np.array([0.311333, 0.165417, 0.242])
    F = np.array([2.870])
    S_old = np.array([0.173196459, 1.06749276, 2])
    dt = 10

    S_out, T_out = unifiedtools.step_temperature(S_old, F, q, d, dt)

    S_compare = np.array([0.19820533, 0.69017337, 0.69455015])
    T_compare = np.array([2.41180904])

    np.testing.assert_allclose(S_out, S_compare)
    np.testing.assert_allclose(T_out, T_compare)


def test_unstep_concentration():
    a1 = np.array([0.2173, 1])
    a2 = np.array([0.224, 0])
    a3 = np.array([0.2824, 0])
    a4 = np.array([0.2763, 0])
    tau1 = np.array([1000000, 9.15])
    tau2 = np.array([394.4, 1])
    tau3 = np.array([36.54, 1])
    tau4 = np.array([4.304, 1])
    a = np.array([a1, a2, a3, a4])
    tau = np.array([tau1, tau2, tau3, tau4])
    alpha = np.array([0.74788084, 2.31332918])
    dt = 10
    R_1 = np.array([63.10423371, 1091.54903])
    R_2 = np.array([42.23272521, 0])
    R_3 = np.array([14.23992555, 0])
    R_4 = np.array([2.040457847, 0])
    R_old = np.array([R_1, R_2, R_3, R_4])
    G_A = np.array([259.3, 3104])

    emissions_out, R_out = unifiedtools.unstep_concentration(
        a=a, dt=dt, alpha=alpha[np.newaxis, ...], tau=tau, R_old=R_old, G_A=G_A
    )

    emissions_compare = np.array([19.15739779, 304.0803832])
    R_compare = np.array(
        [
            [104.73213702, 3104.0],
            [83.01823459, 0.0],
            [55.18263368, 0.0],
            [16.36699471, 0.0],
        ]
    )

    np.testing.assert_allclose(emissions_out, emissions_compare)
    np.testing.assert_allclose(R_out, R_compare)


def test_convert_df_to_numpy():
    test_df = pd.DataFrame(
        data={"a": [7, 4, 10, 1], "c": [9, 6, 12, 3], "b": [8, 5, 11, 2]},
        index=[3, 2, 4, 1],
    )
    res_compare = np.array([[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]])
    res = unifiedtools.convert_df_to_numpy(test_df)
    np.testing.assert_allclose(res_compare, res)


def test_convert_numpy_output_to_df():
    res = np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]])
    column_labels = np.array(["a", "b", "c"])
    column_name = "Gas"
    index_labels = np.array([2019, 2020, 2022, 2023, 2050])
    index_name = "Year"
    df = unifiedtools.convert_numpy_output_to_df(
        res, column_labels, column_name, index_labels, index_name
    )
    df_compare = pd.DataFrame(
        data={"a": [1, 2, 3, 4, 5], "b": [6, 7, 8, 9, 10], "c": [11, 12, 13, 14, 15]},
        index=[2019, 2020, 2022, 2023, 2050],
    )

    assert df.equals(df_compare)
    assert not np.sum(df.index.tolist() != np.array([2019, 2020, 2022, 2023, 2050]))
    assert df.index.name == "Year"
    assert not np.sum(df.columns.tolist() != np.array(["a", "b", "c"]))
    assert df.columns.name == "Gas"

def test_return_np_function_arg_list():
    year_index_np = np.array([2020, 2021, 2023, 2027, 2035])

    SIMPLE_DF = pd.DataFrame(   [
                                ['model_a', 'scen_a', 'World', 'Emissions|CO2', 'GtC/yr', 3.00000000e+00, 1.63829600e+01, 6.36765930e+01, 1.80148890e+02,
        2.43531850e+01],
                                ['model_a', 'scen_a', 'World', 'Emissions|CH4', 'MtCH4/yr', 2.14867530e+02, 4.40019020e+02, 1.70017052e+03, 3.50032202e+03,
        2.35473520e+02],
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
            [True, False],
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

    emission_mode_arg_list = unifiedtools.return_np_function_arg_list(
        inp_df, cfg, concentration_mode=False
    )
    concentration_mode_arg_list = unifiedtools.return_np_function_arg_list(
        inp_df, cfg, concentration_mode=True
    )

    emission_mode_arg_list_compare = [
        np.array(
            [
                [
                    2.14867530e02,
                    4.40019020e02,
                    1.70017052e03,
                    3.50032202e03,
                    2.35473520e02,
                ],
                [
                    3.00000000e00,
                    1.63829600e01,
                    6.36765930e01,
                    1.80148890e02,
                    2.43531850e01,
                ],
            ]
        ),
        np.array([1.00000000e00, 2.17300000e-01]),
        np.array([0.00000000e00, 2.24000000e-01]),
        np.array([0.00000000e00, 2.82400000e-01]),
        np.array([0.00000000e00, 2.76300000e-01]),
        np.array([9.15000000e00, 1.00000000e06]),
        np.array([1.00000000e00, 3.94400000e02]),
        np.array([1.00000000e00, 3.65400000e01]),
        np.array([1.00000000e00, 4.30400000e00]),
        np.array([9.07887400e00, 2.86272960e01]),
        np.array([0.00000000e00, 1.97730000e-02]),
        np.array([0.00000000e00, 4.33443300e00]),
        np.array([-2.87247000e-01, 0.00000000e00]),
        np.array([7.33822081e02, 2.78000000e02]),
        np.array([3.51714000e-01, 4.68952344e-01]),
        np.array([6.17360000e-02, 5.75438900e00]),
        np.array([-4.90000000e-05, 1.21500000e-03]),
        np.array([3.84160000e-02, -6.95980000e-02]),
        np.array([283, 9.88, 0.85]),
        np.array([0.311333, 0.165417, 0.242]),
        np.array([-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]),
        np.array([1, 2, 4, 8, 8]),
    ]
    concentration_mode_arg_list_compare = [
        np.array(
            [
                [
                    2.14867530e02,
                    4.40019020e02,
                    1.70017052e03,
                    3.50032202e03,
                    2.35473520e02,
                ],
                [281.0, 294.38296, 341.676593, 458.14889, 302.353185],
            ]
        ),
        np.array([1.00000000e00, 2.17300000e-01]),
        np.array([0.00000000e00, 2.24000000e-01]),
        np.array([0.00000000e00, 2.82400000e-01]),
        np.array([0.00000000e00, 2.76300000e-01]),
        np.array([9.15000000e00, 1.00000000e06]),
        np.array([1.00000000e00, 3.94400000e02]),
        np.array([1.00000000e00, 3.65400000e01]),
        np.array([1.00000000e00, 4.30400000e00]),
        np.array([9.07887400e00, 2.86272960e01]),
        np.array([0.00000000e00, 1.97730000e-02]),
        np.array([0.00000000e00, 4.33443300e00]),
        np.array([-2.87247000e-01, 0.00000000e00]),
        np.array([7.33822081e02, 2.78000000e02]),
        np.array([3.51714000e-01, 4.68952344e-01]),
        np.array([6.17360000e-02, 5.75438900e00]),
        np.array([-4.90000000e-05, 1.21500000e-03]),
        np.array([3.84160000e-02, -6.95980000e-02]),
        np.array([283, 9.88, 0.85]),
        np.array([0.311333, 0.165417, 0.242]),
        np.array([-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]),
        np.array([1, 2, 4, 8, 8]),
    ]
    for i in range(len(emission_mode_arg_list)):
        np.testing.assert_allclose(
            emission_mode_arg_list[i], emission_mode_arg_list_compare[i]
        )
        np.testing.assert_allclose(
            concentration_mode_arg_list[i], concentration_mode_arg_list_compare[i]
        )

def test_return_np_function_arg_list_no_params():
    year_index_np = np.array([2020, 2021, 2023, 2027, 2035])

    SIMPLE_DF = pd.DataFrame(   [
                                ['model_a', 'scen_a', 'World', 'Emissions|CO2', 'GtC/yr', 3.00000000e+00, 1.63829600e+01, 6.36765930e+01, 1.80148890e+02,
        2.43531850e+01],
                                ['model_a', 'scen_a', 'World', 'Emissions|CH4', 'MtCH4/yr', 2.14867530e+02, 4.40019020e+02, 1.70017052e+03, 3.50032202e+03,
        2.35473520e+02],
                                ],
                                columns=pyam.IAMC_IDX + [2020, 2021, 2023, 2027, 2035],
                            )

    inp_df = pyam.IamDataFrame(SIMPLE_DF)

    ext_forcing_value_np = np.array(
        [-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]
    )

    ext_forcing_df = pd.DataFrame(
        data=ext_forcing_value_np, index=year_index_np, columns=["External Forcing"]
    )

    cfg = {
        "ext_forcing": ext_forcing_df,
    }

    emission_mode_arg_list = unifiedtools.return_np_function_arg_list(
        inp_df, cfg, concentration_mode=False
    )
    concentration_mode_arg_list = unifiedtools.return_np_function_arg_list(
        inp_df, cfg, concentration_mode=True
    )

    emission_mode_arg_list_compare = [
        np.array(
            [
                [
                    2.14867530e02,
                    4.40019020e02,
                    1.70017052e03,
                    3.50032202e03,
                    2.35473520e02,
                ],
                [
                    3.00000000e00,
                    1.63829600e01,
                    6.36765930e01,
                    1.80148890e02,
                    2.43531850e01,
                ],
            ]
        ),
        np.array([1.00000000e00, 2.17300000e-01]),
        np.array([0.00000000e00, 2.24000000e-01]),
        np.array([0.00000000e00, 2.82400000e-01]),
        np.array([0.00000000e00, 2.76300000e-01]),
        np.array([8.8000000e00, 1.00000000e09]),
        np.array([1.00000000e00, 3.94400000e02]),
        np.array([1.00000000e00, 3.65400000e01]),
        np.array([1.00000000e00, 4.30400000e00]),
        np.array([8.8, 30.4]),
        np.array([0.    , 0.0177]),
        np.array([-0.33,  2.64]),
        np.array([0.00032, 0.     ]),
        np.array([720., 278.]),
        np.array([0.351665695, 0.468887594]),
        np.array([0.  , 4.57]),
        np.array([0., 0.]),
        np.array([0.0385, 0.086 ]),
        np.array([283, 9.88, 0.85]),
        np.array([0.311333, 0.165417, 0.242]),
        np.array([-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]),
        np.array([1, 2, 4, 8, 8]),
    ]
    concentration_mode_arg_list_compare = [
        np.array(
            [
                [
                    2.14867530e02,
                    4.40019020e02,
                    1.70017052e03,
                    3.50032202e03,
                    2.35473520e02,
                ],
                [
                    3.00000000e00,
                    1.63829600e01,
                    6.36765930e01,
                    1.80148890e02,
                    2.43531850e01,
                ],
            ]
        ),
        np.array([1.00000000e00, 2.17300000e-01]),
        np.array([0.00000000e00, 2.24000000e-01]),
        np.array([0.00000000e00, 2.82400000e-01]),
        np.array([0.00000000e00, 2.76300000e-01]),
        np.array([8.8000000e00, 1.00000000e09]),
        np.array([1.00000000e00, 3.94400000e02]),
        np.array([1.00000000e00, 3.65400000e01]),
        np.array([1.00000000e00, 4.30400000e00]),
        np.array([8.8, 30.4]),
        np.array([0.    , 0.0177]),
        np.array([-0.33,  2.64]),
        np.array([0.00032, 0.     ]),
        np.array([720., 278.]),
        np.array([0.351665695, 0.468887594]),
        np.array([0.  , 4.57]),
        np.array([0., 0.]),
        np.array([0.0385, 0.086 ]),
        np.array([283, 9.88, 0.85]),
        np.array([0.311333, 0.165417, 0.242]),
        np.array([-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]),
        np.array([1, 2, 4, 8, 8]),
    ]
    for i in range(len(emission_mode_arg_list)):
        np.testing.assert_allclose(
            emission_mode_arg_list[i], emission_mode_arg_list_compare[i]
        )
        np.testing.assert_allclose(
            concentration_mode_arg_list[i], concentration_mode_arg_list_compare[i]
        )

def test_create_output_dataframe_iamc_compliant():
    SIMPLE_DF = pd.DataFrame(   [
                                ['model_a', 'scen_a', 'World', 'Emissions|CO2', 'GtC/yr', 3.00000000e+00, 1.63829600e+01, 6.36765930e+01, 1.80148890e+02],
                                ['model_a', 'scen_a', 'World', 'Emissions|CH4', 'MtCH4/yr', 2.14867530e+02, 4.40019020e+02, 1.70017052e+03, 3.50032202e+03],
                                ],
                                columns=pyam.IAMC_IDX + [dt.date(year = 2020, month = 1, day = 1), dt.date(year = 2021, month = 1, day = 1), dt.date(year = 2023, month = 1, day = 1), dt.date(year = 2027, month = 1, day = 1)],
                            )
    
    inp_df = pyam.IamDataFrame(SIMPLE_DF)

    gas_np = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
    RF_np = np.array([[9, 10, 11, 12], [13, 14, 15, 16]])
    T_np = np.array([17, 18, 19, 20])
    alpha_np = np.array([[21, 22, 23, 24], [25, 26, 27, 28]])
    ext_forcing_np = np.array([0, 0, 0, 0])

    out_df = unifiedtools.create_output_dataframe_iamc_compliant(inp_df, gas_np, RF_np, T_np, alpha_np, ext_forcing_np)

    SIMPLE_DF = pd.DataFrame(   [
                                ['model_a', 'scen_a', 'World', 'Alpha|CO2', 'None', 25, 26, 27, 28],
                                ['model_a', 'scen_a', 'World', 'Alpha|CH4', 'None', 21, 22, 23, 24],
                                ['model_a', 'scen_a', 'World', 'Atmospheric Concentrations|CO2', 'ppm', 5, 6, 7, 8],
                                ['model_a', 'scen_a', 'World', 'Atmospheric Concentrations|CH4', 'ppb', 1, 2, 3, 4],
                                ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing', 'W/m**2', 22, 24, 26, 28],
                                ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing|CH4', 'W/m**2', 9, 10, 11, 12],
                                ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing|CO2', 'W/m**2', 13, 14, 15, 16],
                                ['model_a', 'scen_a', 'World', 'Effective Radiative Forcing|Other', 'W/m**2', 0, 0, 0, 0],
                                ['model_a', 'scen_a', 'World', 'Emissions|CO2', 'GtC/yr', 3.00000000e+00, 1.63829600e+01, 6.36765930e+01, 1.80148890e+02],
                                ['model_a', 'scen_a', 'World', 'Emissions|CH4', 'MtCH4/yr', 2.14867530e+02, 4.40019020e+02, 1.70017052e+03, 3.50032202e+03],
                                ['model_a', 'scen_a', 'World', 'Surface Temperature', 'K', 17, 18, 19, 20]
                                ],
                                columns=pyam.IAMC_IDX + [dt.date(year = 2020, month = 1, day = 1), dt.date(year = 2021, month = 1, day = 1), dt.date(year = 2023, month = 1, day = 1), dt.date(year = 2027, month = 1, day = 1)],
                            )

    compare_df = pyam.IamDataFrame(SIMPLE_DF)

    pd.testing.assert_frame_equal(out_df.timeseries(), compare_df.timeseries())



