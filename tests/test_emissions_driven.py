import numpy as np
import pandas as pd
import pyam as pyam

import fair.emissions_driven as emissions_driven


def test_single_species():
    # emissions are based on exemplar CO2 emissions
    out_dict = emissions_driven._run_numpy(
        inp_ar=np.array(
            [
                [
                    3.00000000e00,
                    1.63829600e01,
                    6.36765930e01,
                    1.80148890e02,
                    2.43531850e01,
                ]
            ]
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
            [-0.256119925, -0.304324144, -0.501633962, -0.904262779, -0.12278275]
        ),
        timestep=np.array([1, 2, 4, 8, 8]),
    )
    C_out = out_dict["C"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]

    C_compare = np.array([[278.575556, 284.650637, 327.415460, 539.092185, 663.849687]])
    T_compare = np.array([-0.022563, -0.048115, 0.031111, 0.534905, 1.357223])
    # This should just be RF from the gas itself, not including external forcing
    RF_compare = np.array([[0.011400, 0.130324, 0.902589, 3.672639, 4.844845]])
    alpha_compare = np.array([[0.125078, 0.123069, 0.121295, 0.109470, 0.056788]])

    np.testing.assert_allclose(alpha_out, alpha_compare, atol=0.000001)
    np.testing.assert_allclose(C_out, C_compare, atol=0.000001)
    np.testing.assert_allclose(T_out, T_compare, atol=0.000001)
    np.testing.assert_allclose(RF_out, RF_compare, atol=0.000001)


def test_dual_species():
    # emissions are based on exemplar CO2 & CH4 emissions
    out_dict = emissions_driven._run_numpy(
        inp_ar=np.array(
            [
                [
                    3.00000000e00,
                    1.63829600e01,
                    6.36765930e01,
                    1.80148890e02,
                    2.43531850e01,
                ],
                [
                    2.1486753e02,
                    4.4001902e02,
                    1.70017052e03,
                    3.50032202e03,
                    2.3547352e02,
                ],
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
    C_out = out_dict["C"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]
    # Given in the same order as inputs (i.e. CO2, CH4)
    C_compare = np.array(
        [
            [278.575556, 284.652403, 327.534824, 543.411648, 672.682666],
            [769.601493, 959.837150, 2306.701551, 3499.109708, 3499.109708],
        ]
    )
    T_compare = np.array([-0.020142, -0.025803, 0.172729, 0.880953, 1.835630])
    # This should just be RF from the gasses themselves, not including external forcing
    RF_compare = np.array(
        [
            [0.011400, 0.130358, 0.904602, 3.717350, 4.919748],
            [0.026254, 0.155020, 0.798028, 1.192708, 1.192708],
        ]
    )
    alpha_compare = np.array(
        [
            [0.125078, 0.123296, 0.123140, 0.119969, 0.065278],
            [9.922729e-01, 8.111710e01, 7.697883e13, 2.245283e106, 2.245283e106],
        ]
    )
    np.testing.assert_allclose(C_out, C_compare, atol=0.000001)
    np.testing.assert_allclose(T_out, T_compare, atol=0.000001)
    np.testing.assert_allclose(RF_out, RF_compare, atol=0.000001)
    np.testing.assert_allclose(alpha_out, alpha_compare, atol=0.000001)


def test_zero_emissions():

    out_dict = emissions_driven._run_numpy(
        inp_ar=np.array([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0,]]),
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
    C_out = out_dict["C"]
    T_out = out_dict["T"]
    RF_out = out_dict["RF"]
    alpha_out = out_dict["alpha"]
    # Given in the same order as inputs (i.e. CO2, CH4)
    C_compare = np.array(
        [
            [278.0, 278.0, 278.0, 278.0, 278.0],
            [733.822081, 733.822081, 733.822081, 733.822081, 733.822081],
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
    np.testing.assert_allclose(C_out, C_compare, atol=0.000001)
    np.testing.assert_allclose(T_out, T_compare, atol=0.000001)
    np.testing.assert_allclose(RF_out, RF_compare, atol=0.000001)
    np.testing.assert_allclose(alpha_out, alpha_compare, atol=0.000001)


def test_run_df():
    year_index_np = np.array([2020, 2021, 2023, 2027, 2035])

    SIMPLE_DF = pd.DataFrame(
        [
            [
                "model_a",
                "scen_a",
                "World",
                "Emissions|CO2",
                "GtC/yr",
                3.00000000e00,
                1.63829600e01,
                6.36765930e01,
                1.80148890e02,
                2.43531850e01,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Emissions|CH4",
                "MtCH4/yr",
                2.14867530e02,
                4.40019020e02,
                1.70017052e03,
                3.50032202e03,
                2.35473520e02,
            ],
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
        data=gas_parameter_value_np, index=gas_parameter_name_np, columns=["CO2", "CH4"]
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

    res_df = emissions_driven.run(inp_df, cfg)

    SIMPLE_DF = pd.DataFrame(
        [
            [
                "model_a",
                "scen_a",
                "World",
                "Alpha|CO2",
                "None",
                0.125078,
                0.123296,
                0.123140,
                0.119969,
                0.065278,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Alpha|CH4",
                "None",
                9.922729e-01,
                8.111710e01,
                7.697883e13,
                2.245283e106,
                2.245283e106,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Atmospheric Concentrations|CO2",
                "ppm",
                278.575556,
                284.652403,
                327.534824,
                543.411648,
                672.682666,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Atmospheric Concentrations|CH4",
                "ppb",
                769.601493,
                959.837150,
                2306.701551,
                3499.109708,
                3499.109708,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Effective Radiative Forcing",
                "W/m**2",
                -0.21846593,
                -0.01894614,
                1.20099604,
                4.00579522,
                5.98967325,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Effective Radiative Forcing|CH4",
                "W/m**2",
                0.026254,
                0.155020,
                0.798028,
                1.192708,
                1.192708,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Effective Radiative Forcing|CO2",
                "W/m**2",
                0.011400,
                0.130358,
                0.904602,
                3.717350,
                4.919748,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Effective Radiative Forcing|Other",
                "W/m**2",
                -0.256119925,
                -0.304324144,
                -0.501633962,
                -0.904262779,
                -0.12278275,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Emissions|CO2",
                "GtC/yr",
                3.00000000e00,
                1.63829600e01,
                6.36765930e01,
                1.80148890e02,
                2.43531850e01,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Emissions|CH4",
                "MtCH4/yr",
                2.14867530e02,
                4.40019020e02,
                1.70017052e03,
                3.50032202e03,
                2.35473520e02,
            ],
            [
                "model_a",
                "scen_a",
                "World",
                "Surface Temperature",
                "K",
                -0.020142,
                -0.025803,
                0.172729,
                0.880953,
                1.835630,
            ],
        ],
        columns=pyam.IAMC_IDX + [2020, 2021, 2023, 2027, 2035],
    )
    compare_df = pyam.IamDataFrame(SIMPLE_DF)

    pd.testing.assert_frame_equal(
        res_df.timeseries(), compare_df.timeseries(), check_exact=False, atol=1e-04
    )
