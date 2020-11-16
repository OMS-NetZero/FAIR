import numexpr as ne
import numpy as np
import pandas as pd
import pyam as pyam
import re

from ..ancil.default_gas_parameters import get_gas_params
from ..ancil.default_thermal_parameters import get_thermal_params
from ..ancil.units import Units


def calculate_alpha(G, G_A, T, r0, rC, rT, rA, g0, g1, iirf100_max=False):
    """
    TODO: docstring
    """
    iirf100_val = ne.evaluate("abs(r0 + rC * (G-G_A) + rT * T + rA * G_A)")
    if iirf100_max:
        iirf100_val = ne.evaluate(  # noqa: F841
            "where(iirf100_val>iirf100_max,iirf100_max,iirf100_val)"
        )
    alpha_val = ne.evaluate("g0 * exp(iirf100_val / g1)")

    return alpha_val


def calculate_g(a, tau):
    """
    TODO: docstring
    """
    g1 = ne.evaluate(
        "sum( a * tau * ( 1. - ( 1. + 100/tau ) * exp(-100/tau) ), axis = 0)"
    )
    g0 = np.exp(-1 * np.sum(a * tau * (1.0 - np.exp(-100 / tau)), axis=0) / g1)
    return g0, g1


def step_concentration(
    emissions, a, dt, alpha, tau, R_old, G_A_old, PI_conc, emis2conc
):
    """
    TODO: docstring
    """
    decay_rate = ne.evaluate("1/(alpha*tau)")  # noqa: F841
    decay_factor = ne.evaluate("exp(-dt*decay_rate)")  # noqa: F841
    R = ne.evaluate(
        "emissions * a / decay_rate *\
            ( 1. - decay_factor ) + R_old * decay_factor"
    )
    G_A = ne.evaluate("sum(R,axis=0)") 
    C = ne.evaluate("PI_conc + emis2conc * (G_A + G_A_old) / 2")

    return C, R, G_A


def step_forcing(C, PI_conc, f1, f2, f3):
    """
    TODO: docstring
    """
    logforc = ne.evaluate(
        "f1 * where( (C/PI_conc) <= 0, 0, log(C/PI_conc) )",
        {"f1": f1, "C": C, "PI_conc": PI_conc},
    )
    linforc = ne.evaluate("f2 * (C - PI_conc)", {"f2": f2, "C": C, "PI_conc": PI_conc})
    sqrtforc = ne.evaluate(
        "f3 * ( (sqrt( where(C<0 ,0 ,C ) ) - sqrt(PI_conc)) )",
        {"f3": f3, "C": C, "PI_conc": PI_conc},
    )

    RF = logforc + linforc + sqrtforc

    return RF


def step_temperature(S_old, F, q, d, dt=1):
    """
    TODO: docstring
    """
    decay_factor = ne.evaluate("exp(-dt/d)")  # noqa: F841
    S_new = ne.evaluate("q * F * (1 - decay_factor) + S_old * decay_factor")
    T = ne.evaluate("sum( (S_old + S_new)/2, axis=0 )")

    return S_new, T

def sort_df(inp_df):
    """
    Convert input df to output df sorted alphabetically.
    Note, this method deliberately sorts based on lower case values,
    this is an attempt to make the software easier to use
    (i.e. one might absent mindedly write a gas as 'CO2' in the gas parameter
    dataframe and 'co2' in the concentrations/emissions dataframe)

    Parameters
    ----------

    inp_df : :obj:`pd.DataFrame`
        Input :obj:`pd.DataFrame` to be sorted by rows/columns.
        Column represent e.g. gas species, index represents e.g. time

    Returns
    -------

    :obj:`pd.DataFrame`
        df sorted by rows/columns.

    """
    # sorts the df so ordering is 'correct' within levels/index
    # sort the df columns
    df = inp_df.iloc[:, inp_df.columns.astype("str").str.lower().argsort()]
    # Sort the df index
    df = df.iloc[df.index.astype("str").str.lower().argsort()]

    return df



def convert_df_to_numpy(inp_df):
    """
    Convert input df to numpy array with correct order for running.
    Note, this method deliberately sorts based on lower case values,
    this is an attempt to make the software easier to use
    (i.e. one might absent mindedly write a gas as 'CO2' in the gas parameter
    dataframe and 'co2' in the concentrations/emissions dataframe)

    Parameters
    ----------

    inp_df : :obj:`pd.DataFrame`
        Input :obj:`pd.DataFrame` to be converted to a numpy array.
        Column represent e.g. gas species, index represents e.g. time

    Returns
    -------

    :obj:`np.ndarray`
        df converted to a numpy array with correct order for running.
        The numpy array is ordered as follows: [Column, Index]
        The columns and indexes are returned in a sorted order

    """
    # converts df to a numpy.ndarray [Column, Index]
    res = sort_df(inp_df).to_numpy().T
    return res


def convert_numpy_output_to_df(
    res_numpy, column_labels, column_name, index_labels, index_name
):
    """
    Convert the numpy output to a dataframe i.e. add metadata to the outputs.
    Note that this function assumes that the labels
    have been given in the correct order.

    Parameters
    ----------

    res_numpy : :obj:`np.ndarray`
        Input :obj:`np.ndarray` to be converted to :obj:`pd.DataFrame`.
        Array should be n1xn2,
        where n1 represents columns and n2 represents index.
        e.g. array([[1,2,3,4,5],[6,7,8,9,10]])
        would represent 2 gasses with 5 timesteps

    column_labels : :obj:`np.ndarray`
        Input :obj:`np.ndarray` giving labels for the columns
        (e.g. typically gasses), running with the above example:
        e.g. array(["CO2",["CH4]]) would label the gasses

    column_name : string
        Input string used to label the columns,
        in our example this might be "Gas"

    index_labels : :obj:`np.ndarray`
        Input :obj:`np.ndarray` giving labels for the index,
        e.g. array([2020,2021,2022,2023,2024]) would label our example above
        from 2020 to 2024

    index_name : : string
        Input string used to label the index,
        in our example this might be "Year"

    Returns
    -------
    :obj:`np.ndarray`
        df converted to a numpy array with correct order for running.
        The numpy array is ordered as follows: [Column, Index]
        The columns and indexes are returned in a sorted order

    """
    # here you add metadata so users know timesteps, units etc. that were used
    # in the run
    res = pd.DataFrame(data=res_numpy.T)
    res.columns = column_labels
    res.columns.name = column_name
    res.index = index_labels
    res.index.name = index_name
    return res


def unstep_concentration(a, dt, alpha, tau, R_old, G_A):
    """
    TODO: docstring
    """

    decay_rate = ne.evaluate("1/(alpha*tau)")
    decay_factor = ne.evaluate("exp(-dt*decay_rate)")
    emissions = (G_A - np.sum(R_old * decay_factor, axis=0)) / np.sum(
        a / decay_rate * (1.0 - decay_factor), axis=0
    )
    R_new = ne.evaluate(
        "emissions * a / decay_rate *\
            ( 1. - decay_factor ) + R_old * decay_factor"
    )

    return emissions, R_new

def split_vars(inp_df_ts, i):
    return [var.split("|", 1)[i] for 
        var in inp_df_ts.index.get_level_values('variable')]

def create_output_dataframe_iamc_compliant(
    inp_df, gas_np, RF_np, T_np, alpha_np, forcing_list
):
    """
    TODO: docstring
    """
    inp_df_ts = inp_df.timeseries()
    ext_forc_np = np.zeros(RF_np.shape[1])
    fn_np = np.unique(split_vars(inp_df_ts,0))

    if len(fn_np) == 2 and 'Effective Radiative Forcing' in fn_np:
        ext_forc_np = np.sum(
            inp_df_ts.filter(regex = 'Effective Radiative Forcing', axis = 0)
            .to_numpy(),axis=0
            )
    elif len(fn_np) > 1:
        raise Exception("Error: Too many inputs passed")
    
    [[out_fn,gases]] = \
        [
            [
                f[f[0] in fn_np]+"|",
                split_vars(
                    inp_df_ts.filter(
                        regex=f[f[1] in fn_np],
                        axis=0,
                        )
                    ,1,
                    ),
            ]
        for f in 
        [
            [
                "Emissions",
                "Atmospheric Concentrations",
            ]
        ]
        ]

    model_region_scenario_array = np.unique(
        inp_df_ts.index.droplevel(("unit", "variable")).to_numpy())
    if len(model_region_scenario_array) > 1:
        raise ValueError( "More than one Model, \
            Region + Scenario combination input passed")
    units = Units()
    header = lambda var: [*model_region_scenario_array[0], var, units[var]]
    
    data_array = np.array(
        [
            *[
                [
                    *header("Effective Radiative Forcing|"+forc), *RF_np[i]
                ]
                    for i,forc in enumerate(forcing_list)
            ],
            *[
                [
                    *header(out_fn+gas), *gas_np[i]
                ] 
                    for i,gas in enumerate(gases)
            ],
            *[
                [
                    *header("Alpha|"+gas), *alpha_np[i]
                ]
                    for i,gas in enumerate(gases)
            ],
            [
                *header("Effective Radiative Forcing"),
                *(RF_np.sum(axis=0) + ext_forc_np)
            ],
            [
                *header("Surface Temperature"),*T_np
            ],
        ]
    )

    output_df = pyam.IamDataFrame(
        pd.DataFrame(
            data_array, columns=pyam.IAMC_IDX + inp_df_ts.columns.to_list()
        )
    ).append(inp_df)

    return output_df

def return_np_function_arg_list(inp_df, cfg, concentration_mode=False):
    """
    TODO: docstring
    """
    inp_df_ts = inp_df.timeseries()
    unique_fn_array = np.unique(split_vars(inp_df_ts,0))
    
    if len(unique_fn_array) >2 or (len(unique_fn_array) ==2
    and not 'Effective Radiative Forcing' in unique_fn_array):
        raise Exception("Error: Too many inputs passed")

    gas_df = inp_df_ts.filter(
        regex = [
            'Emissions',
            'Atmospheric Concentrations',][concentration_mode], axis = 0
    )
    
    gas_np = np.array(split_vars(gas_df,1))
    forc_p_df = get_gas_params(gas_np)
    if "gas_params" in cfg:
        forc_p_df[cfg["gas_params"].columns.tolist()] = cfg["gas_params"]
    forc_p_df = sort_df(forc_p_df)

    forc_np = np.array(forc_p_df.columns.tolist())

    forc_p_np = convert_df_to_numpy(forc_p_df).T
    gas_p_np = forc_p_np\
        [
            :,
            [
                np.where(forc_np == g)[0][0] for g in gas_np #matches forcings to gases
            ],
        ][
            [
                0,#a1
                1,#a2
                2,#a3
                3,#a4
                14,#tau1
                15,#tau2
                16,#tau3
                17,#tau4
                10,#r0
                12,#rC
                13,#rT
                11,#rA
                9,#PI_conc
                5,#emis2conc
                4,#aer_conc
            ]
        ]
    
    arg_list = \
    [
        gas_df.to_numpy() + 
        np.multiply(*gas_p_np
            [
                [
                    12,#PI_conc
                    14,#aer_conc
                ]
            ]
        )
            [...,np.newaxis]*concentration_mode
        , #emissions/concentrations, with concentration mode aerosol concentrations converted to PI_conc + aer_conc
        *gas_p_np
        [
            0:14 #a1,a2,a3,a4,tau1,tau2,tau3,tau4,r0,rC,rT,rA,PI_conc,emis2conc
        ]
        , #gas cycle parameters
        *forc_p_np
        [
            [
                6,#f1
                7,#f2
                8,#f3
            ]
        ]
        , #forcing parameters
        *convert_df_to_numpy(
            cfg["thermal_params"]
            if "thermal_params" in cfg #thermal params inputted
            else get_thermal_params() #default thermal params
        ).T
        , #thermal parameters, either inputted or default
        np.sum(
            inp_df_ts.filter(
                regex = 'Effective Radiative Forcing', axis = 0
            ).to_numpy(),
            axis=0
        )
        , #external forcing, if 'Effective Radiative Forcing' is not present, this returns an array of zeros
        np.diff(gas_df.columns)[
            [
                *range(
                    len(gas_df.columns)-1
                ),
                -1,
            ]
        ]
        .astype("timedelta64[Y]").astype("float64")
        , #this gives the time difference between different gas emission/conc timestamps (assumes the final timestep is the same as the penultimate)
        np.array(
            [
                i
                for forc in forc_np
                for i,gas in enumerate(gas_np)
                if re.match(r'^'+gas,forc)
            ],
        )
        , #maps forcings onto gases, i.e. BC|BC on Snow -> BC
    ]


    return arg_list, forc_np