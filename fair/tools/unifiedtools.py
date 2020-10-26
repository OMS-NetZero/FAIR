import numexpr as ne
import numpy as np
import pandas as pd
import pyam as pyam

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
        "emissions * a / decay_rate * ( 1. - decay_factor ) + R_old * decay_factor"
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


def convert_df_to_numpy(inp_df):
    """
    Convert input df to numpy array with correct order for running. Note, this method deliberately sorts based on lower case values, this is an attempt to make the software easier to use (i.e.
    one might absent mindedly write a gas as 'CO2' in the gas parameter dataframe and 'co2' in the concentrations/emissions dataframe)

    Parameters
    ----------

    inp_df : :obj:`pd.DataFrame`
        Input :obj:`pd.DataFrame` to be converted to a numpy array. Column represent e.g. gas species, index represents e.g. time
    Returns
    -------
    :obj:`np.ndarray`
        df converted to a numpy array with correct order for running. The numpy array is ordered as follows: [Column, Index]
        The columns and indexes are returned in a sorted order

    """

    # sorts the df so ordering is 'correct' within levels/index
    # sort the df columns
    df = inp_df.iloc[:, inp_df.columns.astype("str").str.lower().argsort()]
    # Sort the df index
    df = df.iloc[df.index.astype("str").str.lower().argsort()]
    # converts df to a numpy.ndarray [Column, Index]
    res = df.to_numpy().T
    return res


def convert_numpy_output_to_df(
    res_numpy, column_labels, column_name, index_labels, index_name
):
    """
    Convert the numpy output to a dataframe i.e. add metadata to the outputs. Note that this function assumes that the labels have been given in the correct order.

    Parameters
    ----------

    res_numpy : :obj:`np.ndarray`
        Input :obj:`np.ndarray` to be converted to :obj:`pd.DataFrame`. Array should be n1xn2, where n1 represents columns and n2 represents index.
        e.g. array([[1,2,3,4,5],[6,7,8,9,10]]) would represent 2 gasses with 5 timesteps

    column_labels : :obj:`np.ndarray`
        Input :obj:`np.ndarray` giving labels for the columns (e.g. typically gasses), running with the above example:
        e.g. array(["CO2",["CH4]]) would label the gasses

    column_name : string
        Input string used to label the columns, in our example this might be "Gas"

    index_labels : :obj:`np.ndarray`
        Input :obj:`np.ndarray` giving labels for the index,
        e.g. array([2020,2021,2022,2023,2024]) would label our example above from 2020 to 2024

    index_name : : string
        Input string used to label the index, in our example this might be "Year"

    Returns
    -------
    :obj:`np.ndarray`
        df converted to a numpy array with correct order for running. The numpy array is ordered as follows: [Column, Index]
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
        "emissions * a / decay_rate * ( 1. - decay_factor ) + R_old * decay_factor"
    )

    return emissions, R_new

def create_output_dataframe_aimc_compliant(inp_df, gas_np, RF_np, T_np, alpha_np, ext_forcing_np):
    """
    TODO: docstring
    """

    units = Units()

    inp_df_ts = inp_df.timeseries()
    variable_index_level = inp_df_ts.index.names.tolist().index("variable")
    variable_array = inp_df_ts.index.levels[variable_index_level].to_numpy()
    function_array, gas_array = (np.array([variable_array[i].split('|',1) for i in range(len(variable_array))])).T

    unique_function_array = np.unique(function_array)
    if len(unique_function_array) > 1:
        raise Exception('Error: More than one type of input passed')
    output_function = "Atmospheric Concentrations" if str.lower(unique_function_array[0][0]) == 'e' else "Emissions"

    model_region_scenario_array = np.unique(inp_df.timeseries().index.droplevel(('unit','variable')).to_numpy())
    if len(model_region_scenario_array) > 1:
        raise ValueError('More than one Model, Region + Scenario combination input passed')
    model_region_scenario = np.array(model_region_scenario_array[0])

    data_array = np.append(
                            np.append(     
                                        np.repeat(  
                                                    model_region_scenario[...,np.newaxis],3*gas_np.shape[0],axis=1).T,\
                                        np.append(  np.array(   
                                                                [[var, units[var]] for var in\
                                                                    [j + '|' + k for j,k in \
                                                                        zip([output_function, "Effective Radiative Forcing", "Alpha"]*gas_np.shape[0],np.repeat([gas_array],3))]]),\
                                                    np.array(
                                                                [*gas_np.T,*RF_np.T,*alpha_np.T]).T.reshape(3*gas_np.shape[0],gas_np.shape[1])\
                                                    ,axis=1)\
                                        ,axis=1),
                            [   [*model_region_scenario,"Effective Radiative Forcing|External Forcing", units["Effective Radiative Forcing|External Forcing"],*ext_forcing_np],\
                                [*model_region_scenario,"Effective Radiative Forcing", units["Effective Radiative Forcing"],*(RF_np.sum(axis=0) + ext_forcing_np)],\
                                [*model_region_scenario,"Surface Temperature", units["Surface Temperature"],*T_np]], axis = 0)

    output_df = pyam.IamDataFrame(pd.DataFrame(data_array, columns = pyam.IAMC_IDX + inp_df.timeseries().columns.to_list(),)).append(inp_df)

    return output_df


def return_np_function_arg_list(inp_df, cfg, concentration_mode=False):
    """
    TODO: docstring
    """

    gas_params_numpy = convert_df_to_numpy(cfg["gas_params"]).T

    time_index = inp_df.timeseries().columns
    arg_list = [
        inp_df.timeseries().to_numpy() + ((gas_params_numpy[9] * gas_params_numpy[4])[..., np.newaxis] if concentration_mode else 0),
        *gas_params_numpy[[0,1,2,3,14,15,16,17,10,12,13,11,9,5,6,7,8]],
        *convert_df_to_numpy(cfg["thermal_params"]).T,
        *convert_df_to_numpy(cfg["ext_forcing"]),
        np.append(np.diff(time_index), np.diff(time_index)[-1]).astype('timedelta64[Y]').astype('float64'),
    ]
    return arg_list
