import numexpr as ne
import numpy as np
import pandas as pd


def calculate_alpha(G, G_A, T, r0, rC, rT, rA, g0, g1, iirf100_max=False):
    iirf100_val = ne.evaluate("abs(r0 + rC * (G-G_A) + rT * T + rA * G_A)")
    if iirf100_max:
        iirf100_val = ne.evaluate(  # noqa: F841
            "where(iirf100_val>iirf100_max,iirf100_max,iirf100_val)"
        )
    alpha_val = ne.evaluate("g0 * exp(iirf100_val / g1)")

    return alpha_val


def calculate_g(a, tau):
    g1 = ne.evaluate(
        "sum( a * tau * ( 1. - ( 1. + 100/tau ) * exp(-100/tau) ), axis = 0)"
    )
    g0 = np.exp(-1 * np.sum(a * tau * (1.0 - np.exp(-100 / tau)), axis=0) / g1)
    return g0, g1


def step_concentration(
    emissions, a, dt, alpha, tau, R_old, G_A_old, PI_conc, emis2conc
):
    decay_rate = ne.evaluate("1/(alpha*tau)")  # noqa: F841
    decay_factor = ne.evaluate("exp(-dt*decay_rate)")  # noqa: F841
    R = ne.evaluate(
        "emissions * a / decay_rate * ( 1. - decay_factor ) + R_old * decay_factor"
    )
    G_A = ne.evaluate("sum(R,axis=0)")
    C = ne.evaluate("PI_conc + emis2conc * (G_A + G_A_old) / 2")

    return C, R, G_A


def step_forcing(C, PI_conc, f1, f2, f3):
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
    :obj:'np.ndarray'
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
    :obj:'np.ndarray'
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

    decay_rate = ne.evaluate("1/(alpha*tau)")
    decay_factor = ne.evaluate("exp(-dt*decay_rate)")
    emissions = (G_A - np.sum(R_old * decay_factor, axis=0)) / np.sum(
        a / decay_rate * (1.0 - decay_factor), axis=0
    )
    R_new = ne.evaluate(
        "emissions * a / decay_rate * ( 1. - decay_factor ) + R_old * decay_factor"
    )

    return emissions, R_new


def create_output_dataframes(inp_df, gas_np, RF_np, T_np, alpha_np, ext_forcing_np):
    gas_df = convert_numpy_output_to_df(
        gas_np,
        inp_df.columns.values,
        inp_df.columns.name,
        inp_df.index.values,
        inp_df.index.name,
    )
    RF_np = np.append(RF_np, ext_forcing_np[np.newaxis, ...], axis=0)
    RF_np = np.append(RF_np, RF_np.sum(axis=0)[np.newaxis, ...], axis=0)
    RF_df = convert_numpy_output_to_df(
        RF_np,
        np.append(inp_df.columns.values, np.array(["External Forcing", "Total"])),
        inp_df.columns.name,
        inp_df.index.values,
        inp_df.index.name,
    )
    T_df = convert_numpy_output_to_df(
        np.array([T_np]), np.array(["T"]), None, inp_df.index.values, inp_df.index.name
    )
    alpha_df = convert_numpy_output_to_df(
        alpha_np,
        inp_df.columns.values,
        inp_df.columns.name,
        inp_df.index.values,
        inp_df.index.name,
    )
    return gas_df, RF_df, T_df, alpha_df


def return_np_function_arg_list(inp_df, cfg, concentration_mode=False):

    # a1_np, ..., tau4_np are in format [gas]
    (
        a1_np,
        a2_np,
        a3_np,
        a4_np,
        aer_conc_np,
        emis2conc_np,
        f1_np,
        f2_np,
        f3_np,
        PI_conc_np,
        r0_np,
        rA_np,
        rC_np,
        rT_np,
        tau1_np,
        tau2_np,
        tau3_np,
        tau4_np,
    ) = convert_df_to_numpy(cfg["gas_params"]).T

    # d_np, q_np are in format [thermal box]
    d_np, q_np = convert_df_to_numpy(cfg["thermal_params"]).T

    ext_forcing_np = convert_df_to_numpy(cfg["ext_forcing"])[0]

    inp_df = inp_df.iloc[inp_df.index.astype("str").str.lower().argsort()]
    time_index = inp_df.index

    timestep_np = np.append(np.diff(time_index), np.diff(time_index)[-1])

    # Turns aerosol emissions into concentrations, by adding in Pre-Industrial Concentration
    inp_np = convert_df_to_numpy(inp_df) + (
        (PI_conc_np * aer_conc_np)[..., np.newaxis] if concentration_mode else 0
    )

    arg_list = [
        inp_np,
        a1_np,
        a2_np,
        a3_np,
        a4_np,
        tau1_np,
        tau2_np,
        tau3_np,
        tau4_np,
        r0_np,
        rC_np,
        rT_np,
        rA_np,
        PI_conc_np,
        emis2conc_np,
        f1_np,
        f2_np,
        f3_np,
        d_np,
        q_np,
        ext_forcing_np,
        timestep_np,
    ]
    return arg_list
