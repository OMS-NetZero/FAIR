import numpy as np

from .tools import unifiedtools


def run(inp_df, cfg):
    """
    Run FaIR 2.0

    Parameters
    ----------
    inp_df : :obj:`pd.DataFrame`
        Input :obj:`pd.DataFrame` containing the timeseries to run,
        in IAMC compliant DataFrame format
        (i.e. A multiIndex of Model, Region, Scenario,
        Unit, Variable then Columns for time)

    cfg : dict
        Dictionary containing the configuration for this run,
        in format {'gas_params' : :obj:`pd.DataFrame`,
        'thermal_params': :obj:`pd.DataFrame`,
        'ext_forcing' : :obj:`pd.DataFrame`}

    Returns
    -------
    :obj:`pd.DataFrame`
        Results of the run
    """
    arg_list, forcing_list = unifiedtools.return_np_function_arg_list(
        inp_df, cfg, concentration_mode=False
    )

    res_dict = _run_numpy(*arg_list)

    res_df_iamc_compliant = unifiedtools.create_output_dataframe_iamc_compliant(
        inp_df,
        res_dict["C"],
        res_dict["RF"],
        res_dict["T"],
        res_dict["alpha"],
        forcing_list,
    )

    return res_df_iamc_compliant


def _run_numpy(
    inp_ar,
    a1,
    a2,
    a3,
    a4,
    tau1,
    tau2,
    tau3,
    tau4,
    r0,
    rC,
    rT,
    rA,
    PI_conc,
    emis2conc,
    f1,
    f2,
    f3,
    d,
    q,
    ext_forcing,
    timestep,
    mapping_ar,
):
    """
    Run FaIR 2.0 from numpy array

    This function can *only* run one scenario,
    thermal parameter set & gas parameter set at a time

    Parameters
    ----------
    inp_ar : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing the timeseries to run. No checks
        of the column order are performed here.
        format: [[species],[time]]

    a1, a2, a3, a4, tau1, tau2, tau3, tau4, r0, rC,
    rT, rA, PI_conc, emis2conc : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing gas parameters in format:
        [species],
        note: all species contain the same number of gas/thermal pool
        indices (some are simply populated with 0)

    f1, f2, f3 : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing gas forcing parameters
        in format: [forcing]

    d, q : obj:`np.ndarray`
        Input :obj:`np.ndarray` containing thermal parameters in format:
        [response box]

    ext_forcing : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing any other prescribed forcing
        in format: [time]

    timestep : :obj:`np.ndarray`
        Input :obj:`np.ndarray`
        specifying the length of each entry in inp_ar in years.
        For example: if inp_ar were an nx4 array,
        representing times 2020-2021, 2021-2023, 2023-2027 and 2027-2028:
        timestep would be: np.array([1,2,4,1])

    mapping_ar : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing mapping between gases and forcing
        for example: [0,1] would just map gas -> forcing directly
        [0,0,1] would map for two forcings from the gas at index 0


    Returns
    -------
    dict
        Dictionary containing the results of the run.
        Keys are 'C', 'RF', 'T', 'alpha' and 'S'
        (Concentration, Radiative Forcing, Temperature, Alpha and Temperature Boxes)
        Values are in :obj:`np.ndarray` format,
        with the final index representing 'timestep'
    """

    n_species, n_timesteps = inp_ar.shape
    # Concentrations and Alpha
    C, alpha = np.zeros((2, n_species, n_timesteps))

    n_forcing = len(f1)
    # Radiative Forcing
    RF = np.zeros((n_forcing, n_timesteps))
    # Temperature
    T = np.zeros(n_timesteps)
    # S represents the results of the calculations from the thermal boxes,
    # an Impulse Response calculation (T = sum(S))
    S = np.zeros_like(d)
    # G represents cumulative emissions,
    # while G_A represents emissions accumulated since pre-industrial times,
    # both in the same units as emissions
    # So at any point, G - G_A is equal
    # to the amount of a species that has been absorbed
    G_A, G = np.zeros((2, n_species))
    # R in format [[index],[species]]
    R = np.zeros((4, n_species))
    # a,tau in format [[index], [species]]
    a = np.array([a1, a2, a3, a4])
    tau = np.array([tau1, tau2, tau3, tau4])
    # g0, g1 in format [species]
    g0, g1 = unifiedtools.calculate_g(a=a, tau=tau)
    for i, tstep in enumerate(timestep):
        alpha[..., i] = unifiedtools.calculate_alpha(
            G=G, G_A=G_A, T=np.sum(S, axis=0), r0=r0, rC=rC, rT=rT, rA=rA, g0=g0, g1=g1
        )
        C[..., i], R, G_A = unifiedtools.step_concentration(
            emissions=inp_ar[np.newaxis, ..., i],
            a=a,
            dt=tstep,
            alpha=alpha[np.newaxis, ..., i],
            tau=tau,
            R_old=R,
            G_A_old=G_A,
            PI_conc=PI_conc,
            emis2conc=emis2conc,
        )
        RF[..., i] = unifiedtools.step_forcing(
            C=C[mapping_ar, i], PI_conc=PI_conc[mapping_ar], f1=f1, f2=f2, f3=f3
        )
        S, T[i] = unifiedtools.step_temperature(
            S_old=S, F=np.sum(RF[..., i], axis=0) + ext_forcing[i], q=q, d=d, dt=tstep
        )
        G += inp_ar[..., i]
    res = {"C": C, "RF": RF, "T": T, "alpha": alpha, "S": S}
    return res
