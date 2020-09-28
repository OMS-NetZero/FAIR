import pandas as pd
import numpy as np
from tools import unifiedtools

def run(inp_df, cfg):
    """
    Run FaIR 2.0

    Parameters
    ----------
    inp_df : :obj:`pd.DataFrame`
        Input :obj:`pd.DataFrame` containing the timeseries to run

    cfg : dict
        Dictionary containing the configuration for this run, in format {'gas_params' : :obj:`pd.DataFrame`, 'thermal_params': :obj:`pd.DataFrame`, 'ext_forcing' : :obj:`pd.DataFrame`}

    Returns
    -------
    :obj:`pd.DataFrame`
        Results of the run
    """
    raise NotImplementedError
    emissions_np = unifiedtools.convert_df_to_numpy(inp_df)

    #a1_np, ..., tau4_np are in format [gas]
    a1_np, a2_np, a3_np, a4_np, aer_conc_np, emis2conc_np, f1_np, f2_np, f3_np, PI_conc_np, r0_np, rA_np, rC_np, rT_np, tau1_np, tau2_np, tau3_np, tau4_np = unifiedtools.convert_df_to_numpy(cfg['gas_params'])

    #d_np, q_np are in format [thermal box]
    d_np, q_np = unifiedtools.convert_df_to_numpy(cfg['thermal_params'])

    ext_forcing_np = unifiedtools.convert_df_to_numpy(cfg['ext_forcing'])

    inp_df.columns.values.sort()
    inp_df.index.values.sort()
    time_index = inp_df.index

    timestep_np = np.append(np.diff(time_index),np.diff(time_index)[-1])

    res_dict = _run_numpy(  emissions_np,\
                            a1 = a1_np,\
                            a2 = a2_np,\
                            a3 = a3_np,\
                            a4 = a4_np,\
                            tau1 = tau1_np,\
                            tau2 = tau2_np,\
                            tau3 = tau3_np,\
                            tau4 = tau4_np,\
                            r0 = r0_np,\
                            rC = rC_np,\
                            rT = rT_np,\
                            rA = rA_np,\
                            PI_conc = PI_conc_np,\
                            emis2conc = emis2conc_np,\
                            f1 = f1_np,\
                            f2 = f2_np,\
                            f3 = f3_np,\
                            d = d_np,\
                            q = q_np,\
                            ext_forcing = ext_forcing_np,\
                            timestep = timestep_np)
    C_df, RF_df, T_df, alpha_df = unifiedtools.create_output_dataframes(inp_df,\
                                                                        res_dict["C"],\
                                                                        res_dict["RF"],\
                                                                        res_dict["T"],\
                                                                        res_dict["alpha"],\
                                                                        ext_forcing_np)
    res_df_dict = {'emissions':inp_df, 'C':C_df, 'RF' : RF_df, 'T' : T_df, 'alpha':alpha_df}
    return res_df_dict

def _run_numpy( inp_ar,\
                a1,\
                a2,\
                a3,\
                a4,\
                tau1,\
                tau2,\
                tau3,\
                tau4,\
                r0,\
                rC,\
                rT,\
                rA,\
                PI_conc,\
                emis2conc,\
                f1,\
                f2,\
                f3,\
                d,\
                q,\
                ext_forcing,\
                timestep):
    """
    Run FaIR 2.0 from numpy array

    This function can *only* run one scenario, thermal parameter set & gas parameter set at a time

    Parameters
    ----------
    inp_ar : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing the timeseries to run. No checks
        of the column order are performed here.
        format: [[species],[time]]

    a1, a2, a3, a4, tau1, tau2, tau3, tau4, r0, rC, rT, rA, PI_conc, emis2conc, f1, f2, f3 : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing gas parameters in format: [species], 
        note: all species contain the same number of gas/thermal pool indices (some are simply populated with 0)

    d, q : obj:`np.ndarray`
        Input :obj:`np.ndarray` containing thermal parameters in format: [response box]
    
    ext_forcing : :obj:`np.ndarray`
        Input :obj:`np.ndarray` containing any other prescribed forcing in format: [time]
    
    timestep : :obj:`np.ndarray`
        Input :obj:`np.ndarray` specifying the length of each entry in inp_ar in years.
        For example: if inp_ar were an nx4 array, representing times 2020-2021, 2021-2023, 2023-2027 and 2027-2028:
        timestep would be: np.array([1,2,4,1])
    

    Returns
    -------
    dict
        Dictionary containing the results of the run. 
        Keys are 'C', 'RF', 'T', and 'alpha'
        (Concentration, Radiative Forcing, Temperature and Alpha)
        Values are in :obj:`np.ndarray` format, with the final index representing 'timestep'
    """

    n_species, n_timesteps = inp_ar.shape
    #Concentration, Radiative Forcing and Alpha
    C, RF, alpha = np.zeros((3,n_species,n_timesteps))
    #Temperature
    T = np.zeros(n_timesteps)
    #S represents the results of the calculations from the thermal boxes, an Impulse Response calculation (T = sum(S))
    S = np.zeros_like(d)
    #G represents cumulative emissions, while G_A represents emissions accumulated since pre-industrial times, both in the same units as emissions
    #So at any point, G - G_A is equal to the amount of a species that has been absorbed
    G_A, G = np.zeros((2,n_species))
    #R in format [[index],[species]]
    R = np.zeros((4,n_species))
    #a,tau in format [[index], [species]]
    a = np.array([a1,a2,a3,a4])
    tau = np.array([tau1,tau2,tau3,tau4])
    #g0, g1 in format [species]
    g0,g1 = unifiedtools.calculate_g(a = a,tau = tau)
    for i, tstep in enumerate(timestep):
        alpha[...,i] = unifiedtools.calculate_alpha(G=G,\
                                                    G_A=G_A,\
                                                    T=T[i-1],\
                                                    r0=r0,\
                                                    rC=rC,\
                                                    rT=rT,\
                                                    rA=rA,\
                                                    g0=g0,\
                                                    g1=g1)
        C[...,i], R, G_A = unifiedtools.step_concentration( emissions = inp_ar[np.newaxis,...,i],\
                                                            a = a,\
                                                            dt = tstep,\
                                                            alpha = alpha[np.newaxis,...,i],\
                                                            tau = tau,\
                                                            R_old = R,\
                                                            G_A_old = G_A,\
                                                            PI_conc = PI_conc,\
                                                            emis2conc = emis2conc)
        RF[...,i] = unifiedtools.step_forcing(  C=C[...,i],\
                                                PI_conc=PI_conc,\
                                                f1=f1,\
                                                f2=f2,\
                                                f3=f3)
        S,T[i] = unifiedtools.step_temperature( S_old=S,\
                                                F=np.sum(RF[...,i],axis=0) + ext_forcing[i],\
                                                q=q,\
                                                d=d,\
                                                dt=tstep)
        G += inp_ar[i]
    res = {'C':C, 'RF' : RF, 'T' : T, 'alpha':alpha}
    return res


