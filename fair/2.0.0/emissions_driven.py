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
        Dictionary containing the configuration for this run

    Returns
    -------
    :obj:`pd.DataFrame`
        Results of the run
    """
    raise NotImplementedError
    np_input = unifiedtools.convert_df_to_numpy(inp_df)
    res_dict = _run_numpy( np_input,\
                            a1 = cfg['a1'],\
                            a2 = cfg['a2'],\
                            a3 = cfg['a3'],\
                            a4 = cfg['a4'],\
                            tau1 = cfg['tau1'],\
                            tau2 = cfg['tau2'],\
                            tau3 = cfg['tau3'],\
                            tau4 = cfg['tau4'],\
                            r0 = cfg['r0'],\
                            rC = cfg['rC'],\
                            rT = cfg['rT'],\
                            rA = cfg['rA'],\
                            PI_conc = cfg['PI_conc'],\
                            emis2conc = cfg['emis2conc'],\
                            f1 = cfg['f1'],\
                            f2 = cfg['f2'],\
                            f3 = cfg['f3'],\
                            d = cfg['d'],\
                            q = cfg['q'],\
                            ext_forcing = cfg['ext_forcing'],\
                            timestep = cfg['timestep'])
    inp_df.columns.values.sort()
    inp_df.index.values.sort()
    C_array = res_dict["C"]
    RF_array = res_dict["RF"]
    T_array = res_dict["T"]
    alpha_array = res_dict["alpha"]
    C_df = unifiedtools.convert_numpy_output_to_df( C_array,\
                                                    inp_df.columns.values,\
                                                    inp_df.columns.name,\
                                                    inp_df.index.values,\
                                                    inp_df.index.name)
    RF_array = np.append(RF_array, np.array([cfg['ext_forcing']]))
    RF_array = np.append(RF_array, np.array([RF_array.sum(axis=0)]))
    RF_df = unifiedtools.convert_numpy_output_to_df(RF_array,\
                                                    np.append(inp_df.columns.values,np.array(['External Forcing', 'Total'])),\
                                                    inp_df.columns.name,\
                                                    inp_df.index.values,\
                                                    inp_df.index.name)
    T_df = unifiedtools.convert_numpy_output_to_df( np.array([T_array]),\
                                                    np.array(['T']),\
                                                    None,\
                                                    inp_df.index.values,\
                                                    inp_df.index.name)
    alpha_df = unifiedtools.convert_numpy_output_to_df( alpha_array,\
                                                        inp_df.columns.values,\
                                                        inp_df.columns.name,\
                                                        inp_df.index.values,\
                                                        inp_df.index.name)
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
    raise NotImplementedError
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
        alpha[...,i] = unifiedtools.calculate_alpha( G=G,\
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
        RF[...,i] = unifiedtools.step_forcing( C=C[...,i],\
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


