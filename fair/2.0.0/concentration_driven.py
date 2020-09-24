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
    res_numpy = _run_numpy( np_input,\
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
    res = unifiedtools.convert_numpy_output_to_df(res_numpy)
    return res

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
        of the column order are performed here. This array contains gasses as concentrations (NOT Aerosol Emissions), i.e. PI_conc + cumulative emissions
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
    #Emissions, Radiative Forcing and Alpha
    emissions, RF, alpha = np.zeros((3,n_species,n_timesteps))
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

    G_A_list = np.append((inp_ar[...,:-1] + inp_ar[...,1:])/2,(3*inp_ar[...,-1] - inp_ar[...,-2])/2)
    G_A_list = (G_A-PI_conc)/emis2conc

    for i, tstep in enumerate(timestep):
        G_A = G_A_list[i]
        alpha[...,i] = unifiedtools.calculate_alpha(G=G,\
                                                    G_A=G_A,\
                                                    T=T,\
                                                    r0=r0,\
                                                    rC=rC,\
                                                    rT=rT,\
                                                    rA=rA,\
                                                    g0=g0,\
                                                    g1=g1)
        emissions[...,i], R = unifiedtools.unstep_concentration(a = a,\
                                                                dt = tstep,\
                                                                alpha = alpha[np.newaxis,...,i],\
                                                                tau = tau,\
                                                                R_old = R,\
                                                                G_A = G_A)
        G +=emissions[...,i]
        S,T[i] = unifiedtools.step_temperature( S_old=S,\
                                                F=np.sum(RF[...,i],axis=0) + ext_forcing[i],\
                                                q=q,\
                                                d=d,\
                                                dt=tstep)

    res = {'emissions' : emissions, 'RF' : RF, 'T' : T, 'alpha':alpha}
    return res
