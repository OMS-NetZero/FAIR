import pytest
import numpy as np
import pandas as pd

from tools import unifiedtools

def calculate_alpha_test():
    #Based on CO2
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

    #Based on CH4
    G_2 = 34000
    G_A_2 = 3100
    T_2 = 1.24068
    r0_2 = 8.444641
    rC_2 = 0
    rT_2 = -0.287247
    rA_2 = 0.000343
    g0_2 = 0.850699166
    g1_2 = 9.148042797

    alpha_no_max_out = unifiedtools.calculate_alpha( np.array([G_1, G_2]),\
                                                     np.array([G_A_1, G_A_2]),\
                                                     np.array([T_1, T_2]),\
                                                     np.array([r0_1, r0_2]),\
                                                     np.array([rC_1, rC_2]),\
                                                     np.array([rT_1, rT_2]),\
                                                     np.array([rA_1, rA_2]),\
                                                     np.array([g0_1, g0_2]),\
                                                     np.array([g1_1, g1_2]))
    alpha_with_max_out = unifiedtools.calculate_alpha(  np.array([G_1, G_2]),\
                                                        np.array([G_A_1, G_A_2]),\
                                                        np.array([T_1, T_2]),\
                                                        np.array([r0_1, r0_2]),\
                                                        np.array([rC_1, rC_2]),\
                                                        np.array([rT_1, rT_2]),\
                                                        np.array([rA_1, rA_2]),\
                                                        np.array([g0_1, g0_2]),\
                                                        np.array([g1_1, g1_2]),\
                                                        iirf100_max)

    alpha_no_max_compare = np.array([0.74788084, 2.31332918])
    alpha_with_max_compare = np.array([0.07581137, 2.31332918])

    assert np.allclose(alpha_no_max_out, alpha_no_max_compare)
    assert np.allclose(alpha_with_max_out, alpha_with_max_compare)


def calculate_g_test():
    a1 = np.array([0.2173,1])
    a2 = np.array([0.224,0])
    a3 = np.array([0.2824,0])
    a4 = np.array([0.2763,0])
    tau1 = np.array([1000000,9.15])
    tau2 = np.array([394.4,1])
    tau3 = np.array([36.54,1])
    tau4 = np.array([4.304,1])
    a = np.array([a1, a2, a3, a4])
    tau = np.array([tau1, tau2, tau3, tau4])

    g0_out, g1_out = unifiedtools.calculate_g(a,tau)

    g0_compare = np.array([0.0101837, 0.36780734])
    g1_compare = np.array([11.4137078, 9.1480428])
    assert np.allclose(g0_out, g0_compare)
    assert np.allclose(g1_out, g1_compare)

def step_concentration_test():

    emissions = np.array([10.2123, 361])
    a1 = np.array([0.2173,1])
    a2 = np.array([0.224,0])
    a3 = np.array([0.2824,0])
    a4 = np.array([0.2763,0])
    tau1 = np.array([1000000,9.15])
    tau2 = np.array([394.4,1])
    tau3 = np.array([36.54,1])
    tau4 = np.array([4.304,1])
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
    emis2conc = np.array([0.468952344,0.351714258])


    C_out, R_out, G_A_out = unifiedtools.step_concentration(emissions,a,dt,alpha,tau,R_old, G_A_old, PI_conc,emis2conc)

    C_compare = np.array([383.68002648, 1891.49581172])
    R_compare = np.array([  [85.29456948, 3557.63389782],\
                            [63.31706123, 0.],\
                            [34.02781502, 0.],\
                            [8.767445190, 0.]
                         ]
                         )
    G_A_compare = np.array([191.40689091, 3557.63389782])

    assert np.allclose(C_out, C_compare)
    assert np.allclose(R_out, R_compare)
    assert np.allclose(G_A_out, G_A_compare)


def step_forcing_test():
    C = np.array([399.6,1811])
    PI_conc = np.array([278, 720])
    f1 = np.array([5.754389, 0.061736])
    f2 = np.array([0.001215, -0.000049])
    f3 = np.array([-0.069598, 0.038416])

    RF_out = unifiedtools.step_forcing(C,PI_conc,f1,f2,f3)
    RF_compare = np.array([2.0048501 , 0.60750117])

    assert np.allclose(RF_out, RF_compare)

def step_temperature_test():
    d = np.array([  283, 9.88, 0.85])
    q = np.array([  0.311333, 0.165417, 0.242])
    F = np.array([2.870])
    S_old = np.array([0.173196459, 1.06749276, 2])
    dt = 10

    S_out, T_out = unifiedtools.step_temperature(S_old,F,q,d,dt)

    S_compare = np.array([0.19820533, 0.69017337, 0.69455015])
    T_compare = np.array([2.41180904])

    assert np.allclose(S_out, S_compare)
    assert np.allclose(T_out, T_compare)

def unstep_concentration_test():
    a1 = np.array([0.2173,1])
    a2 = np.array([0.224,0])
    a3 = np.array([0.2824,0])
    a4 = np.array([0.2763,0])
    tau1 = np.array([1000000,9.15])
    tau2 = np.array([394.4,1])
    tau3 = np.array([36.54,1])
    tau4 = np.array([4.304,1])
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

    emissions_out, R_out = unifiedtools.unstep_concentration(a,dt,alpha,tau,R_old, G_A)

    emissions_compare = np.array([19.15739779, 304.0803832])
    R_compare = np.array([  [104.73213702, 3104.0],\
                            [83.01823459, 0.],\
                            [55.18263368, 0.],\
                            [16.36699471, 0.]
                         ]
                         )

    assert np.allclose(emissions_out, emissions_compare)
    assert np.allclose(R_out, R_compare)

def convert_df_to_numpy_test():
    test_df = pd.DataFrame(data = {"a":[7,4,10,1],"c":[9,6,12,3],"b":[8,5,11,2]}, index = [3,2,4,1])
    res_compare = np.array([[ 1,  4,  7, 10], [ 2,  5,  8, 11], [ 3,  6,  9, 12]])
    res = unifiedtools.convert_df_to_numpy_test(test_df)
    assert np.allclose(res_compare, res)

def convert_numpy_output_to_df_test():
    res = np.array([[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]])
    column_labels = np.array(['a','b','c'])
    column_name = 'Gas'
    index_labels = np.array([2019,2020,2022,2023,2050])
    index_name = 'Year'
    df = unifiedtools.convert_numpy_output_to_df(res, column_labels, column_name, index_labels, index_name)
    df_compare = pd.DataFrame(data = {"a":[1,2,3,4,5],"b":[6,7,8,9,10],"c":[11,12,13,14,15]}, index = [2019,2020,2022,2023,2050])
    
    assert df.equals(df_compare)
    assert not np.sum(df.index.values != np.array([2019, 2020, 2022, 2023, 2050]))
    assert df.index.name == 'Year'
    assert not np.sum(df.columns.values != np.array(['a','b','c']))
    assert df.columns.name == 'Gas'