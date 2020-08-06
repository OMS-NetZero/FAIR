from __future__ import division

import numpy as np

from ..constants import radeff


def meinshausen(
    C,
    Cpi=np.array([277.15, 731.41, 273.87]),
    a1=-2.4785e-07, b1=0.00075906, c1=-0.0021492, d1 = 5.2488,
    a2=-0.00034197, b2 = 0.00025455, c2 =-0.00024357, d2 = 0.12173,
    a3 =-8.9603e-05, b3 = -0.00012462, d3 = 0.045194,
    F2x=3.71, scale_F2x=True
    ):
    """Modified Etminan relationship from Meinshausen et al 2019
    https://gmd.copernicus.org/preprints/gmd-2019-222/gmd-2019-222.pdf
    table 3

    Inputs:
        C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]

    Keywords:
        Cpi: pre-industrial [CO2, CH4, N2O] concentrations. Should use defaults
        a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, d3: coefficients
        F2x: radiative forcing from a doubling of CO2.
        scale_F2x: boolean. Scale the calculated value to the specified F2x?

    Returns:
        3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]

    """
    # Tune the coefficient of CO2 forcing to acheive desired F2x, using
    # pre-industrial CO2 and N2O. F2x_etminan ~= 3.801.
    scaleCO2 = 1
    if scale_F2x:
        F2x_etminan = (
          -2.4e-7*Cpi[0]**2 + 7.2e-4*Cpi[0] - 2.1e-4*Cpi[2] + 5.36) * np.log(2)
        scaleCO2 = F2x/F2x_etminan

    F = np.zeros(3)

    # CO2
    Camax = Cpi[0] - b1/(2*a1)
    if Cpi[0] < C[0] <= Camax: # the most likely case
        alphap = d1 + a1*(C[0] - Cpi[0])**2 + b1*(C[0] - Cpi[0])
    elif C[0] <= Cpi[0]:
        alphap = d1
    else:
        alphap = d1 - b1**2/(4*a1)
    alphaN2O = c1*np.sqrt(C[2])
    F[0] = (alphap + alphaN2O) * np.log(C[0]/Cpi[0]) * scaleCO2

    # CH4
    F[1] = (a3*np.sqrt(C[1]) + b3*np.sqrt(C[2]) + d3) * (np.sqrt(C[1]) - np.sqrt(Cpi[1]))

    # N2O
    F[2] = (a2*np.sqrt(C[0]) + b2*np.sqrt(C[2]) + c2*np.sqrt(C[1]) + d2) * (np.sqrt(C[2]) - np.sqrt(Cpi[2]))

    return F


def etminan(C, Cpi, F2x=3.71, scale_F2x=True):
    """Calculate the radiative forcing from CO2, CH4 and N2O.

    This function uses the updated formulas of Etminan et al. (2016),
    including the overlaps between CO2, methane and nitrous oxide.

    Reference: Etminan et al, 2016, JGR, doi: 10.1002/2016GL071930

    Inputs:
        C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]
        Cpi: pre-industrial [CO2, CH4, N2O] concentrations

    Keywords:
        F2x: radiative forcing from a doubling of CO2.
        scale_F2x: boolean. Scale the calculated value to the specified F2x?

    Returns:
        3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]
    """

    Cbar = 0.5 * (C[0] + Cpi[0])
    Mbar = 0.5 * (C[1] + Cpi[1])
    Nbar = 0.5 * (C[2] + Cpi[2])

    # Tune the coefficient of CO2 forcing to acheive desired F2x, using 
    # pre-industrial CO2 and N2O. F2x_etminan ~= 3.801.
    scaleCO2 = 1
    if scale_F2x:
        F2x_etminan = (
          -2.4e-7*Cpi[0]**2 + 7.2e-4*Cpi[0] - 2.1e-4*Cpi[2] + 5.36) * np.log(2)
        scaleCO2 = F2x/F2x_etminan

    F = np.zeros(3)
    F[0] = (-2.4e-7*(C[0] - Cpi[0])**2 + 7.2e-4*np.fabs(C[0]-Cpi[0]) - \
      2.1e-4 * Nbar + 5.36) * np.log(C[0]/Cpi[0]) * scaleCO2
    F[1] = (-1.3e-6*Mbar - 8.2e-6*Nbar + 0.043) * (np.sqrt(C[1]) - \
      np.sqrt(Cpi[1]))
    F[2] = (-8.0e-6*Cbar + 4.2e-6*Nbar - 4.9e-6*Mbar + 0.117) * \
      (np.sqrt(C[2]) - np.sqrt(Cpi[2]))

    return F


def MN(M, N):
    return 0.47 * np.log(1 + 2.01e-5*(M*N)**(0.75) + 5.31e-15*M*(M*N)**(1.52))


def co2_log(C, Cpi, F2x=3.71):
    return F2x/np.log(2) * np.log(C/Cpi)


def myhre(C, Cpi, F2x=3.71, scale_F2x=None):
# TODO: remove scale_F2x in v1.6
    """Calculate the radiative forcing from CO2, CH4 and N2O.

    This uses the Myhre et al. (1998) relationships including the band
    overlaps between CH4 and N2O. It is also used in AR5.

    Reference: Myhre et al, 1998, JGR, doi: 10.1029/98GL01908

    Inputs:
        C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]
        Cpi: pre-industrial [CO2, CH4, N2O] concentrations

    Keywords:
        F2x: radiative forcing from a doubling of CO2.
        scale_F2x: redundant; included for compatibility on import
            in fair_scm

    Returns:
        3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]
    """

    F = np.zeros(3)

    F[0] = co2_log(C[0], Cpi[0], F2x)
    F[1] = 0.036 * (np.sqrt(C[1]) - np.sqrt(Cpi[1])) - (
      MN(C[1],Cpi[2]) - MN(Cpi[1],Cpi[2]))
    F[2] = 0.12 * (np.sqrt(C[2]) - np.sqrt(Cpi[2])) - (
      MN(Cpi[1],C[2]) - MN(Cpi[1],Cpi[2])) 

    return F


def minor_gases(C, Cpi):
    """
    Calculate radiative forcing from minor gas species.

    Inputs:
        C: concentration of minor GHGs (in order of MAGICC RCP concentration
            spreadsheets)
        Cpi: Pre-industrial concentrations of GHGs
    Returns:
        28 element array of minor GHG forcings
    """

    return (C - Cpi) * radeff.aslist[3:] * 0.001
