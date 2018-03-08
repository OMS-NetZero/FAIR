from __future__ import division

import numpy as np

def etminan(C, Cpi, F2x=3.71):
  """Calculate the radiative forcing from CO2, CH4 and N2O.

  This function uses the updated formulas of Etminan et al. (2016),
  including the overlaps between CO2, methane and nitrous oxide.

  Reference: Etminan et al, 2016, JGR, doi: 10.1002/2016GL071930

  Inputs:
    C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]
    Cpi: pre-industrial [CO2, CH4, N2O] concentrations

  Keywords:
    F2x: radiative forcing from a doubling of CO2.

  Returns:
    3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]
  """

  Cbar = 0.5 * (C[0] + Cpi[0])
  Mbar = 0.5 * (C[1] + Cpi[1])
  Nbar = 0.5 * (C[2] + Cpi[2])

  # Tune the coefficient of CO2 forcing to acheive desired F2x, using 
  # pre-industrial CO2 and N2O. F2x_etminan ~= 3.801.
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


def myhre(C, Cpi, F2x=3.71):
  """Calculate the radiative forcing from CO2, CH4 and N2O.

  This uses the Myhre et al. (1998) relationships including the band
  overlaps between CH4 and N2O. It is also used in AR5.

  Reference: Myhre et al, 1998, JGR, doi: 10.1029/98GL01908

  Inputs:
    C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]
    Cpi: pre-industrial [CO2, CH4, N2O] concentrations

  Keywords:
    F2x: radiative forcing from a doubling of CO2.

  Returns:
    3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]
  """

  F = np.zeros(3)

  F[0] = F2x/np.log(2) * np.log(C[0]/Cpi[0])
  F[1] = 0.036 * (np.sqrt(C[1]) - np.sqrt(Cpi[1])) - (
    MN(C[1],Cpi[2]) - MN(Cpi[1],Cpi[2]))
  F[2] = 0.12 * (np.sqrt(C[2]) - np.sqrt(Cpi[2])) - (
    MN(Cpi[1],C[2]) - MN(Cpi[1],Cpi[2])) 

  return F
