# Fractional stratospheric release for ODSs.
# References:
# Daniel, J. and Velders, G.: A focus on information and options for 
#   policymakers, in: Scientific Assessment of Ozone Depletion, WMO, 2011
# Newman et al., 2007: A new formulation of equivalent effective stratospheric
#   chlorine (EESC)

import numpy as np

CFC11     = 0.47
CFC12     = 0.23
CFC113    = 0.29
CFC114    = 0.12
CFC115    = 0.04
CARB_TET  = 0.56
MCF       = 0.67
HCFC22    = 0.13
HCFC141B  = 0.34
HCFC142B  = 0.17
HALON1211 = 0.62
HALON1202 = 0.62
HALON1301 = 0.28
HALON2402 = 0.65
CH3BR     = 0.60
CH3CL     = 0.44
CH2CL2    = np.nan # no literature value available
CHCL3     = np.nan # no literature value available

# This is the list of gases included in the RCPs/AR5/CMIP5.
aslist    = [CFC11, CFC12, CFC113, CFC114, CFC115, CARB_TET, MCF, HCFC22,
             HCFC141B, HCFC142B, HALON1211, HALON1202, HALON1301, HALON2402,
             CH3BR, CH3CL]
