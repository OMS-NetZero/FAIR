# all lifetimes in years

import numpy as np

CO2         = np.nan # dummy - to preserve indexing order
CH4         = 9.3   # for constant lifetime runs
N2O         = 121.
CF4         = 50000.
C2F6        = 10000.
C3F8        = 2600. 
C4F10       = 2600.
C5F12       = 4100.
C6F14       = 3100.
C7F16       = 3000.
C8F18       = 3000.
C_C4F8      = 3200.
HFC23       = 222.
HFC32       = 5.2
HFC43_10    = 16.1
HFC43_10MEE = HFC43_10
HFC125      = 28.2
HFC134A     = 13.4
HFC143A     = 47.1
HFC152A     = 1.5
HFC227EA    = 38.9
HFC236FA    = 242.
HFC245FA    = 7.7
HFC365MFC   = 8.7
SF6         = 3200.
NF3         = 500.
SO2F2       = 36.
CFC11       = 45.
CFC12       = 100.
CFC113      = 85.
CFC114      = 190.
CFC115      = 1020.
CARB_TET    = 26.
CCL4        = CARB_TET
MCF         = 5.
CH3CCL3     = MCF
HCFC22      = 11.9
HCFC141B    = 9.2
HCFC142B    = 17.2
HALON1211   = 16.
HALON1202   = 2.9
HALON1301   = 65.
HALON2402   = 20.
CH3BR       = 0.8
CH3CL       = 1.
CH2CL2      = 0.3945  # from Hodnebrog et al., 2013
CHCL3       = 0.4082  # from Hodnebrog et al., 2013

# This is the list of gases included in the RCPs/AR5/CMIP5.
aslist   = [CO2, CH4, N2O, CF4, C2F6, C6F14, HFC23, HFC32, HFC43_10, HFC125,
            HFC134A, HFC143A, HFC227EA, HFC245FA, SF6, CFC11, CFC12, CFC113,
            CFC114, CFC115, CARB_TET, MCF, HCFC22, HCFC141B, HCFC142B,
            HALON1211, HALON1202, HALON1301, HALON2402, CH3BR, CH3CL]
