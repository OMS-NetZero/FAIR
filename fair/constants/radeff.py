# all efficiencies in W/m2/ppb from AR5 WG1 Ch8 SM
# except those marked * from Hodnebrog et al., 2013
CO2         = 1.37E-5     # Hodnebrog paper: 1.3665E-5
CH4         = 3.63E-4
N2O         = 3.00E-3
CF4         = 0.09
C2F6        = 0.25
C3F8        = 0.27730479  # *
C4F10       = 0.36298028  # *
C5F12       = 0.40512062  # *
C6F14       = 0.44
C7F16       = 0.50171525  # *
C8F18       = 0.55174034  # *
C_C4F8      = 0.31526465  # *
HFC23       = 0.18
HFC32       = 0.11
HFC43_10    = 0.42
HFC43_10MEE = HFC43_10
HFC125      = 0.23
HFC134A     = 0.16
HFC143A     = 0.16
HFC152A     = 0.10
HFC227EA    = 0.26
HFC236FA    = 0.24309707   # *
HFC245FA    = 0.24
HFC365MFC   = 0.22261699   # *
SF6         = 0.57
NF3         = 0.20467908   # *
SO2F2       = 0.20119879   # *
CFC11       = 0.26
CFC12       = 0.32
CFC113      = 0.30
CFC114      = 0.31
CFC115      = 0.20
CARB_TET    = 0.17
CCL4        = CARB_TET
MCF         = 0.07
CH3CCL3     = MCF
HCFC22      = 0.21
HCFC141B    = 0.16
HCFC142B    = 0.19
HALON1211   = 0.29
HALON1202   = 0.27
HALON1301   = 0.30
HALON2402   = 0.31
CH3BR       = 0.004
CH3CL       = 0.01
CH2CL2      = 0.0312       # *
CHCL3       = 0.07796969   # *

# This is the list of gases included in the RCPs/AR5/CMIP5.
aslist   = [CO2, CH4, N2O, CF4, C2F6, C6F14, HFC23, HFC32, HFC43_10, HFC125,
            HFC134A, HFC143A, HFC227EA, HFC245FA, SF6, CFC11, CFC12, CFC113,
            CFC114, CFC115, CARB_TET, MCF, HCFC22, HCFC141B, HCFC142B,
            HALON1211, HALON1202, HALON1301, HALON2402, CH3BR, CH3CL]

