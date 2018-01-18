from __future__ import division

import numpy as np
from ..constants import molwt

def from_aviNOx(emissions, frac, E_ref=2.946, F_ref=0.0448, ref_isNO2=True):
    """Calculates contrail radiative forcing from emissions of aviation NOx.

    Inputs:
        emissions:  Raw emissions data.
        frac:       fraction of total NOx emissions due to aviation.
    Keywords:
        E_ref:      reference-year emissions of aviation NOx, Mt/yr. 2.946 
                    is 2005 emissions of aviation NOx in RCP4.5 measured in
                    Mt NO2/yr.
        F_ref:      Forcing from linear persistent contrails + contrail 
                    induced cirrus. The default of 0.0448 W/m2 for 2005 is
                    from Lee et al, 2009 (Atmos. Environ.,
                    doi:10.1016/j.atmosenv.2009.04.024).
        ref_isNO2:  True if E_ref is in units of NO2 rather than N.
    """

    em_NOx   = emissions[:,8]
    factor = 1

    if ref_isNO2:
        factor = molwt.NO2/molwt.N
    return em_NOx*frac * F_ref/E_ref * factor
