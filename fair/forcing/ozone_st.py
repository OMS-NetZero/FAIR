from __future__ import division

import numpy as np
from ..constants import cl_atoms, br_atoms, fracrel

def magicc(C_ODS,
           C0, 
           eta1=-1.46030698e-5,
           eta2=2.05401270e-3,
           eta3=1.03143308):

    Cl = np.array(cl_atoms.aslist)
    Br = np.array(br_atoms.aslist)
    FC = np.array(fracrel.aslist)

    EESC = (np.sum(Cl * 1000.*(C_ODS-C0) * FC/FC[0]) +
             45*np.sum(Br * 1000.*(C_ODS-C0) * FC/FC[0])) * FC[0]

    # EESC takes ODS concentrations in ppb, we provide ppt.
    EESC = np.max((EESC,0))

    F = eta1 * (eta2 * EESC) ** eta3
    return F

