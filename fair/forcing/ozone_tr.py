import numpy as np

def regress(emissions,
            aCH4   =  2.8249e-4,
            aCO    =  1.0695e-4,
            aNMVOC = -9.3604e-4,
            aNOx   = 99.7831e-4):

    em_CH4   = emissions[:,3]
    em_CO    = emissions[:,6]
    em_NMVOC = emissions[:,7]
    em_NOx   = emissions[:,8]

    em       = np.array([em_CH4, em_CO, em_NMVOC, em_NOx]).T
    beta     = np.array([aCH4, aCO, aNMVOC, aNOx])

    F = np.sum(beta * em, axis=-1)
    return F

