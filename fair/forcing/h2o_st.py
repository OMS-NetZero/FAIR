from __future__ import division

import numpy as np

def linear(F_CH4, ratio=0.15):
    """Calculates radiative forcing from oxidation of methane to H2O.

    Stratospheric water vapour forcing follows a practically linear
    relationship with the CH4 radiative forcing in MAGICC and AR5.
    """

    F_H2O = ratio * F_CH4
    return F_H2O

