from __future__ import division

def linear(emissions, E_ref=8.09, F_ref=0.04):
    E_BC = emissions[:,9]
    return E_BC * F_ref/E_ref
