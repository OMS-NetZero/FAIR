from __future__ import division

import pickle
import os
import sys
import numpy as np
from scipy.interpolate import Rbf
from ..constants import molwt


def Stevens(emissions, E_SOx_nat=60, alpha=0.001875, beta=0.634, ref_isSO2=True):
    """Calculates aerosol forcing based on Stevens (2015) that relates sulphate
    aerosol forcing to SOx emissions in a logarithmic fashion.

    Input:
        emissions:   anthropogenic emissions database
    Keywords:
        E_SOx_nat:   natural emissions of SOx in Mt/yr
        alpha:       scaling parameter for ERFari
        beta:        scaling parameter for ERFaci
        ref_isSO2:   True if E_SOx_nat is in units of SO2 rather than S.
    Output:
        F:           aerosol effective radiative forcing
    """

    factor = 1
    if ref_isSO2:
        factor = molwt.SO2/molwt.S
    em_SOx = emissions[:,5] * factor

    ERFari = -alpha * em_SOx
    ERFaci = -beta * np.log(em_SOx/E_SOx_nat + 1)
    F = ERFari + ERFaci
    return F


def aerocom_direct(emissions, beta = np.array(
    [-35.29e-4, 0.0, -5.034e-4, -5.763e-4, 453e-4, -37.83e-4, -10.35e-4]),
    scale_AR5=False):

    """Calculates direct aerosol forcing based on linear relationships between
    emissions and forcing in Aerocom models.

    Reference: Myhre et al., 2013: https://www.atmos-chem-phys.net/13/1853/2013

    If inputs from an RCPs SCEN file are used, the units will be correct.

    Inputs: 
        emissions: (nt x 40) emissions array
    Keywords:
        beta: 7-element array of forcing efficiencies in W m-2 (Mt yr-1)-1 for
            SOx, CO, NMVOC, NOx, BC, OC, NH3 (in that order)
        scale_AR5:       If True, scale the forcing output so that the best
                         estimate forcing in 2011 is -0.45 W/m2 based on 2011
                         emissions from the RCPs. The discrepancy between the
                         Aerocom and RCP values comes from two sources; (1)
                         there is no dust forcing considered in Aerocom/FAIR,
                         and (2) the aerosol semi-direct is accounted for by
                         applying this scaling factor.
    Outputs:
        Forcing time series
    """

    em_SOx, em_CO, em_NMVOC, em_NOx, em_BC, em_OC, em_NH3 = \
        emissions[:,[5, 6, 7, 8, 9, 10, 11]].T

    F_SOx    = beta[0] * em_SOx    * molwt.SO2 / molwt.S
    F_CO     = beta[1] * em_CO
    F_NMVOC  = beta[2] * em_NMVOC
    F_NOx    = beta[3] * em_NOx    * molwt.NO / molwt.S
    F_BC     = beta[4] * em_BC
    F_OC     = beta[5] * em_OC
    F_NH3    = beta[6] * em_NH3

    F = F_SOx+F_CO+F_NMVOC+F_NOx+F_BC+F_OC+F_NH3

    return 1.500 * scale_AR5 * F


def ghan_indirect_emulator(emissions, fix_pre1850_RCP=True,
    scale_AR5=False):
    """Estimates the aerosol indirect effect based on the simple model in
    Ghan et al., (2013), doi:10.1002/jgrd.50567.

    This function is just an emulator - a full implementation in Python of the
    Ghan routine (originally coded in Fortran) exists, but will require
    optimisation before it can be used in FAIR. I hope to make the full version
    available in a future version.

    A 100-member Latin Hypercube sample of emissions of SOx, NMVOC, BC and OC
    was prepared offline and run through the Ghan simple model. A radial basis
    function interpolator then estimates the radiative forcing based on the 
    sampled points.

    A word of caution: this function may be unreliable outside of the points
    on which it has been trained - generally anything greater than about 5x
    present day emissions (of course there is no guarantee that the underlying
    model will perform well in such high emission scenarios, either) - but may
    also struggle for other parameter combinations as it has not been
    extensively tested. Use common sense at all times!

    Inputs:
        emissions: (nt x 40) numpy emissions array
    Keywords:
        fix_pre1850_RCP: Use different relationship for 1750/65 to 1850 based
                         on anthropogenic emissions from Skeie et al (2011)
                         for 1750 (atmos-chem-phys.net/11/11827/2011)
        scale_AR5:       If True, scale the forcing output so that the best
                         estimate forcing in 2011 is -0.45 W/m2 based on 2011
                         emissions from the RCPs. The Ghan emulator is built on
                         results from the CAM5 GCM. As reported in AR5 WG1 Ch7,
                         GCMs tend to overestimate forcing from aerosol-cloud
                         interactions.
    Outputs:
        Forcing timeseries
    """

    year, em_SOx, em_NMVOC, em_BC, em_OC = emissions[:,[0, 5, 7, 9, 10]].T

    # load up prepared radial basis function. Thanks to
    # http://atucla.blogspot.co.uk/2016/01/save-and-load-rbf-object-fromto-file.html
    RBFfile = open(os.path.join(os.path.dirname(__file__),
        'ghan_emulator.pickle'),'rb')
    if sys.version_info > (3,0):
        RBFunpickler = pickle.Unpickler(RBFfile, encoding='latin1')
    else:
        RBFunpickler = pickle.Unpickler(RBFfile)
    RBFdict = RBFunpickler.load()
    RBFfile.close()

    # This is a dummy but creates an Rbf with 4 predictors and a response
    ghan_emulator = Rbf(np.array([1,2,3]), np.array([10,20,30]), 
        np.array([100,200,300]), np.array([1000,2000,3000]),
        function = RBFdict['function'])

    for key, value in RBFdict.iteritems():
        ghan_emulator.__setattr__(key,value) 

    # PI forcing was not zero as there were some emissions. Use estimates
    # from Skeie et al, 2011 for 1750 forcing. NMVOC is estimated as half of
    # 1850 - different sources use different measurement units. As for 
    # tropospheric ozone a fix can be applied
    F_1750 = ghan_emulator(1, 5, 1.2, 10)
    if isinstance(year, np.ndarray):
        F_pdtotal = np.zeros_like(year)
        for i in range(len(year)):
            if year[i]>=1850 or fix_pre1850_RCP==False:
                F_pdtotal[i] = ghan_emulator(em_SOx[i], 
                                             em_NMVOC[i], 
                                             em_BC[i], 
                                             em_OC[i])
            else:
                F_1850 = ghan_emulator(2.34586, 68.6829, 3.09885, 22.0414)
                F_pdtotal[i] = ghan_emulator(1+em_SOx[i]/2.34586*(2.34586-1),
                                             5+em_NMVOC[i]/68.6829*(68.6829-5),
                                             1.2+em_BC[i]/3.09855*(3.09855-1.2),
                                             10+em_OC[i]/22.0414*(22.0414-10))
    else:
        if year>=1850 or fix_pre1850_RCP==False:
            F_pdtotal = ghan_emulator(em_SOx, em_NMVOC, em_BC, em_OC)
        else:
            F_pdtotal = ghan_emulator(1+em_SOx/2.34586*(2.34586-1),
                                      5+em_NMVOC/68.6829*(68.6829-5),
                                      1.2+em_BC/3.09855*(3.09855-1.2),
                                      10+em_OC/22.0414*(22.0414-10))

    return 0.392 * scale_AR5 * (F_pdtotal-F_1750)
