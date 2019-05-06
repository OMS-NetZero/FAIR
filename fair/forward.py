from __future__ import division

import inspect
import numpy as np
import warnings
from scipy.optimize import root
from .ancil import natural, cmip6_volcanic, cmip6_solar, historical_scaling
from .constants import molwt, lifetime, radeff
from .constants.general import M_ATMOS, ppm_gtc
from .defaults import carbon, thermal
from .forcing import ozone_tr, ozone_st, h2o_st, contrails, aerosols, bc_snow,\
                                         landuse
from .forcing.ghg import co2_log


def iirf_interp(alp_b,a,tau,iirf_h,targ_iirf):
    """Interpolation function for finding alpha, the CO2 decay time constant
    scaling factor, in iirf_h equation. See Eq. (7) of Millar et al ACP (2017).

    Inputs:
        alp_b    : Guess for alpha, the scale factor, for tau 
        a        : partition fractions for CO2 boxes
        tau      : time constants for CO2 boxes
        iirf_h   : time horizon for time-integrated airborne fraction
        targ_iirf: iirf_h calculated using simple parameterisation (Eq. (8),
                   Millar et al (2017)).
    """

    iirf_arr = alp_b*(np.sum(a*tau*(1.0 - np.exp(-iirf_h/(tau*alp_b)))))
    return iirf_arr - targ_iirf

    
def iirf_simple(c_acc, temp, r0, rc, rt, iirf_max):
    """Simple linear iIRF relationship. Eq. (8) of Millar et al ACP (2017).
    
    Inputs:
        c_acc    : cumulative airborne carbon anomaly (GtC) since
                   pre-industrial
        temp     : temperature anomaly since pre-industrial
        r0       : pre-industrial time-integrated airborne fraction (yr)
        rc       : sensitivity of time-integrated airborne fraction to airborne
                   carbon (yr/GtC)
        rt       : sensitivity of time-integrated airborne fraction to
                   temperature (yr/K)
        iirf_max : maximum value of time-integrated airborne fraction (yr)
    
    Outputs:
        iirf     : time-integrated airborne fraction of carbon (yr)
    """
    
    return np.min([r0 + rc * c_acc + rt * temp, iirf_max])
    
    
def calculate_q(tcrecs, d, f2x, tcr_dbl, nt):
    """If TCR and ECS are supplied, calculate the q model coefficients.
    See Eqs. (4) and (5) of Millar et al ACP (2017).
    
    Inputs:
        tcrecs  : 2-element array of transient climate response (TCR) and
                  equilibrium climate sensitivity (ECS).
        d       : The slow and fast thermal response time constants
        f2x     : Effective radiative forcing from a doubling of CO2
        tcr_dbl : time to a doubling of CO2 under 1% per year CO2 increase, yr
        nt      : number of timesteps
        
    Outputs:
        q       : coefficients of slow and fast temperature change in each
                  timestep ((nt, 2) array).
    """
    
    # TODO:
    # error checking before call
    # benchmark one call per timestep and if not slower do not convert to 2D
    #  - will make code cleaner
    
    k = 1.0 - (d/tcr_dbl)*(1.0 - np.exp(-tcr_dbl/d))
    # if ECS and TCR are not time-varying, expand them to 2D array anyway
    if tcrecs.ndim==1:
        if len(tcrecs)!=2:
            raise ValueError(
              "Constant TCR and ECS should be a 2-element array")
        tcrecs = np.ones((nt, 2)) * tcrecs
    elif tcrecs.ndim==2:
        if tcrecs.shape!=(nt, 2):
            raise ValueError(
              "Transient TCR and ECS should be a nt x 2 array")
    q  = (1.0 / f2x) * (1.0/(k[0]-k[1])) * np.array([
        tcrecs[:,0]-tcrecs[:,1]*k[1],tcrecs[:,1]*k[0]-tcrecs[:,0]]).T
    return q
    

def carbon_cycle(e0, c_acc0, temp, r0, rc, rt, iirf_max, time_scale_sf0, a, tau,
    iirf_h, carbon_boxes0, c_pi, c0, e1):
    """Calculates CO2 concentrations from emissions.
    
    Inputs:
        e0            : emissions of CO2 (GtC) in timestep t-1
        c_acc0        : cumulative airborne carbon anomaly (GtC) since
                        pre-industrial, timestep t-1
        temp          : temperature anomaly above pre-industrial (K)
        r0            : pre-industrial time-integrated airborne fraction (yr)
        rc            : sensitivity of time-integrated airborne fraction to 
                        airborne carbon (yr/GtC)
        rt            : sensitivity of time-integrated airborne fraction to
                        temperature (yr/K)
        iirf_max      : maximum value of time-integrated airborne fraction (yr)
        time_scale_sf0: initial guess of alpha scaling factor
        a             : partition coefficient of carbon boxes
        tau           : present-day decay time constants of CO2 (yr)
        iirf_h        : time horizon for time-integrated airborne fraction (yr)
        carbon_boxes0 : carbon stored in each atmospheric reservoir at timestep
                        t-1 (GtC)
        c_pi          : pre-industrial concentration of CO2, ppmv
        c0            : concentration of CO2 in timestep t-1, ppmv
        e1            : emissions of CO2 in timestep t, GtC

    Outputs:
        c1            : concentrations of CO2 in timestep t, ppmv
        c_acc1        : cumulative airborne carbon anomaly (GtC) since
                        pre-industrial, timestep t
        carbon_boxes1 : carbon stored in each atmospheric reservoir at timestep
                        t (GtC)
        time_scale_sf : scale factor for CO2 decay constants
    """
    iirf = iirf_simple(c_acc0, temp, r0, rc, rt, iirf_max)
    time_scale_sf = root(iirf_interp, time_scale_sf0,
      args=(a, tau, iirf_h, iirf))['x']
    tau_new = tau * time_scale_sf
    carbon_boxes1 = carbon_boxes0*np.exp(-1.0/tau_new) + a*e1 / ppm_gtc
    c1 = np.sum(carbon_boxes1) + c_pi
    c_acc1 = c_acc0 + 0.5*(e1 + e0) - (c1 - c0)*ppm_gtc
    return c1, c_acc1, carbon_boxes1, time_scale_sf
    
    
def emis_to_conc(c0, e0, e1, ts, lt, vm):
    """Calculate concentrations of well mixed GHGs from emissions for simple
    one-box model.
    
    Inputs (all can be scalar or 1D arrays for multiple species):
        c0: concentrations in timestep t-1
        e0: emissions in timestep t-1
        e1: emissions in timestep t
        ts: length of timestep. Use 1 for sensible results in FaIR 1.3.
        lt: atmospheric (e-folding) lifetime of GHG
        vm: conversion from emissions units (e.g. Mt) to concentrations units
            (e.g. ppb)
            
    Outputs:
        c1: concentrations in timestep t
    """
    c1 = c0 - c0 * (1.0 - np.exp(-ts/lt)) + 0.5 * ts * (e1 + e0) * vm
    return c1


def forc_to_temp(t0, q, d, f, e=1.0):
    """Calculate temperature from a given radiative forcing.

    Inputs:
        t0: Temperature in timestep t-1
        q: The matrix contributions to slow and fast temperature change
           calculated from ECS and TCR (2 element array)
        d: The slow and fast thermal response time constants (2 element array)
        f: radiative forcing (can be scalar or 1D array representing multiple
           species)

    Keywords:
        e: efficacy factor (default 1); if f is an array, e should be an array
           of the same length.

    Outputs:
        t1: slow and fast contributions to total temperature (2 element array)
        in timestep t
    """
    t1 = t0*np.exp(-1.0/d) + q*(1.0-np.exp((-1.0)/d))*np.sum(f*e)
    return t1


def fair_scm(
    emissions=False,
    emissions_driven=True,
    C=None,
    other_rf=0.0,
    q        = thermal.q,
    tcrecs   = thermal.tcrecs,
    d        = thermal.d,
    F2x      = thermal.f2x,
    tcr_dbl  = thermal.tcr_dbl,
    a        = carbon.a,
    tau      = carbon.tau,
    r0       = carbon.r0,
    rc       = carbon.rc,
    rt       = carbon.rt,
    iirf_max = carbon.iirf_max,
    iirf_h   = carbon.iirf_h,
    C_pi=np.array([278., 722., 273., 34.497] + [0.]*25 + [13.0975, 547.996]),
    restart_in=False,
    restart_out=False,
    F_tropO3 = 0.,
    F_aerosol = 0.,
    F_volcanic=cmip6_volcanic.Forcing.volcanic,
    F_solar=cmip6_solar.Forcing.solar,
    F_contrails=0.,
    F_bcsnow=0.,
    F_landuse=0.,
    aviNOx_frac=0.,
    fossilCH4_frac=0.,
    natural=natural.Emissions.emissions,
    efficacy=np.array([1.]*9 + [3.] + [1.]*3),
    scale=None,
    oxCH4_frac=0.61,
    ghg_forcing="Etminan",
    stwv_from_ch4=None,
    b_aero = np.array([-6.2227e-3, 0.0, -3.8392e-4, -1.16551e-3, 1.601537e-2,
      -1.45339e-3, -1.55605e-3]),
    b_tro3 = np.array([2.8249e-4, 1.0695e-4, -9.3604e-4, 99.7831e-4]),
    ghan_params = np.array([-1.95011431, 0.01107147, 0.01387492]),
    stevens_params = np.array([0.001875, 0.634, 60.]),
    useMultigas=True,
    useStevenson=True,
    lifetimes=False,
    aerosol_forcing="aerocom+ghan",
    scaleAerosolAR5=True,
    fixPre1850RCP=True,
    useTropO3TFeedback=True,
    scaleHistoricalAR5=False,
    contrail_forcing='NOx',
    kerosene_supply=0.,
    landuse_forcing='co2',
    ):

    # is iirf_h < iirf_max? Don't stop the code, but warn user
    if iirf_h < iirf_max:
        warnings.warn('iirf_h=%f, which is less than iirf_max (%f)'
          % (iirf_h, iirf_max), RuntimeWarning)

    # Conversion between ppb/ppt concentrations and Mt/kt emissions
    # in the RCP databases ppb = Mt and ppt = kt so factor always 1e18
    emis2conc = M_ATMOS/1e18*np.asarray(molwt.aslist)/molwt.AIR

    # Funny units for nitrogen emissions - N2O is expressed in N2 equivalent
    n2o_sf = molwt.N2O/molwt.N2
    emis2conc[2] = emis2conc[2] / n2o_sf

    # Convert any list to a numpy array for (a) speed and (b) consistency.
    # Goes through all variables in scope and converts them.
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    for arg_to_check in args:
        if type(values[arg_to_check]) is list:
            exec(arg_to_check + '= np.array(' + arg_to_check + ')')

    # Set up the output timeseries variables depending on options and perform
    # basic sense checks
    if useMultigas:
        ngas = 31
        nF   = 13
        if emissions_driven:
            if type(emissions) is not np.ndarray or emissions.shape[1] != 40:
                raise ValueError(
                  "emissions timeseries should be a nt x 40 numpy array")
            carbon_boxes_shape = (emissions.shape[0], a.shape[0])
            thermal_boxes_shape = (emissions.shape[0], d.shape[0])
            nt = emissions.shape[0]
        else:
            if type(C) is not np.ndarray or C.shape[1] != ngas:
                raise ValueError(
                  "C timeseries should be a nt x %d numpy array" % ngas)
            thermal_boxes_shape = (C.shape[0], d.shape[0])
            nt = C.shape[0]
        if np.isscalar(fossilCH4_frac):
            fossilCH4_frac = np.ones(nt) * fossilCH4_frac
        # If custom gas lifetimes are supplied, use them, else import defaults
        if type(lifetimes) is np.ndarray:
            if len(lifetimes)!=ngas:
                raise ValueError(
                  "custom GHG lifetime array must have " + str(ngas) + 
                  " elements")
        else:
            lifetimes = lifetime.aslist
        # Select the desired GHG forcing relationship and populate 
        # stratospheric water vapour from methane scale factor if not specified
        # by user
        if ghg_forcing.lower()=="etminan":
            from .forcing.ghg import etminan as ghg
            if stwv_from_ch4==None: stwv_from_ch4=0.12
        elif ghg_forcing.lower()=="myhre":
            from .forcing.ghg import myhre as ghg
            if stwv_from_ch4==None: stwv_from_ch4=0.15
        else:
            raise ValueError(
              "ghg_forcing should be 'etminan' (default) or 'myhre'")
            
        # Check natural emissions and convert to 2D array if necessary
        if type(natural) in [float,int]:
            natural = natural * np.ones((nt,2))
        elif type(natural) is np.ndarray:
            if natural.ndim==1:
                if natural.shape[0]!=2:
                    raise ValueError(
                      "natural emissions should be a 2-element or nt x 2 " +
                      "array")
                natural = np.tile(natural, nt).reshape((nt,2))
            elif natural.ndim==2:
                if natural.shape[1]!=2 or natural.shape[0]!=nt:
                    raise ValueError(
                      "natural emissions should be a 2-element or nt x 2 " +
                      "array")
        else:
            raise ValueError(
              "natural emissions should be a scalar, 2-element, or nt x 2 " +
              "array")

        # check scale factor is correct shape. If 1D inflate to 2D
        if scale is None:
            scale = np.ones((nt,nF))
        elif scale.shape[-1]==nF:
            if scale.ndim==2 and scale.shape[0]==nt:
                pass
            elif scale.ndim==1:
                scale = np.tile(scale, nt).reshape((nt,nF))
        else:
            raise ValueError("in multi-gas mode, scale should be None, or a "+
              "(13,) or (nt, 13) array")

        # if scaling the historical time series to match AR5, apply these
        # factors to whatever the user specifies
        if scaleHistoricalAR5:
            scale=scale*historical_scaling.all[:nt,:]

    else:
        ngas = 1
        nF   = 1

        if emissions_driven:
            if type(emissions) is np.ndarray:
                if emissions.ndim != 1:
                    raise ValueError(
                      "In CO2-only mode, emissions should be a 1D array")
                nt = emissions.shape[0]
                carbon_boxes_shape = (nt, a.shape[0])
                thermal_boxes_shape = (nt, d.shape[0])
            elif type(other_rf) is np.ndarray:
                if other_rf.ndim != 1:
                    raise ValueError(
                      "In CO2-only mode, other_rf should be a 1D array")
                nt = other_rf.shape[0]
                carbon_boxes_shape = (nt, a.shape[0])
                thermal_boxes_shape = (nt, d.shape[0])
                emissions = np.zeros(nt)
            else:
                raise ValueError(
                  "Neither emissions or other_rf is defined as a timeseries")

        else:
            if type(C) is not np.ndarray or C.ndim != 1:
                raise ValueError(
                  "In CO2-only mode, concentrations should be a 1D array")
            nt = C.shape[0]
            thermal_boxes_shape = (nt, d.shape[0])
            # expand C to 2D array for consistency with other calcs
            C = C.reshape((nt, 1))

        # check scale factor is correct shape - either scalar or 1D
        # needs try/except really
        if scale is None:
            scale = np.ones(nt)
        elif np.isscalar(scale):
            scale = np.ones(nt) * scale
        elif scale.ndim==1 and scale.shape[0]==nt:
            pass
        else:
            raise ValueError("in CO2-only mode, scale should be None, a "+
              "scalar or a (nt,) array")

        # if scaling the historical time series to match AR5, apply these
        # factors to whatever the user specifies
        if scaleHistoricalAR5:
            scale=scale*historical_scaling.co2[:nt]

    # If TCR and ECS are supplied, calculate q coefficients
    if type(tcrecs) is np.ndarray:
        q = calculate_q(tcrecs, d, F2x, tcr_dbl, nt)

    # Check a and tau are same size
    if a.ndim != 1:
        raise ValueError("a should be a 1D array")
    if tau.ndim != 1:
        raise ValueError("tau should be a 1D array")
    if len(a) != len(tau):
        raise ValueError("a and tau should be the same size")
    if not np.isclose(np.sum(a), 1.0):
        raise ValueError("a should sum to one")

    # Allocate intermediate and output arrays
    F = np.zeros((nt, nF))
    C_acc = np.zeros(nt)
    T_j = np.zeros(thermal_boxes_shape)
    T = np.zeros(nt)
    C_0 = np.copy(C_pi)
    if emissions_driven:
        C = np.zeros((nt, ngas))
        R_i = np.zeros(carbon_boxes_shape)

    if restart_in:
        R_minus1 = restart_in[0]
        T_j_minus1 = restart_in[1]
        C_acc_minus1 = restart_in[2]
        E_minus1 = restart_in[3]
        C_minus1 = np.sum(R_minus1,axis=-1) + C_0[0]

        C[0,0], C_acc[0], R_i[0,:], time_scale_sf = carbon_cycle(
          E_minus1,
          C_acc_minus1,
          np.sum(T_j_minus1),
          r0,
          rc,
          rt,
          iirf_max,
          0.16,
          a,
          tau,
          iirf_h,
          R_minus1,
          C_pi[0],
          C_minus1,
          emissions[0]
        )

        if np.isscalar(other_rf):
            F[0,0] = co2_log(C[0,0], C_pi[0], F2x) + other_rf
        else:
            F[0,0] = co2_log(C[0,0], C_pi[0], F2x) + other_rf[0]

        F[0,0] = F[0,0] * scale[0]

        T_j[0,:] = forc_to_temp(T_j_minus1, q[0,:], d, F[0,:])
        T[0]=np.sum(T_j[0,:],axis=-1)

    else:
        # Initialise the carbon pools to be correct for first timestep in
        # numerical method
        if emissions_driven:
            if useMultigas:
                R_i[0,:] = a * (np.sum(emissions[0,1:3])) / ppm_gtc
                C[0,1:] = C_0[1:]
            else:
                R_i[0,:] = a * emissions[0,np.newaxis] / ppm_gtc
            C[0,0] = np.sum(R_i[0,:],axis=-1) + C_0[0]

    if useMultigas:
        # CO2, CH4 and N2O are co-dependent
        F[0,0:3] = ghg(C[0,0:3], C_pi[0:3], F2x=F2x)
        # Minor (F- and H-gases) are linear in concentration
        # the factor of 0.001 here is because radiative efficiencies are given
        # in W/m2/ppb and concentrations of minor gases are in ppt.
        F[0,3] = np.sum((C[0,3:] - C_pi[3:]) * radeff.aslist[3:] * 0.001)

        # Tropospheric ozone:
        if emissions_driven:
            if useStevenson:
                F[0,4] = ozone_tr.stevenson(emissions[0,:], C[0,1],
                  T=np.sum(T_j[0,:]), 
                  feedback=useTropO3TFeedback,
                  fix_pre1850_RCP=fixPre1850RCP)
            else:
                F[0,4] = ozone_tr.regress(emissions[0,:], beta=b_tro3)
        else:
            F[:,4] = F_tropO3

        # Stratospheric ozone depends on concentrations of ODSs (index 15-30)
        F[0,5] = ozone_st.magicc(C[0,15:], C_pi[15:])

        # Stratospheric water vapour is a function of the methane ERF
        F[0,6] = h2o_st.linear(F[0,1], ratio=stwv_from_ch4)

        # Forcing from contrails. No climate feedback so can live outside
        # of forward model in this version
        if emissions_driven:
            if contrail_forcing.lower()[0]=='n':   # from NOx emissions
                F[:,7] = contrails.from_aviNOx(emissions, aviNOx_frac)
            elif contrail_forcing.lower()[0]=='f': # from kerosene production
                F[:,7] = contrails.from_fuel(kerosene_supply)
            elif contrail_forcing.lower()[0]=='e': # external forcing timeseries
                F[:,7] = F_contrails
            else:
                raise ValueError("contrails must be one of 'NOx' (estimated "+
                 "from NOx emissions), 'fuel' (estimated from annual jet fuel "+
                 "supplied) or 'external' (an external forcing time series).")
        else:
            F[:,7] = F_contrails

        # Forcing from aerosols - again no feedback dependence
        if emissions_driven:
            if aerosol_forcing.lower()=='stevens':
                F[:,8] = aerosols.Stevens(emissions, stevens_params=stevens_params)
            elif 'aerocom' in aerosol_forcing.lower():
                F[:,8] = aerosols.aerocom_direct(emissions, beta=b_aero)
                if 'ghan' in aerosol_forcing.lower():
                    F[:,8] = F[:,8] + aerosols.ghan_indirect(emissions,
                      scale_AR5=scaleAerosolAR5,
                      fix_pre1850_RCP=fixPre1850RCP,
                      ghan_params=ghan_params)
            elif aerosol_forcing.lower()[0] == 'e':
                F[:,8] = F_aerosol
            else:
                raise ValueError("aerosol_forcing should be one of 'stevens', " +
                  "aerocom, aerocom+ghan or external")
        else:
            F[:,8] = F_aerosol

        # Black carbon on snow - no feedback dependence
        if emissions_driven:
            F[:,9] = bc_snow.linear(emissions)
        else:
            F[:,9] = F_bcsnow

        # Land use change - either use a scaling with cumulative CO2 emissions
        # or an external time series
        if emissions_driven:
            if landuse_forcing.lower()[0]=='c':
                F[:,10] = landuse.cumulative(emissions)
            elif landuse_forcing.lower()[0]=='e':
                F[:,10] = F_landuse
            else:
                raise ValueError(
                 "landuse_forcing should be one of 'co2' or 'external'")
        else:
            F[:,10] = F_landuse
            
        # Volcanic and solar copied straight to the output arrays
        F[:,11] = F_volcanic
        F[:,12] = F_solar

        # multiply by scale factors
        F[0,:] = F[0,:] * scale[0,:]

    else:
        if np.isscalar(other_rf):
            F[0,0] = co2_log(C[0,0], C_pi[0], F2x) + other_rf
        else:
            F[0,0] = co2_log(C[0,0], C_pi[0], F2x) + other_rf[0]
        F[0,0] = F[0,0] * scale[0]

    if restart_in == False:
        # Update the thermal response boxes
        T_j[0,:] = (q[0,:]/d)*(np.sum(F[0,:]))

    # Sum the thermal response boxes to get the total temperature anomaly
    T[0]=np.sum(T_j[0,:],axis=-1)

    for t in range(1,nt):

        if emissions_driven:
            if useMultigas:
                if t == 1:
                    time_scale_sf = 0.16
                # Calculate concentrations
                # a. CARBON DIOXIDE
                # Firstly add any oxidised methane from last year to the CO2
                # pool
                oxidised_CH4 = ((C[t-1,1]-C_pi[1]) *
                  (1.0 - np.exp(-1.0/lifetimes[1])) * 
                  (molwt.C/molwt.CH4 * 0.001 * oxCH4_frac * fossilCH4_frac[t]))
                oxidised_CH4 = np.max((oxidised_CH4, 0))

                C[t,0], C_acc[t], R_i[t,:], time_scale_sf = carbon_cycle(
                  np.sum(emissions[t-1,1:3]),
                  C_acc[t-1],
                  T[t-1],
                  r0,
                  rc,
                  rt,
                  iirf_max,
                  time_scale_sf,
                  a,
                  tau,
                  iirf_h,
                  R_i[t-1,:] + oxidised_CH4,
                  C_pi[0],
                  C[t-1,0],
                  np.sum(emissions[t,1:3])
                )

                # b. METHANE
                C[t,1] = emis_to_conc(
                    C[t-1,1],
                    emissions[t-1,3]+natural[t,0], 
                    emissions[t,3]+natural[t,0],
                    1.0,
                    lifetimes[1],
                    1.0/emis2conc[1]
                    )

                # c. NITROUS OXIDE
                C[t,2] = emis_to_conc(
                    C[t-1,2],
                    emissions[t-1,4]+natural[t,1], 
                    emissions[t,4]+natural[t,1],
                    1.0,
                    lifetimes[2],
                    1.0/emis2conc[2]
                    )

                # d. OTHER WMGHGs
                C[t,3:] = emis_to_conc(
                    C[t-1,3:],
                    emissions[t-1,12:], 
                    emissions[t,12:],
                    1.0,
                    np.array(lifetimes[3:]),
                    1.0/emis2conc[3:]
                    )

                # 2. Radiative forcing
                F[t,0:3] = ghg(C[t,0:3], C_pi[0:3], F2x=F2x)
                F[t,3] = np.sum((C[t,3:] - C_pi[3:]) * radeff.aslist[3:]
                  * 0.001)
                if useStevenson:
                    F[t,4] = ozone_tr.stevenson(emissions[t,:],
                      C[t,1],
                      T=T[t-1], 
                      feedback=useTropO3TFeedback,
                      fix_pre1850_RCP=fixPre1850RCP)
                else:
                    F[t,4] = ozone_tr.regress(emissions[t,:], beta=b_tro3)
                F[t,5] = ozone_st.magicc(C[t,15:], C_pi[15:])
                F[t,6] = h2o_st.linear(F[t,1], ratio=stwv_from_ch4)

                # multiply by scale factors
                F[t,:] = F[t,:] * scale[t,:]

                # 3. Temperature
                # Update the thermal response boxes
                T_j[t,:] = forc_to_temp(
                  T_j[t-1,:], q[t,:], d, F[t,:], e=efficacy)
                # Sum the thermal response boxes to get the total temperature
                T[t]=np.sum(T_j[t,:],axis=-1)

            else:
                if t == 1:
                    time_scale_sf = 0.16
                C[t,0], C_acc[t], R_i[t,:], time_scale_sf = carbon_cycle(
                  emissions[t-1],
                  C_acc[t-1],
                  T[t-1],
                  r0,
                  rc,
                  rt,
                  iirf_max,
                  time_scale_sf,
                  a,
                  tau,
                  iirf_h,
                  R_i[t-1,:],
                  C_pi[0],
                  C[t-1,0],
                  emissions[t]
                )
                if np.isscalar(other_rf):
                    F[t,0] = co2_log(C[t,0], C_pi[0], F2x) + other_rf
                else:
                    F[t,0] = co2_log(C[t,0], C_pi[0], F2x) + other_rf[t]

                F[t,0] = F[t,0] * scale[t]

                T_j[t,:] = forc_to_temp(T_j[t-1,:], q[t,:], d, F[t,:])
                T[t]=np.sum(T_j[t,:],axis=-1)

        else:

            if useMultigas:
                F[t,0:3] = ghg(C[t,0:3], C_pi[0:3], F2x=F2x)
                F[t,3] = np.sum((C[t,3:] - C_pi[3:]) * radeff.aslist[3:]
                  * 0.001)
                F[t,5] = ozone_st.magicc(C[t,15:], C_pi[15:])
                F[t,6] = h2o_st.linear(F[t,1], ratio=stwv_from_ch4)

                # multiply by scale factors
                F[t,:] = F[t,:] * scale[t,:]

                # 3. Temperature
                # Update the thermal response boxes
                T_j[t,:] = T_j[t,:] = forc_to_temp(
                  T_j[t-1,:], q[t,:], d, F[t,:], e=efficacy)
                # Sum the thermal response boxes to get the total temperature
                T[t]=np.sum(T_j[t,:],axis=-1)

            else:
                if np.isscalar(other_rf):
                    F[t,0] = co2_log(C[t,0], C_pi[0], F2x) + other_rf
                else:
                    F[t,0] = co2_log(C[t,0], C_pi[0], F2x) + other_rf[t]

                F[t,0] = F[t,0] * scale[t]

                T_j[t,:] = forc_to_temp(T_j[t-1,:], q[t,:], d, F[t,:])
                T[t]=np.sum(T_j[t,:],axis=-1)

    if not useMultigas:
        C = np.squeeze(C)
        F = np.squeeze(F)

    if restart_out:
        if useMultigas:
            E_minus1 = np.sum(emissions[-1,1:3])
        else:
            E_minus1 = emissions[-1]
        restart_out_val=(R_i[-1],T_j[-1],C_acc[-1],E_minus1)
        return C, F, T, restart_out_val
    else:
        return C, F, T
