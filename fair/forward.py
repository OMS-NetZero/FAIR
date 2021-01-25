from __future__ import division

import inspect
import numpy as np
import warnings

from .ancil import natural, cmip6_volcanic, cmip6_solar, historical_scaling
from .constants import molwt, lifetime, radeff
from .constants.general import M_ATMOS, ppm_gtc
from .defaults import carbon, thermal
from .forcing import (ozone, ozone_tr, ozone_st, h2o_st, contrails, aerosols,
    bc_snow, landuse)
from .gas_cycle.gir import calculate_alpha, step_concentration
from .gas_cycle.fair1 import carbon_cycle
from .forcing.ghg import co2_log, minor_gases
from .temperature.millar import calculate_q



# TODO: unified interface to the different carbon cycles


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
    E_pi=np.zeros(40),
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
    F_ref_aviNOx=0.0448,
    E_ref_aviNOx=2.946,
    F_ref_BC=0.04,
    E_ref_BC=8.09,
    fossilCH4_frac=0.,
    natural=natural.Emissions.emissions,
    efficacy=np.array([1.]*9 + [3.] + [1.]*3),
    scale=None,
    oxCH4_frac=0.61,
    ghg_forcing="Etminan",
    scale_F2x=True,
    stwv_from_ch4=None,
    b_aero = np.array([-6.2227e-3, 0.0, -3.8392e-4, -1.16551e-3, 1.601537e-2,
      -1.45339e-3, -1.55605e-3]),
    b_tro3 = np.array([2.8249e-4, 1.0695e-4, -9.3604e-4, 99.7831e-4]),
    pi_tro3 =np.array([722, 170, 10, 4.29]),
    ghan_params = np.array([-1.95011431, 0.01107147, 0.01387492]),
    stevens_params = np.array([0.001875, 0.634, 60.]),
    ref_isSO2=True, # is Stevens SO2 emissions in units SO2 (T) or S (F)
    useMultigas=True,
    tropO3_forcing='stevenson',
    ozone_feedback=-0.037, # W m-2 K-1, only for Thornhill-Skeie
    lifetimes=False,
    aerosol_forcing="aerocom+ghan",
    scaleAerosolAR5=True,
    fixPre1850RCP=True,
    useTropO3TFeedback=True,
    scaleHistoricalAR5=False,
    contrail_forcing='NOx',
    kerosene_supply=0.,
    landuse_forcing='co2',
    aCO2land=-0.00113789,
    ariaci_out=False,
    bcsnow_forcing='emissions',
    diagnostics=None,
    gir_carbon_cycle=False,
    temperature_function='Millar',
    lambda_global=1.18,  # this and the below only used in two-layer model
    ocean_heat_capacity=np.array([8.2, 109.0]),
    ocean_heat_exchange=0.67,
    deep_ocean_efficacy=1.28,
    ):

    # Prevents later errors when SLCFs not specified
    if type(emissions) is bool and not emissions_driven:
        tropO3_forcing='external'

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

    # Initialise simplified carbon cycle parameters
    if gir_carbon_cycle:
        g1 = np.sum(a*tau * (1 - (1 + iirf_h/tau) * np.exp(-iirf_h/tau)))
        g0 = 1/(np.sinh(np.sum(a*tau*(1 - np.exp(-iirf_h/tau)) , axis=-1)/g1))
        if useMultigas:
            cumulative_emissions = np.cumsum(emissions[:,1:3].sum(axis=1))
        else:
            cumulative_emissions = np.cumsum(emissions)
        airborne_emissions = np.zeros_like(cumulative_emissions)

    # import correct conversion
    if temperature_function=='Millar':
        from .temperature.millar import forcing_to_temperature
    elif temperature_function=='Geoffroy':
        from .temperature.geoffroy import forcing_to_temperature
    else:
        raise ValueError('temperature_function must be "Millar" or "Geoffroy"')
        
    # Set up the output timeseries variables depending on options and perform
    # basic sense checks
    if useMultigas:
        ngas = 31
        if diagnostics=='AR6':
            iF_tro3 = 31
            iF_sto3 = 32
            iF_ch4h = 33
            iF_cont = 34
            iF_adso = 35
            iF_advo = 36
            iF_adni = 37
            iF_adbc = 38
            iF_adoc = 39
            iF_aeri = 40
            iF_bcsn = 41
            iF_luch = 42
            iF_volc = 43
            iF_solr = 44
            nF = 45
        else:
            iF_tro3 = 4
            iF_sto3 = 5
            iF_ch4h = 6
            iF_cont = 7
            iF_aer  = 8
            iF_bcsn = 9
            iF_luch = 10
            iF_volc = 11
            iF_solr = 12
            nF = 13
        if emissions_driven:
            if type(emissions) is not np.ndarray or emissions.shape[1] != 40:
                raise ValueError(
                  "emissions timeseries should be a nt x 40 numpy array")
            carbon_boxes_shape = (emissions.shape[0], a.shape[0])
            if temperature_function=='Millar':
                thermal_boxes_shape = (emissions.shape[0], d.shape[0])
            else:
                thermal_boxes_shape = (emissions.shape[0], d.shape[0], 2)
            nt = emissions.shape[0]
        else:
            if type(C) is not np.ndarray or C.shape[1] != ngas:
                raise ValueError(
                  "C timeseries should be a nt x %d numpy array" % ngas)
            if temperature_function=='Millar':
                thermal_boxes_shape = (C.shape[0], d.shape[0])
            else:
                thermal_boxes_shape = (C.shape[0], d.shape[0], 2)
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
        elif ghg_forcing.lower()=="meinshausen":
            from .forcing.ghg import meinshausen as ghg
            if stwv_from_ch4==None: stwv_from_ch4=0.12
        else:
            raise ValueError(
              "ghg_forcing should be 'etminan' (default), 'meinshausen' or 'myhre'")
        # aerosol breakdown
        ariaci = np.zeros((nt,2))
            
        # Check natural emissions and convert to 2D array if necessary
        if emissions_driven: # don't check for conc runs
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
              "(%d,) or (%d, 13) array" % (nF, nF))

        # if scaling the historical time series to match AR5, apply these
        # factors to whatever the user specifies
        if scaleHistoricalAR5:
            scale=scale*historical_scaling.all[:nt,:]

        # if tropospheric ozone is directly specified and scalar, inflate to
        # 1D array. Raise ValueError if wrong shape
        if tropO3_forcing[0].lower()=='e':
            if type(F_tropO3) is np.ndarray:
                if F_tropO3.shape[0]!=nt or F_tropO3.ndim!=1:
                    raise ValueError("F_tropO3 should be a scalar or (nt,) "+
                    "array")
            elif type(F_tropO3) in [float,int]:
                F_tropO3 = F_tropO3 * np.ones(nt)
            else:
                raise ValueError("F_tropO3 should be a scalar or (nt,) array")

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
                if temperature_function=='Millar':
                    thermal_boxes_shape = (nt, d.shape[0])
                else:
                    thermal_boxes_shape = (nt, d.shape[0], 2)
            elif type(other_rf) is np.ndarray:
                if other_rf.ndim != 1:
                    raise ValueError(
                      "In CO2-only mode, other_rf should be a 1D array")
                nt = other_rf.shape[0]
                carbon_boxes_shape = (nt, a.shape[0])
                if temperature_function=='Millar':
                    thermal_boxes_shape = (nt, d.shape[0])
                else:
                    thermal_boxes_shape = (nt, d.shape[0], 2)
                emissions = np.zeros(nt)
            else:
                raise ValueError(
                  "Neither emissions or other_rf is defined as a timeseries")

        else:
            if type(C) is not np.ndarray or C.ndim != 1:
                raise ValueError(
                  "In CO2-only mode, concentrations should be a 1D array")
            nt = C.shape[0]
            if temperature_function=='Millar':
                thermal_boxes_shape = (nt, d.shape[0])
            else:
                thermal_boxes_shape = (nt, d.shape[0], 2)
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
    if type(tcrecs) is np.ndarray and temperature_function=='Millar':
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

    if temperature_function!='Millar':
        heatflux = np.zeros(nt)
        ohc = np.zeros(nt)
        lambda_eff = np.zeros(nt)

    if restart_in:
        R_minus1 = restart_in[0]
        T_j_minus1 = restart_in[1]
        C_acc_minus1 = restart_in[2]
        E_minus1 = restart_in[3]
        C_minus1 = np.sum(R_minus1,axis=-1) + C_0[0]

        if gir_carbon_cycle:
            raise NotImplementedError('GIR carbon cycle not configured to ' +
                'work with restarts')
        else:
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

        if temperature_function=='Millar':
            T_j[0,:] = forcing_to_temperature(T_j_minus1, q[0,:], d, F[0,:])
            T[0]=np.sum(T_j[0,:])
        else:
            # leave unimplemented unless somebody invents a use case
            raise(NotImplementedError('Restarts not implemented with Geoffroy '+
                'temperature function'))

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
        F[0,0:3] = ghg(C[0,0:3], C_pi[0:3], F2x=F2x, scale_F2x=scale_F2x)
        # Minor (F- and H-gases) are linear in concentration
        # the factor of 0.001 here is because radiative efficiencies are given
        # in W/m2/ppb and concentrations of minor gases are in ppt.
        if diagnostics=='AR6':
            F[0,3:31] = minor_gases(C[0,3:], C_pi[3:])
        else:
            F[0,3] = np.sum(minor_gases(C[0,3:], C_pi[3:]))

        # Tropospheric ozone:
        # v1.5 update: don't require emissions driven runs here for GHGs
        # because SLCFs can still be given as emissions with GHGs as
        # concentrations
        if type(emissions) is not bool:
            if tropO3_forcing[0].lower()=='s':  # stevenson
                F[0,iF_tro3] = ozone_tr.stevenson(emissions[0,:], C[0,1],
                  T=np.sum(T_j[0,:]), 
                  feedback=useTropO3TFeedback,
                  fix_pre1850_RCP=fixPre1850RCP,
                  PI=pi_tro3)
            elif tropO3_forcing[0].lower()=='c':  # cmip6 stevenson
                F[0,iF_tro3] = ozone_tr.cmip6_stevenson(emissions[0,:], C[0,1],
                  T=np.sum(T_j[0,:]),
                  feedback=useTropO3TFeedback,
                  PI=np.array([C_pi[1],E_pi[6],E_pi[7],E_pi[8]]),
                  beta=b_tro3)
            elif tropO3_forcing[0].lower()=='r':  # regression
                F[0,iF_tro3] = ozone_tr.regress(emissions[0,:]-E_pi[:], beta=b_tro3)
            elif tropO3_forcing[0].lower()=='t':  # thornhill-skeie
                F[0,iF_tro3] = ozone.thornhill_skeie(
                    emissions=emissions[0,:],
                    concentrations=C[0,:],
                    temperature=T[0],
                    feedback=ozone_feedback,
                    beta=b_tro3,
                    emissions_pi=E_pi,
                    concentrations_pi=C_pi,
                )
            else:
                F[0,iF_tro3] = F_tropO3[0]
        else:
            F[0,iF_tro3] = F_tropO3[0]

        if tropO3_forcing[0].lower()!='t':
            # Stratospheric ozone depends on concentrations of ODSs (index 15-30)
            F[0,iF_sto3] = ozone_st.magicc(C[0,15:], C_pi[15:])

        # Stratospheric water vapour is a function of the methane ERF
        F[0,iF_ch4h] = h2o_st.linear(F[0,1], ratio=stwv_from_ch4)

        # Forcing from contrails. No climate feedback so can live outside
        # of forward model in this version
        # v1.5 update: don't require emissions driven runs here for GHGs
        # because SLCFs can still be given as emissions with GHGs as
        # concentrations
        if type(emissions) is not bool:
            if contrail_forcing.lower()[0]=='n':   # from NOx emissions
                F[:,iF_cont] = contrails.from_aviNOx(emissions, aviNOx_frac,
                  F_ref=F_ref_aviNOx, E_ref=E_ref_aviNOx)
            elif contrail_forcing.lower()[0]=='f': # from kerosene production
                F[:,iF_cont] = contrails.from_fuel(kerosene_supply)
            elif contrail_forcing.lower()[0]=='e': # external forcing timeseries
                F[:,iF_cont] = F_contrails
            else:
                raise ValueError("contrails must be one of 'NOx' (estimated "+
                 "from NOx emissions), 'fuel' (estimated from annual jet fuel "+
                 "supplied) or 'external' (an external forcing time series).")
        else:
            F[:,iF_cont] = F_contrails

        # Forcing from aerosols - again no feedback dependence
        # v1.5 update: don't require emissions driven runs here for GHGs
        # because SLCFs can still be given as emissions with GHGs as
        # concentrations

        # TODO: this whole code block is a mess!
        if type(emissions) is not bool:
            if aerosol_forcing.lower()=='stevens':
                ariaci[:,0], ariaci[:,1] = aerosols.Stevens(
                  emissions, stevens_params=stevens_params, E_pi=E_pi[5], 
                  ref_isSO2=ref_isSO2)
                if diagnostics=='AR6':
                    F[:,iF_adso] = ariaci[:,0]
                    F[:,iF_aeri] = ariaci[:,1]
                else:
                    F[:,iF_aer] = np.sum(ariaci, axis=1)
            elif 'aerocom' in aerosol_forcing.lower():
                aerosol_direct = aerosols.aerocom_direct(emissions, beta=b_aero,
                  E_pi=E_pi, diagnostics=diagnostics)
                if diagnostics=='AR6':
                    F[:,iF_adso] = aerosol_direct[:,0]
                    F[:,iF_advo] = aerosol_direct[:,1] + aerosol_direct[:,2]
                    F[:,iF_adni] = aerosol_direct[:,3] + aerosol_direct[:,6]
                    F[:,iF_adbc] = aerosol_direct[:,4]
                    F[:,iF_adoc] = aerosol_direct[:,5]
                    ariaci[:,0] = np.sum(aerosol_direct, axis=1)
                else:
                    ariaci[:,0] = aerosol_direct
                if 'ghan2' in aerosol_forcing.lower():
                    ariaci[:,1] = aerosols.ghan2(emissions, E_pi, ghan_params)
                elif 'ghan' in aerosol_forcing.lower():
                    ariaci[:,1] = aerosols.ghan_indirect(emissions,
                      scale_AR5=scaleAerosolAR5,
                      fix_pre1850_RCP=fixPre1850RCP,
                      ghan_params=ghan_params,
                      E_pi=E_pi)
                elif 'stevens' in aerosol_forcing.lower():
                    _, ariaci[:,1] = aerosols.Stevens(
                      emissions, stevens_params=stevens_params, E_pi=E_pi[5],
                      ref_isSO2=ref_isSO2)
                if diagnostics=='AR6':
                    F[:,iF_aeri] = ariaci[:,1]
                else:
                    F[:,iF_aer] = np.sum(ariaci, axis=1)
            elif aerosol_forcing.lower()[0] == 'e':
                if diagnostics!='AR6':
                    F[:,iF_aer] = F_aerosol
                    ariaci[:] = np.nan
                else:
                    raise ValueError('AR6 diagnostics not compatible with ' +
                'externally forced aerosols')
            else:
                raise ValueError("aerosol_forcing should be one of 'stevens', " +
                  "aerocom, aerocom+ghan, aerocom+stevens or external")
        else:
            if diagnostics!='AR6':
                F[:,iF_aer] = F_aerosol
                ariaci[:] = np.nan
            else:
                raise ValueError('AR6 diagnostics not compatible with ' +
            'externally forced aerosols')

        # Black carbon on snow - no feedback dependence
        # v1.5 update: don't require emissions driven runs here for GHGs
        # because SLCFs can still be given as emissions with GHGs as
        # concentrations
        if type(emissions) is not bool:
           if bcsnow_forcing.lower()[0]=='e':
               F[:,iF_bcsn] = bc_snow.linear(emissions-E_pi, F_ref=F_ref_BC,
                   E_ref=E_ref_BC)
           else:
               F[:,iF_bcsn] = F_bcsnow
        else:
            F[:,iF_bcsn] = F_bcsnow

        # Land use change - either use a scaling with cumulative CO2 emissions
        # or an external time series
        # v1.5 update: don't require emissions driven runs here for GHGs
        # because SLCFs can still be given as emissions with GHGs as
        # concentrations
        if type(emissions) is not bool:
            if landuse_forcing.lower()[0]=='c':
                F[:,iF_luch] = landuse.cumulative(emissions-E_pi, aCO2land=aCO2land)
            elif landuse_forcing.lower()[0]=='e':
                F[:,iF_luch] = F_landuse
            else:
                raise ValueError(
                "landuse_forcing should be one of 'co2' or 'external'")
        else:
            F[:,iF_luch] = F_landuse

        # Volcanic and solar copied straight to the output arrays
        F[:,iF_volc] = F_volcanic
        F[:,iF_solr] = F_solar

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
        if temperature_function=='Millar':
            T_j[0,:] = (q[0,:]/d)*(np.sum(F[0,:]))
            T[0] = np.sum(T_j[0,:])
        else:
            T_j[0,:,:], heatflux[0], del_ohc, lambda_eff[0] = forcing_to_temperature(
                T_j[0,:,:],
                np.sum(F[0,:]),
                np.sum(F[0,:]),
                lambda_global=lambda_global,
                ocean_heat_capacity=ocean_heat_capacity,
                ocean_heat_exchange=ocean_heat_exchange,
                deep_ocean_efficacy=deep_ocean_efficacy,
                dt=1
            )
            T[0] = np.sum(T_j[0,:,:], axis=1)[0]
            ohc[0] = ohc[0] + del_ohc


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

                if gir_carbon_cycle:
                    time_scale_sf = calculate_alpha(
                        cumulative_emissions[t-1],
                        airborne_emissions[t-1],
                        T[t-1],
                        r0, rc, rt, g0, g1)
                    C[t,0], R_i[t,:], airborne_emissions[t] = step_concentration(
                        R_i[t-1,:] + oxidised_CH4,
                        np.sum(emissions[t-1,1:3]),
                        time_scale_sf,
                        a,
                        tau,
                        C_pi[0],
                    )
                else:
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
                F[t,0:3] = ghg(C[t,0:3], C_pi[0:3], F2x=F2x, scale_F2x=scale_F2x)
                if diagnostics=='AR6':
                    F[t,3:31] = minor_gases(C[t,3:], C_pi[3:])
                else:
                    F[t,3] = np.sum(minor_gases(C[t,3:], C_pi[3:]))

                if tropO3_forcing[0].lower()=='s':
                    F[t,iF_tro3] = ozone_tr.stevenson(emissions[t,:],
                      C[t,1],
                      T=T[t-1], 
                      feedback=useTropO3TFeedback,
                      fix_pre1850_RCP=fixPre1850RCP,
                      PI=pi_tro3)
                elif tropO3_forcing[0].lower()=='c':
                    F[t,iF_tro3] = ozone_tr.cmip6_stevenson(emissions[t,:], C[t,1],
                      T=np.sum(T_j[t,:]),
                      feedback=useTropO3TFeedback,
                      PI=np.array([C_pi[1],E_pi[6],E_pi[7],E_pi[8]]),
                      beta=b_tro3)
                elif tropO3_forcing[0].lower()=='r':
                    F[t,iF_tro3] = ozone_tr.regress(emissions[t,:]-E_pi, beta=b_tro3)
                elif tropO3_forcing[0].lower()=='t':
                    F[t,iF_tro3] = ozone.thornhill_skeie(
                        emissions=emissions[t,:],
                        concentrations=C[t,:],
                        temperature=T[t-1],
                        feedback=ozone_feedback,
                        beta=b_tro3,
                        emissions_pi=E_pi,
                        concentrations_pi=C_pi,
                    )
                else:
                    F[t,iF_tro3] = F_tropO3[t]

                if tropO3_forcing[0].lower()!='t':
                    F[t,iF_sto3] = ozone_st.magicc(C[t,15:], C_pi[15:])
                F[t,iF_ch4h] = h2o_st.linear(F[t,1], ratio=stwv_from_ch4)

                # multiply by scale factors
                F[t,:] = F[t,:] * scale[t,:]

                # 3. Temperature
                # Update the thermal response boxes
                if temperature_function=='Millar':
                    T_j[t,:] = forcing_to_temperature(
                      T_j[t-1,:], q[t,:], d, F[t,:], e=efficacy)
                    T[t] = np.sum(T_j[t,:])
                else:
                    T_j[t,:,:], heatflux[t], del_ohc, lambda_eff[t] = forcing_to_temperature(
                        T_j[t-1,:,:],
                        np.sum(F[t-1,:]),
                        np.sum(F[t,:]),
                        lambda_global=lambda_global,
                        ocean_heat_capacity=ocean_heat_capacity,
                        ocean_heat_exchange=ocean_heat_exchange,
                        deep_ocean_efficacy=deep_ocean_efficacy,
                        dt=1
                    )
                    T[t] = np.sum(T_j[t,:,:], axis=1)[0]
                    ohc[t] = ohc[t-1] + del_ohc

            else:
                if t == 1:
                    time_scale_sf = 0.16
                if gir_carbon_cycle:
                    time_scale_sf = calculate_alpha(
                        cumulative_emissions[t-1],
                        airborne_emissions[t-1],
                        T[t-1],
                        r0, rc, rt, g0, g1)
                    C[t,0], R_i[t,:], airborne_emissions[t] = step_concentration(
                        R_i[t-1,:] + oxidised_CH4,
                        emissions[t-1],
                        time_scale_sf,
                        a,
                        tau,
                        C_pi[0],
                    )
                else:
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

                if temperature_function=='Millar':
                    T_j[t,:] = forcing_to_temperature(
                      T_j[t-1,:], q[t,:], d, F[t,:])
                    T[t] = np.sum(T_j[t,:])
                else:
                    T_j[t,:,:], heatflux[t], del_ohc, lambda_eff[t] = forcing_to_temperature(
                        T_j[t-1,:,:],
                        np.sum(F[t-1,:]),
                        np.sum(F[t,:]),
                        lambda_global=lambda_global,
                        ocean_heat_capacity=ocean_heat_capacity,
                        ocean_heat_exchange=ocean_heat_exchange,
                        deep_ocean_efficacy=deep_ocean_efficacy,
                        dt=1
                    )
                    T[t] = np.sum(T_j[t,:,:], axis=1)[0]
                    ohc[t] = ohc[t-1] + del_ohc


        else:

            if useMultigas:
                F[t,0:3] = ghg(C[t,0:3], C_pi[0:3], F2x=F2x)
                if diagnostics=='AR6':
                    F[t,3:31] = minor_gases(C[t,3:], C_pi[3:])
                else:
                    F[t,3] = np.sum(minor_gases(C[t,3:], C_pi[3:]))
                if type(emissions) is not bool:
                    if tropO3_forcing[0].lower()=='s':
                        F[t,iF_tro3] = ozone_tr.stevenson(emissions[t,:]-E_pi,
                          C[t,1],
                          T=T[t-1],
                          feedback=useTropO3TFeedback,
                          fix_pre1850_RCP=fixPre1850RCP)
                    elif tropO3_forcing[0].lower()=='c':
                        F[t,iF_tro3] = ozone_tr.cmip6_stevenson(emissions[t,:], C[t,1],
                          T=np.sum(T_j[t,:]),
                          feedback=useTropO3TFeedback,
                          PI=np.array([C_pi[1],E_pi[6],E_pi[7],E_pi[8]]),
                          beta=b_tro3)
                    elif tropO3_forcing[0].lower()=='r':
                        F[t,iF_tro3] = ozone_tr.regress(emissions[t,:]-E_pi, beta=b_tro3)
                    elif tropO3_forcing[0].lower()=='t':
                        F[t,iF_tro3] = ozone.thornhill_skeie(
                            emissions=emissions[t,:],
                            concentrations=C[t,:],
                            temperature=T[t-1],
                            feedback=ozone_feedback,
                            beta=b_tro3,
                            emissions_pi=E_pi,
                            concentrations_pi=C_pi,
                        )
                    else:
                        F[t,iF_tro3] = F_tropO3[t]
                else:
                    F[t,iF_tro3] = F_tropO3[t]

                if tropO3_forcing[0].lower()!='t':
                    F[t,iF_sto3] = ozone_st.magicc(C[t,15:], C_pi[15:])
                F[t,iF_ch4h] = h2o_st.linear(F[t,1], ratio=stwv_from_ch4)

                # multiply by scale factors
                F[t,:] = F[t,:] * scale[t,:]

                # 3. Temperature
                # Update the thermal response boxes
                if temperature_function=='Millar':
                    T_j[t,:] = forcing_to_temperature(
                      T_j[t-1,:], q[t,:], d, F[t,:], e=efficacy)
                    T[t] = np.sum(T_j[t,:])
                else:
                    T_j[t,:,:], heatflux[t], del_ohc, lambda_eff[t] = forcing_to_temperature(
                        T_j[t-1,:,:],
                        np.sum(F[t-1,:]),
                        np.sum(F[t,:]),
                        lambda_global=lambda_global,
                        ocean_heat_capacity=ocean_heat_capacity,
                        ocean_heat_exchange=ocean_heat_exchange,
                        deep_ocean_efficacy=deep_ocean_efficacy,
                        dt=1
                    )
                    T[t] = np.sum(T_j[t,:,:], axis=1)[0]
                    ohc[t] = ohc[t-1] + del_ohc

            else:
                if np.isscalar(other_rf):
                    F[t,0] = co2_log(C[t,0], C_pi[0], F2x) + other_rf
                else:
                    F[t,0] = co2_log(C[t,0], C_pi[0], F2x) + other_rf[t]

                F[t,0] = F[t,0] * scale[t]

                if temperature_function=='Millar':
                    T_j[t,:] = forcing_to_temperature(T_j[t-1,:], q[t,:], d, F[t,:])
                    T[t] = np.sum(T_j[t,:])
                else:
                    T_j[t,:,:], heatflux[t], del_ohc, lambda_eff[t] = forcing_to_temperature(
                        T_j[t-1,:,:],
                        np.sum(F[t-1], axis=1),
                        np.sum(F[t], axis=1),
                        lambda_global=lambda_global,
                        ocean_heat_capacity=ocean_heat_capacity,
                        ocean_heat_exchange=ocean_heat_exchange,
                        deep_ocean_efficacy=deep_ocean_efficacy,
                        dt=1
                    )
                    T[t] = np.sum(T_j[t,:,:], axis=1)[0]
                    ohc[t] = ohc[t-1] + del_ohc

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

    if ariaci_out:
        if temperature_function=='Geoffroy':
            return C, F, T, ariaci, lambda_eff, ohc, heatflux
        else:
            return C, F, T, ariaci
    else:
        if temperature_function=='Geoffroy':
            if gir_carbon_cycle:
                return C, F, T, lambda_eff, ohc, heatflux, airborne_emissions/cumulative_emissions
            else:
                return C, F, T, lambda_eff, ohc, heatflux
        else:
            return C, F, T

