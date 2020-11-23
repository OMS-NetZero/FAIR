from ..ancil.one_point_six import natural, cmip6_volcanic, cmip6_solar, historical_scaling
from ..constants.one_point_six import molwt, lifetime, radeff, cl_atoms, br_atoms, fracrel
from ..constants.one_point_six.general import M_ATMOS, ppm_gtc, EARTH_RADIUS, SECONDS_PER_YEAR
from ..defaults.one_point_six import carbon, thermal
import numpy as np
from .. import emissions_driven

#TODO: Add in support for 'scale'

def _forward_interpretter(
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

    if useMultigas:
        if emissions_driven:
            if type(emissions) is not np.ndarray or emissions.shape[1] != 40:
                raise ValueError(
                  "emissions timeseries should be a nt x 40 numpy array")
            nt = emissions.shape[0]
        else:
            #TODO: Implement Concentration Driven Mode
            raise ValueError("Concentration Driven Mode is Currently Unsupported")

        if type(fossilCH4_frac) is np.ndarray or fossilCH4_frac is not 0:
            Raise ValueError("Oxidisation of CH4 is not currently supported in FaIR 2.0")
        if type(lifetimes) is np.ndarray:
            if len(lifetimes)!=31:
                raise ValueError(
                  "custom GHG lifetime array must have 31 elements"
                )
        else:
            lifetimes = lifetime.aslist
        # Select the desired GHG forcing relationship and populate 
        # stratospheric water vapour from methane scale factor if not specified
        # by user
        
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

        if scale is not None or scaleHistoricalAR5 is not None:
            #TODO: Implement Scale
            Raise ValueError("Scale is not supported in FaIR 2.0")

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

        if emissions_driven:
            if type(emissions) is np.ndarray:
                if emissions.ndim != 1:
                    raise ValueError(
                      "In CO2-only mode, emissions should be a 1D array")
                nt = emissions.shape[0]
            elif type(other_rf) is np.ndarray:
                if other_rf.ndim != 1:
                    raise ValueError(
                      "In CO2-only mode, other_rf should be a 1D array")
                nt = other_rf.shape[0]
                emissions = np.zeros(nt)
            else:
                raise ValueError(
                  "Neither emissions or other_rf is defined as a timeseries")

        else:
            raise ValueError("Concentration Driven Mode is Currently Unsupported")

        if scale is not None or scaleHistoricalAR5 is not None:
            #TODO: Implement Scale
            Raise ValueError("Scale is not supported in FaIR 2.0")

    # Check a and tau are same size
    if a.ndim != 1:
        raise ValueError("a should be a 1D array")
    if tau.ndim != 1:
        raise ValueError("tau should be a 1D array")
    if len(a) != len(tau):
        raise ValueError("a and tau should be the same size")
    if not np.isclose(np.sum(a), 1.0):
        raise ValueError("a should sum to one")

    # Initialise GIR gas cycle parameters
    if gir_carbon_cycle:
        #TODO: Shift into either MultiGas or SingleGas modes seperately
        ext_forcing_np = np.zeros(nt)
        if np.isscalar(other_rf):
            ext_forcing_np += nt*[other_rf]
        else:
            ext_forcing_np += other_rf
        

        inp_ar = np.zeros((nt,ngas))
        if useMultigas:
            emis_i = {
                'time' : 0,
                'CO2' : 1,
                'CO2 Land Use' : 2,
                'CH4' : 3,
                'N2O' : 4,
                'Sox' : 5,
                'CO' : 6,
                'NMVOC' : 7,
                'NOx' : 8,
                'BC' : 9,
                'OC' : 10,
                'NH3' : 11,
                'CF4' : 12,
                'C2F6' : 13,
                'C6F14' : 14,
                'HFC23' : 15,
                'HFC32' : 16,
                'HFC43_10MEE' : 17,
                'HFC125' : 18,
                'HFC134A' : 19,
                'HFC143A' : 20,
                'HFC227EA' : 21,
                'HFC245FA' : 22,
                'SF6' : 23,
                'CFC11' : 24,
                'CFC12' : 25,
                'CFC113' : 26,
                'CFC114' : 27,
                'CFC115' : 28,
                'CCL4' : 29,
                'CH3CCL3' : 30,
                'HCFC22' : 31,
                'HCFC141B' : 32,
                'HCFC142B' : 33,
                'HALON1211' : 34,
                'HALON1202' : 35,
                'HALON1301' : 36,
                'HALON2402' : 37,
                'CH3BR' : 38,
                'CH3CL' : 39,
            }

            inp_ar_i = {
                'CO2' : 0,
                'CH4' : 1,
                'N2O' : 2,
                'Sox' : 3,
                'CO' : 4,
                'NMVOC' : 5,
                'NOx' : 6,
                'NOx|avi' : 7,
                'BC' : 8,
                'OC' : 9,
                'NH3' : 10,
                'CF4' : 11,
                'C2F6' : 12,
                'C6F14' : 13,
                'HFC23' : 14,
                'HFC32' : 15,
                'HFC43_10MEE' : 16,
                'HFC125' : 17,
                'HFC134A' : 18,
                'HFC143A' : 19,
                'HFC227EA' : 20,
                'HFC245FA' : 21,
                'SF6' : 22,
                'CFC11' : 23,
                'CFC12' : 24,
                'CFC113' : 25,
                'CFC114' : 26,
                'CFC115' : 27,
                'CCL4' : 28,
                'CH3CCL3' : 29,
                'HCFC22' : 30,
                'HCFC141B' : 31,
                'HCFC142B' : 32,
                'HALON1211' : 33,
                'HALON1202' : 34,
                'HALON1301' : 35,
                'HALON2402' : 36,
                'CH3BR' : 37,
                'CH3CL' : 38,
                'ghan_forcing' : 39,
                'cumulative_co2_land' : 40,
            }

            forc_i = {
                'CO2' : 0,
                'CH4' : 1,
                'N2O' : 2,
                'SOx' : 3,
                'CO' : 4,
                'NMVOC' : 5,
                'NOx' : 6,
                'NOx|avi' : 7,
                'BC' : 8,
                'OC' : 9,
                'NH3' : 10,
                'CF4' : 11,
                'C2F6' : 12,
                'C6F14' : 13,
                'HFC23' : 14,
                'HFC32' : 15,
                'HFC43_10MEE' : 16,
                'HFC125' : 17,
                'HFC134A' : 18,
                'HFC143A' : 19,
                'HFC227EA' : 20,
                'HFC245FA' : 21,
                'SF6' : 22,
                'CFC11' : 23,
                'CFC12' : 24,
                'CFC113' : 25,
                'CFC114' : 26,
                'CFC115' : 27,
                'CCL4' : 28,
                'CH3CCL3' : 29,
                'HCFC22' : 30,
                'HCFC141B' : 31,
                'HCFC142B' : 32,
                'HALON1211' : 33,
                'HALON1202' : 34,
                'HALON1301' : 35,
                'HALON2402' : 36,
                'CH3BR' : 37,
                'CH3CL' : 38,
                'ghan_forcing' : 39,
                'cumulative_co2_land' : 40,
                'CH4|Trop_o3' : 41,
                'CH4|Strat_h2o' : 42,
                'CO|Trop_o3' : 43,
                'NOx|Trop_o3' : 44,
                'NMVOC|Trop_o3' : 45,
                'BC|BC_on_snow' : 46,
                'NOx_avi|contrails' : 47,
                'SOx|aci' : 48,
                'BC|aci' : 49,
                'OC|aci' : 50,
                'CCL4|strat_o3' : 51,
                'CFC113|strat_o3' : 52,
                'CFC114|strat_o3' : 53,
                'CFC115|strat_o3' : 54,
                'CFC11|strat_o3' : 55,
                'CFC12|strat_o3' : 56,
                'CH3CCL3|strat_o3' : 57,
                'HALON1211|strat_o3' : 58,
                'HALON1202|strat_o3' : 59,
                'HALON1301|strat_o3' : 60,
                'HALON2402|strat_o3' : 61,
                'HCFC141B|strat_o3' : 62,
                'HCFC142B|strat_o3' : 63,
                'HCFC22|strat_o3' : 64,
                'CH3BR|strat_o3' : 65,
                'CH3CL|strat_o3' : 66,
            }

            #TODO: Oxidisation of CH4 is currently unsupported by FaIR 2.0
            a1_np, a2_np, a3_np, a4_np = np.array([len(inp_ar_i)*[1],*3*[len(inp_ar_i)*[0]]])
            a1_np[inp_ar_i['CO2']] = a[0]
            a2_np[inp_ar_i['CO2']] = a[1]
            a3_np[inp_ar_i['CO2']] = a[2]
            a4_np[inp_ar_i['CO2']] = a[3]

            tau1_np, tau2_np, tau3_np, tau4_np = np.ones((4,len(inp_ar_i)))
            tau1_np[inp_ar_i['CO2']] = tau[0]
            tau1_np[[inp_ar_i['CH4'],inp_ar_i['N2O']]] = lifetimes[[1,2]]
            tau1_np[inp_ar_i['CF4']:inp_ar_i['CH3CL']+1] = lifetimes[3:]
            tau1_np[[inp_ar_i['ghan_forcing'],inp_ar_i['cumulative_co2_land']]] = 1
            tau2_np[inp_ar_i['CO2']] = tau[1]
            tau3_np[inp_ar_i['CO2']] = tau[2]
            tau4_np[inp_ar_i['CO2']] = tau[3]

            r0_np, rC_np, rT_np, rA_np = np.zeros((4,len(inp_ar_i)))
            r0_np[inp_ar_i['CO2']] = r0
            r0_np[inp_ar_i['CO2']+1:] = tau1_np[1:]*(1-np.exp(-100/tau1_np[1:]))
            rC_np[inp_ar_i['CO2']] = rC
            rT_np[inp_ar_i['CO2']] = rT
            rA_np[inp_ar_i['CO2']] = rA

            inp_ar = np.zeros((len(inp_ar_i), nt))
            inp_ar[inp_ar_i['CO2'],:] = emissions[:,emis_i['CO2']] + emissions[:,emis_i['CO2 Land Use']] 
            inp_ar[inp_ar_i['CH4'],:] = emissions[:,emis_i['CH4']] + natural[:,0]
            inp_ar[inp_ar_i['N2O'],:] = emissions[:,emis_i['N2O']] + natural[:,1]
            inp_ar[inp_ar_i['SOx']:inp_ar_i['NOx']+1,:] = emissions[:,emis_i['SOx']:emis_i['NOx']+1]
            inp_ar[inp_ar_i['NOx|avi'],:] = inp_ar[:,inp_ar_i['NOx']]*aviNOx_frac
            inp_ar[inp_ar_i['BC']:inp_ar_i['CH3CL']+1,:] = emissions[:,emis_i['BC']:]

            PI_conc_np = np.zeros(len(inp_ar_i))
            PI_conc_np[inp_ar_i['CO2']:inp_ar_i['N2O']+1] = C_pi[0:3]
            PI_conc_np[inp_ar_i['C2F6']:inp_ar_i['CH3CL']+1,:] = C_pi[3:]

            F2x_etminan = (-2.4e-7*PI_conc_np[inp_ar_i['CO2']**2 + 7.2e-4*PI_conc_np[inp_ar_i['CO2'] - 2.1e-4*PI_conc_np[inp_ar_i['N2O']] + 5.36) * np.log(2)

            f1_np = np.zeroes(len(forc_i))
            f2_np = np.zeroes(len(forc_i))
            f3_np = np.zeroes(len(forc_i))
            f2_np[forc_i['C2F6']:forc_i['CH3CL']+1] = radeff.aslist[4:]
            if ghg_forcing.lower()=="etminan":
                if stwv_from_ch4==None: stwv_from_ch4=0.12
                f1_np[forc_i['CO2']] = (-2.1e-4 * PI_conc_np[inp_ar_i['N2O']] + 5.36)*(F2x/F2x_etminan if scale_F2x else 1)
                f3_np[forc_i['CH4']] = -1/3e-6 * PI_conc_np[inp_ar_i['CH4']] - 8.2e-6 * PI_conc_np[inp_ar_i['N2O']] + 0.043
                f3_np[forc_i['N2O']] =  -8e-6 * PI_conc_np[inp_ar_i['CO2']] + 4.2e-6 * PI_conc_np[inp_ar_i['N2O']] - 4.9e-6 * PI_conc_np[inp_ar_i['CH4']] + 0.117
            elif ghg_forcing.lower()=="myhre":
                if stwv_from_ch4==None: stwv_from_ch4=0.15
                f1_np[forc_i['CO2']] = F2x/np.log(2)
                f3_np[forc_i['CH4']] = 0.036
                f3_np[forc_i['N2O']] = 0.12
            elif ghg_forcing.lower()=="meinshausen":
                if stwv_from_ch4==None: stwv_from_ch4=0.12
                f1_np[forc_i['CO2']] = (5.2488 - 0.0021492 * sqrt(PI_conc_np[inp_ar_i['N2O']]))*(F2x/F2x_etminan if scale_F2x else 1)
                f2_np[forc_i['CH4']] = -8.9603e-05
                f2_np[forc_i['N2O']] = 0.00025455
                f3_np[forc_i['CH4']] = 8.9603e-05*sqrt(PI_conc_np[inp_ar_i['CH4']]) - 0.00012462*sqrt(PI_conc_np[inp_ar_i['N2O']]) + 0.045194
                f3_np[forc_i['N2O']] = -0.00034197*sqrt(PI_conc_np[inp_ar_i['CO2']) + -0.00024357*sqrt(PI_conc_np[inp_ar_i['CH4']]) - 0.00025455*sqrt(PI_conc_np[inp_ar_i['N2O']]) + 0.12173)
            else:
                raise ValueError(
                "ghg_forcing should be 'etminan' (default), 'meinshausen' or 'myhre'")
            
            #TODO: Ensure Correct Indices for these

            Cl = np.array(cl_atoms.aslist)
            Br = np.array(br_atoms.aslist)
            FC = np.array(fracrel.aslist)
            f2[ 
                [
                    forc_i['CFC11|strat_o3'],
                    forc_i['CFC12|strat_o3'],
                    forc_i['CFC113|strat_o3'],
                    forc_i['CFC114|strat_o3'],
                    forc_i['CFC115|strat_o3'],
                    forc_i['CCL4|strat_o3'],
                    forc_i['CH3CCL3|strat_o3'],
                    forc_i['HCFC22|strat_o3'],
                    forc_i['HCFC141B|strat_o3'],
                    forc_i['HCFC142B|strat_o3'],
                    forc_i['HALON1211|strat_o3'],
                    forc_i['HALON1202|strat_o3'],
                    forc_i['HALON1301|strat_o3'],
                    forc_i['HALON2402|strat_o3'],
                    forc_i['CH3BR|strat_o3'],
                    forc_i['CH3CL|strat_o3'],        
                ]
            ] = 1000(45*Br + Cl/FC.CFC11)*FC

            #TODO:CH4 Strat H2O Done
            f1[
                forc_i['CH4|Strat_h2o']
                ] = f1[
                        forc_i['CH4']
                    ]*stwv_from_ch4
            f2[
                forc_i['CH4|Strat_h2o']
                ] = f2[
                        forc_i['CH4']
                    ]*stwv_from_ch4
            f3[
                forc_i['CH4|Strat_h2o']
                ] = f3[
                        forc_i['CH4']
                    ]*stwv_from_ch4

            #TODO: Stevenson and CMIP6 Stevenson Functions DONE
            if useTropO3TFeedback:
                Raise ValueError('Temperature Feedback is not an option with FaIR 2.0')
            if fixPre1850RCP:
                Raise ValueError('Fixing Pre 1850 RCP is not implemented with FaIR 2.0')
            
            if tropO3_forcing[0].lower() == 's' or tropO3_forcing[0].lower() == 'c'
                #Set PI_CONC[CO] to PI *Emissions*
                PI_conc_np[
                    [
                        inp_ar_i['CO'],
                        inp_ar_i['NMVOC'],
                        inp_ar_i['NOx']
                    ]
                ] = E_pi[
                        [
                            6,
                            7,
                            8]
                    ]
                if tropO3_forcing[0].lower()=='s':
                    f2_np[
                        [
                            forc_i['CH4|Trop_o3'],
                            forc_i['CO|Trop_o3'],
                            forc_i['NMVOC|Trop_o3'],
                            forc_i['NOx|Trop_o3']
                        ]
                    ] = [
                            0.166/960,
                            0.058/681.8,
                            0.035/155.84,
                            0.119/61.16
                        ]
                    #Adjust inp_ar of NOx according to molwt.NO/molwt.N
                    inp_ar[
                        inp_ar_i['NOx']
                    ] *= (molwt.NO/molwt.N)
                else:
                    f2_np[
                        [
                            forc_i['CH4|Trop_o3'],
                            forc_i['CO|Trop_o3'],
                            forc_i['NMVOC|Trop_o3'],
                            forc_i['NOx|Trop_o3']
                        ]
                    ] = b_tro3
            elif tropO3_forcing[0].lower()=='r':
                Raise ValueError('Regress Tropospheric Ozone Function is deprecated')
            else:
                ext_forcing_np += F_tropO3
            
            #TODO: DONE Code in contrail forcing from_aviNOx and from_fuel
            if contrail_forcing.lower()[0]=='n':   # from NOx emissions
                f2_np[
                    forc_i[
                        'NOx_avi|contrails'
                    ]
                ] = (F_ref_aviNOx/E_ref_aviNOx)*(molwt.NO2/molwt.N)
            elif contrail_forcing.lower()[0]=='f': # from kerosene production
                F_kerosene_contrails = kerosene_supply*(0.0448/236.063)
                ext_forcing_np += F_kerosene_contrails
            elif contrail_forcing.lower()[0]=='e': # external forcing timeseries
                ext_forcing_np += F_contrails
            else:
                raise ValueError("contrails must be one of 'NOx' (estimated "+
                 "from NOx emissions), 'fuel' (estimated from annual jet fuel "+
                 "supplied) or 'external' (an external forcing time series).")

            #TODO: DONE Code in Aerosol Direct and Cloud Forcings
            if aerosol_forcing.lower()=='stevens' or ('aercom' in aerosol_forcing.lower() and not 'ghan' in aerosol_forcing.lower() and 'stevens' in aerosol_forcing.lower()):
                alpha, beta, natural_SOx = stevens_params
                factor = 1
                if ref_isSO2:
                    factor = molwt.SO2/molwt.S
                inp_ar[
                    inp_ar_i['SOx']
                ] += natural_SOx
                inp_ar[
                    inp_ar_i['SOx']
                ] *= factor
                PI_conc_np[
                    inp_ar_i['SOx']
                ] = E_pi[5] + natural_SOx
                PI_conc_np[
                    inp_ar_i['SOx']
                ] *= factor
                f2_np[
                    forc_i['SOx']
                ] = -alpha
                #Set f1[SOx|ari] to -beta
                f1[
                    forc_i['SOx|aci']
                ] = -beta
            if 'aerocom' in aerosol_forcing.lower():
                #These are amalgamated into one for non-AR6
                #Set f2[SOx, CO, NMVOC, NOx, BC, OC, NH3] to beta
                f2_np[
                    [
                        forc_i['SOx'],
                        forc_i['CO'],
                        forc_i['NMVOC'],
                        forc_i['NOx'],
                        forc_i['BC'],
                        forc_i['OC'],
                        forc_i['NH3'],
                    ]
                ] = b_aero
                #Set PI_conc[SOx, CO, NMVOc, NOx, BC, OC, NH3] to E_pi[5:12]
                PI_conc_np[
                    [
                        inp_ar_i['SOx'],
                        inp_ar_i['CO'],
                        inp_ar_i['NMVOC'],
                        inp_ar_i['NOx'],
                        inp_ar_i['BC'],
                        inp_ar_i['OC'],
                        inp_ar_i['NH3'],
                    ]
                ] = E_pi[5:12]
                if 'ghan' in aerosol_forcing.lower():
                    #inp_ar[ghan] = Aghan + Bghan * inp_ar[SOx] + Cghan * (inp_ar[BC] + inp_ar[OC])
                    #PI_conc_np[ghan] = Aghan + Bghan * PI_conc_np[SOx] + Cghan * (PI_conc_np[BC] + PI_conc_np[OC])
                    if 'ghan2' in aerosol_forcing.lower():
                        beta, n_so2, n_pom = ghan_params
                        Aghan = n_pom + n_so2
                        Bghan = n_pom
                        Cghan = n_so2
                        f1_np[
                            forc_i['ghan_forcing']
                        ] = -beta
                        PI_conc_np[
                        inp_ar_i['ghan_forcing']
                        ] = Aghan + Bghan * PI_conc_np[
                                inp_ar_i['SOx']
                            ] + Cghan * (
                                PI_conc_np[
                                    inp_ar_i['BC']
                                ] +
                                PI_conc_np[
                                    inp_ar_i['OC']
                                ]
                            )
                    else:
                        scale_ghan = -1.95011431
                        b_SOx = 0.01107147
                        b_POM = 0.01387492
                        F_1765 = -0.3002836449793625
                        F_2011 = -1.5236182344467388
                        Aghan = 1
                        Bghan = b_SOx
                        Cghan = b_POM
                        scaleAerosol = 1
                        if scaleAerosolAR5:
                            scaleAerosol = -0.45/(F_2011-F_1765)
                        PI_conc_np[
                        inp_ar_i['ghan_forcing']
                        ] = np.exp(F_1765/scale_ghan)
                        f1_np[
                            forc_i['ghan_forcing']
                        ] = scaleAerosol*scale_ghan
                    inp_ar[
                        inp_ar_i['ghan_forcing']
                    ] = Aghan + Bghan * inp_ar[
                            inp_ar_i['SOx']
                        ] + Cghan * (
                            inp_ar[
                                inp_ar_i['BC']
                            ] +
                            inp_ar[
                                inp_ar_i['OC']
                            ]
                        )
            elif aerosol_forcing.lower()[0] == 'e':
                if diagnostics!='AR6':
                    ext_forcing_np += F_aerosol
                else:
                    raise ValueError('AR6 diagnostics not compatible with ' +
                'externally forced aerosols')
            elif aerosol_forcing.lower()!='stevens':
                raise ValueError("aerosol_forcing should be one of 'stevens', " +
                  "aerocom, aerocom+ghan, aerocom+stevens or external")

            if bcsnow_forcing.lower()[0]=='e':
                PI_conc_np[
                    inp_ar_i['BC']
                ] = E_pi[9]
                f2_np[
                    forc_i['BC|BC_on_snow']
                ] = (F_ref_BC/E_ref_BC)
            else:
               ext_forcing_np += F_bcsnow

            #TODO: Code in iF_Landuse
            if landuse_forcing.lower()[0]=='c':
                inp_ar[
                    inp_ar_i['cumulative_co2_land']
                ] = np.cumsum(
                    emissions[:,
                        emis_i['CO2 Land Use']
                    ] - E_pi[2]
                )
                f2[
                    forc_i['cumulative_co2_land']
                ] = aCO2land
            elif landuse_forcing.lower()[0]=='e':
                ext_forcing_np += F_landuse
            else:
                raise ValueError(
                "landuse_forcing should be one of 'co2' or 'external'")

            ext_forcing_np += F_volcanic + F_solar

            mapping_ar = np.zeros(len(forc_i))
            
            for i,key in enumerate(forc_i.keys()):
                gas_string = key
                if '|' in gas_string:
                    gas_string = gas_string[0:gas_string.index('|')]
                mapping_ar[i] = forc_i[gas_string],
        else:
            f1_np[0] = F2x/np.log(2)

    else:
        #TODO: Implement Non-GIR Carbon Cycle
        raise ValueError('Non-GIR Carbon Cycle is Unsupported in FaIR 2.0')

    # import correct conversion
    if temperature_function=='Millar':
        if type(tcrecs) is np.ndarray:
            if tcrecs.ndim != 1:
                raise ValueError('Transient TCR and ECS are not supported in FaIR 2.0')
            elif len(tcres) != 2:
                raise ValueError("Constant TCR and ECS should be a 2-element array")
            k = 1.0 - (d/tcr_dbl)*(1.0 - np.exp(-tcr_dbl/d))
            q_np  = (1.0 / f2x) * (1.0/(k[0]-k[1])) * np.array([tcrecs[0]-tcrecs[1]*k[1],tcrecs[1]*k[0]-tcrecs[0]])
        else:
            q_np = q
        d_np = d
        from .temperature.millar import forcing_to_temperature
    elif temperature_function=='Geoffroy':
        #With thanks to Zeb Nicholls for doing the heavy lifting here

        C = ocean_heat_capacity[0]
        C_D = ocean_heat_capacity[1]
        lambda0 = lambda_global
        efficacy = deep_ocean_efficacy
        eta = ocean_heat_exchange

        b_pt1 = (lambda0 + efficacy * eta) / (C)
        b_pt2 = (eta) / (C_D)
        b = b_pt1 + b_pt2
        b_star = b_pt1 - b_pt2
        delta = b ** 2 - (4 * lambda0 * eta) / (C * C_D)

        dcoeff = C * C_D / (2 * lambda0 * eta)
        d1 = dcoeff * (b - delta ** 0.5)
        d2 = dcoeff * (b + delta ** 0.5)

        qdenom = C * (2 * delta ** 0.5)
        q1 = d1 * (delta ** 0.5 + b_star) / qdenom
        q2 = d2 * (delta ** 0.5 - b_star) / qdenom

        d_np = np.array([d1, d2])
        q_np = np.array([q1, q2])
    else:
        raise ValueError('temperature_function must be "Millar" or "Geoffroy"')
        
    # Set up the output timeseries variables depending on options and perform
    # basic sense checks

    timeseries = emissions[:,0]
    timestep = np.diff(timeseries)[
            [
                *range(
                    len(timeseries)-1
                ),
                -1,
            ]
        ]

    res_dict = emissions_driven._run_numpy(
        inp_ar,
        a1_np,
        a2_np,
        a3_np,
        a4_np,
        tau1_np,
        tau2_np,
        tau3_np,
        tau4_np,
        r0_np,
        rC_np,
        rT_np,
        rA_np,
        PI_conc_np,
        emis2conc_np,
        f1_np,
        f2_np,
        f3_np,
        d_np,
        q_np,
        ext_forcing_np,
        timestep,
        mapping_ar,
    )

    C_res = res_dict['C']
    T_out = res_dict['T']
    F_res = res_dict['RF']
    S_res = res_dict['S']

    if restart_out:
        Raise ValueError('Restarts are not supported in FaIR 2.0')

    


    if diagnostics == 'AR6':
        C_out = C_res[[0,1,2,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38]]
        F_out = np.array(
            [
                *F_res[[0,1,2,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38]],
                np.sum(F_res[[41, 43, 44, 45, 46]]),
                np.sum(F_res[[51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66]])
                F_res[42],
                F_res[47],
                F_res[3],
                F_res[4]+F_res[5],
                F_res[6] + F_res[10],
                F_res[8],
                F_res[9],
                F_res[48] + F_res[49] + F_res[50],
                F_res[46],
                F_landuse if landuse_forcing.lower()[0]=='e' else F_res[40],
                F_volcanic,
                F_solar
            ]
            )
        T_out = T_res

            if ariaci_out:
                ariaci = np.zeros((nt,2))
                ariaci[:,0] = F_res[3] + F_res[4] + F_res[5]+ F_res[6]+ F_res[8]+ F_res[9]+ F_res[10]
                ariaci[:,1] = F_res[48] + F_res[49] + F_res[50] + F_res[39]
                if temperature_function=='Geoffroy':
                    c_dtemp = ocean_heat_capacity[0]*[S_out[0,0], *np.diff(S_out[0,:])] + ocean_heat_capacity[1]*[S_out[1,0], *np.diff(S_out[1,:])]
                    heatflux = c_dtemp/timestep
                    ntoa_joule = 4 * np.pi * EARTH_RADIUS**2 * SECONDS_PER_YEAR
                    del_ohc  = ntoa_joule * c_dtemp
                    ohc = np.cumsum(del_ohc)

                    factor_lambda_eff = (deep_ocean_efficacy-1.0)*ocean_heat_exchange
                    lambda_eff = np.zeros(nt)
                    for i in range(nt):
                        if S[0,i] > 1e-6:
                            ratio = (S[0,i] - S[1,i])/S[1,i]
                            lambda_eff[i] = lambda_global + factor_lambda_eff*ratio
                        else:
                            lambda_eff[i] = lambda_global + factor_lambda_eff
                    return C_out, F_out, T_out, ariaci, lambda_eff, ohc, heatflux
                    

                else:
                    return C_out, F_out, T_out, ariaci
            else:
                if temperature_function=='Geoffroy':
                    Raise ValueError('Geoffroy is not currently implemented in FaIR 2.0')
                    cumulative_emissions = np.cumsum(emissions[:,1:3].sum(axis=1))
                    airborne_emissions = C_res[0]/emis2conc[0]
                    c_dtemp = ocean_heat_capacity[0]*[S_out[0,0], *np.diff(S_out[0,:])] + ocean_heat_capacity[1]*[S_out[1,0], *np.diff(S_out[1,:])]
                    heatflux = c_dtemp/timestep
                    ntoa_joule = 4 * np.pi * EARTH_RADIUS**2 * SECONDS_PER_YEAR
                    del_ohc  = ntoa_joule * c_dtemp
                    ohc = np.cumsum(del_ohc)

                    factor_lambda_eff = (deep_ocean_efficacy-1.0)*ocean_heat_exchange
                    lambda_eff = np.zeros(nt)
                    for i in range(nt):
                        if S[0,i] > 1e-6:
                            ratio = (S[0,i] - S[1,i])/S[1,i]
                            lambda_eff[i] = lambda_global + factor_lambda_eff*ratio
                        else:
                            lambda_eff[i] = lambda_global + factor_lambda_eff
                    return C_out, F_out, T_out, lambda_eff, ohc, heatflux, airborne_emissions/cumulative_emissions
                else:
                    return C_out, F_out, T_out