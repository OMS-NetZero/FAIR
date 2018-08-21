# Convenience module for loading in RCP emissions datasets
#
# Usage:
#
# import rcp45
# rcp45.Emissions.co2

import numpy as np
import os
emissions_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP45_EMISSIONS.csv')
concentrations_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP45_MIDYEAR_CONCENTRATIONS.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP45_MIDYEAR_RADFORCING.csv')
aviNOx_filename = os.path.join(
    os.path.dirname(__file__), 'data/aviNOx_fraction.csv')
fossilCH4_filename = os.path.join(
    os.path.dirname(__file__), 'data/fossilCH4_fraction.csv')

aviNOx_frac = np.loadtxt(aviNOx_filename, skiprows=5, usecols=(2,),
    delimiter=',')
fossilCH4_frac = np.loadtxt(fossilCH4_filename, skiprows=5, usecols=(2,),
    delimiter=',')

class Emissions:
    emissions = np.loadtxt(emissions_filename, skiprows=37, delimiter=',')
    year      = emissions[:,0]
    co2_fossil= emissions[:,1]
    co2_land  = emissions[:,2]
    co2       = np.sum(emissions[:,1:3],axis=1)
    ch4       = emissions[:,3]
    n2o       = emissions[:,4]
    sox       = emissions[:,5]
    co        = emissions[:,6]
    nmvoc     = emissions[:,7]
    nox       = emissions[:,8]
    bc        = emissions[:,9]
    oc        = emissions[:,10]
    nh3       = emissions[:,11]
    cf4       = emissions[:,12]
    c2f6      = emissions[:,13]
    c6f14     = emissions[:,14]
    hfc23     = emissions[:,15]
    hfc32     = emissions[:,16]
    hfc43_10  = emissions[:,17]
    hfc125    = emissions[:,18]
    hfc134a   = emissions[:,19]
    hfc143a   = emissions[:,20]
    hfc227ea  = emissions[:,21]
    hfc245fa  = emissions[:,22]
    sf6       = emissions[:,23]
    cfc11     = emissions[:,24]
    cfc12     = emissions[:,25]
    cfc113    = emissions[:,26]
    cfc114    = emissions[:,27]
    cfc115    = emissions[:,28]
    carb_tet  = emissions[:,29]
    mcf       = emissions[:,30]
    hcfc22    = emissions[:,31]
    hcfc141b  = emissions[:,32]
    hcfc142b  = emissions[:,33]
    halon1211 = emissions[:,34]
    halon1202 = emissions[:,35]
    halon1301 = emissions[:,36]
    halon2402 = emissions[:,37]
    ch3br     = emissions[:,38]
    ch3cl     = emissions[:,39]


class Concentrations:
    concentrations = np.loadtxt(concentrations_filename, skiprows=38,
      delimiter=',')
    gas_indices= np.concatenate(([3,4,5], np.arange(8,36)))
    gases      = concentrations[:,gas_indices]
    year       = concentrations[:,0]
    co2eq      = concentrations[:,1]
    kyotoco2eq = concentrations[:,2]
    co2        = concentrations[:,3]
    ch4        = concentrations[:,4]
    n2o        = concentrations[:,5]
    fgassum    = concentrations[:,6]
    mhalosum   = concentrations[:,7]
    cf4        = concentrations[:,8]
    c2f6       = concentrations[:,9]
    c6f14      = concentrations[:,10]
    hfc23      = concentrations[:,11]
    hfc32      = concentrations[:,12]
    hfc43_10   = concentrations[:,13]
    hfc125     = concentrations[:,14]
    hfc134a    = concentrations[:,15]
    hfc143a    = concentrations[:,16]
    hfc227ea   = concentrations[:,17]
    hfc245fa   = concentrations[:,18]
    sf6        = concentrations[:,19]
    cfc11      = concentrations[:,20]
    cfc12      = concentrations[:,21]
    cfc113     = concentrations[:,22]
    cfc114     = concentrations[:,23]
    cfc115     = concentrations[:,24]
    carb_tet   = concentrations[:,25]
    mcf        = concentrations[:,26]
    hcfc22     = concentrations[:,27]
    hcfc141b   = concentrations[:,28]
    hcfc142b   = concentrations[:,29]
    halon1211  = concentrations[:,30]
    halon1202  = concentrations[:,31]
    halon1301  = concentrations[:,32]
    halon2402  = concentrations[:,33]
    ch3br      = concentrations[:,34]
    ch3cl      = concentrations[:,35]

    
class Forcing:
    forcing   = np.loadtxt(forcing_filename, skiprows=59, delimiter=',')
    year      = forcing[:,0]
    total     = forcing[:,1]
    volcanic  = forcing[:,2]
    solar     = forcing[:,3]
    ghg       = forcing[:,5]
    co2       = forcing[:,8]
    ch4       = forcing[:,9]
    n2o       = forcing[:,10]
    fgas      = forcing[:,11]
    halo      = forcing[:,12]
    aero      = forcing[:,41]
    cloud     = forcing[:,48]
    strato3   = forcing[:,49]
    tropo3    = forcing[:,50]
    stwv      = forcing[:,51]
    dust      = forcing[:,47]
    landuse   = forcing[:,52]
    bcsnow    = forcing[:,53]

