# Convenience module for loading in RCP emissions datasets
#
# Usage:
#
# import rcp6
# rcp6.Emissions.co2

import numpy as np
import os
emissions_filename = os.path.join(
    os.path.dirname(__file__), 'RCP6_EMISSIONS.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'RCP6_MIDYEAR_RADFORCING.csv')
aviNOx_filename = os.path.join(
    os.path.dirname(__file__), 'aviNOx_fraction.csv')
fossilCH4_filename = os.path.join(
    os.path.dirname(__file__), 'fossilCH4_fraction.csv')

aviNOx_frac = np.loadtxt(aviNOx_filename, skiprows=5, usecols=(3,),
    delimiter=',')
fossilCH4_frac = np.loadtxt(fossilCH4_filename, skiprows=5, usecols=(3,),
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

