import numpy as np

class Emissions:
    def __init__(self, emissions):
        self.emissions = emissions
        self.year      = emissions[:,0]
        self.co2_fossil= emissions[:,1]
        self.co2_land  = emissions[:,2]
        self.co2       = np.sum(emissions[:,1:3],axis=1)
        self.ch4       = emissions[:,3]
        self.n2o       = emissions[:,4]
        self.sox       = emissions[:,5]
        self.co        = emissions[:,6]
        self.nmvoc     = emissions[:,7]
        self.nox       = emissions[:,8]
        self.bc        = emissions[:,9]
        self.oc        = emissions[:,10]
        self.nh3       = emissions[:,11]
        self.cf4       = emissions[:,12]
        self.c2f6      = emissions[:,13]
        self.c6f14     = emissions[:,14]
        self.hfc23     = emissions[:,15]
        self.hfc32     = emissions[:,16]
        self.hfc43_10  = emissions[:,17]
        self.hfc125    = emissions[:,18]
        self.hfc134a   = emissions[:,19]
        self.hfc143a   = emissions[:,20]
        self.hfc227ea  = emissions[:,21]
        self.hfc245fa  = emissions[:,22]
        self.sf6       = emissions[:,23]
        self.cfc11     = emissions[:,24]
        self.cfc12     = emissions[:,25]
        self.cfc113    = emissions[:,26]
        self.cfc114    = emissions[:,27]
        self.cfc115    = emissions[:,28]
        self.carb_tet  = emissions[:,29]
        self.mcf       = emissions[:,30]
        self.hcfc22    = emissions[:,31]
        self.hcfc141b  = emissions[:,32]
        self.hcfc142b  = emissions[:,33]
        self.halon1211 = emissions[:,34]
        self.halon1202 = emissions[:,35]
        self.halon1301 = emissions[:,36]
        self.halon2402 = emissions[:,37]
        self.ch3br     = emissions[:,38]
        self.ch3cl     = emissions[:,39]


class Concentrations:
    def __init__(self, concentrations):
        self.concentrations = concentrations
        self.gas_indices= np.concatenate(([3,4,5], np.arange(8,36)))
        self.gases      = concentrations[:,self.gas_indices]
        self.year       = concentrations[:,0]
        self.co2eq      = concentrations[:,1]
        self.kyotoco2eq = concentrations[:,2]
        self.co2        = concentrations[:,3]
        self.ch4        = concentrations[:,4]
        self.n2o        = concentrations[:,5]
        self.fgassum    = concentrations[:,6]
        self.mhalosum   = concentrations[:,7]
        self.cf4        = concentrations[:,8]
        self.c2f6       = concentrations[:,9]
        self.c6f14      = concentrations[:,10]
        self.hfc23      = concentrations[:,11]
        self.hfc32      = concentrations[:,12]
        self.hfc43_10   = concentrations[:,13]
        self.hfc125     = concentrations[:,14]
        self.hfc134a    = concentrations[:,15]
        self.hfc143a    = concentrations[:,16]
        self.hfc227ea   = concentrations[:,17]
        self.hfc245fa   = concentrations[:,18]
        self.sf6        = concentrations[:,19]
        self.cfc11      = concentrations[:,20]
        self.cfc12      = concentrations[:,21]
        self.cfc113     = concentrations[:,22]
        self.cfc114     = concentrations[:,23]
        self.cfc115     = concentrations[:,24]
        self.carb_tet   = concentrations[:,25]
        self.mcf        = concentrations[:,26]
        self.hcfc22     = concentrations[:,27]
        self.hcfc141b   = concentrations[:,28]
        self.hcfc142b   = concentrations[:,29]
        self.halon1211  = concentrations[:,30]
        self.halon1202  = concentrations[:,31]
        self.halon1301  = concentrations[:,32]
        self.halon2402  = concentrations[:,33]
        self.ch3br      = concentrations[:,34]
        self.ch3cl      = concentrations[:,35]
    
    
class Forcing:
    def __init__(self, forcing):
        self.forcing   = forcing
        self.year      = forcing[:,0]
        self.total     = forcing[:,1]
        self.volcanic  = forcing[:,2]
        self.solar     = forcing[:,3]
        self.ghg       = forcing[:,5]
        self.co2       = forcing[:,8]
        self.ch4       = forcing[:,9]
        self.n2o       = forcing[:,10]
        self.fgas      = forcing[:,11]
        self.halo      = forcing[:,12]
        self.aero      = forcing[:,41]
        self.cloud     = forcing[:,48]
        self.strato3   = forcing[:,49]
        self.tropo3    = forcing[:,50]
        self.stwv      = forcing[:,51]
        self.dust      = forcing[:,47]
        self.landuse   = forcing[:,52]
        self.bcsnow    = forcing[:,53]
        
        
def load_emissions_data(filepath, skiprows):
    emissions = np.loadtxt(filepath, skiprows=skiprows, delimiter=',')
    out = Emissions(emissions)
    return out

def load_concentrations_data(filepath, skiprows):
    concentrations = np.loadtxt(filepath, skiprows=skiprows, delimiter=',')
    out = Concentrations(concentrations)
    return out

def load_forcing_data(filepath, skiprows):
    forcing = np.loadtxt(filepath, skiprows=skiprows, delimiter=',')
    out = Forcing(forcing)
    return out