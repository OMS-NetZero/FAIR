# Convenience classes for loading in rcp datasets

import numpy as np


class Emissions(object):

    def __init__(self, skiprows=37, **kwargs):
        self._filename =  kwargs.get('filename')
        self.emissions = np.loadtxt(self._filename, skiprows=skiprows, delimiter=',')

    def __repr__(self):
        return '<{} with source \'{}\'>'.format(self.__class__.__name__, self._filename)

    @property
    def year(self):
        return self.emissions[:,0]

    @property
    def co2_fossil(self):
        return self.emissions[:,1]

    @property
    def co2_land(self):
        return self.emissions[:,2]

    @property
    def co2(self):
        return np.sum(self.emissions[:,1:3], axis=1)

    @property
    def ch4(self):
        return self.emissions[:,3]

    @property
    def n2o(self):
        return self.emissions[:,4]

    @property
    def sox(self):
        return self.emissions[:,5]

    @property
    def co(self):
        return self.emissions[:,6]

    @property
    def nmvoc(self):
        return self.emissions[:,7]

    @property
    def nox(self):
        return self.emissions[:,8]

    @property
    def bc(self):
        return self.emissions[:,9]

    @property
    def oc(self):
        return self.emissions[:,10]

    @property
    def nh3(self):
        return self.emissions[:,11]

    @property
    def cf4(self):
        return self.emissions[:,12]

    @property
    def c2f6(self):
        return self.emissions[:,13]

    @property
    def c6f14(self):
        return self.emissions[:,14]

    @property
    def hfc23(self):
        return self.emissions[:,15]

    @property
    def hfc32(self):
        return self.emissions[:,16]

    @property
    def hfc43_10(self):
        return self.emissions[:,17]

    @property
    def hfc125(self):
        return self.emissions[:,18]

    @property
    def hfc134a(self):
        return self.emissions[:,19]

    @property
    def hfc143a(self):
        return self.emissions[:,20]

    @property
    def hfc227ea(self):
        return self.emissions[:,21]

    @property
    def hfc245fa(self):
        return self.emissions[:,22]

    @property
    def sf6(self):
        return self.emissions[:,23]

    @property
    def cfc11(self):
        return self.emissions[:,24]

    @property
    def cfc12(self):
        return self.emissions[:,25]

    @property
    def cfc113(self):
        return self.emissions[:,26]

    @property
    def cfc114(self):
        return self.emissions[:,27]

    @property
    def cfc115(self):
        return self.emissions[:,28]

    @property
    def carb_tet(self):
        return self.emissions[:,29]

    @property
    def mcf(self):
        return self.emissions[:,30]

    @property
    def hcfc22(self):
        return self.emissions[:,31]

    @property
    def hcfc141b(self):
        return self.emissions[:,32]

    @property
    def hcfc142b(self):
        return self.emissions[:,33]

    @property
    def halon1211(self):
        return self.emissions[:,34]

    @property
    def halon1202(self):
        return self.emissions[:,35]

    @property
    def halon1301(self):
        return self.emissions[:,36]

    @property
    def halon2402(self):
        return self.emissions[:,37]

    @property
    def ch3br(self):
        return self.emissions[:,38]

    @property
    def ch3cl(self):
        return self.emissions[:,39]


class Concentrations(object):

    def __init__(self, **kwargs):
        self._filename = kwargs.get('filename')
        skiprows = kwargs.get('skiprows', 38)
        self.concentrations = np.loadtxt(self._filename, skiprows=skiprows, delimiter=',')

    def __repr__(self):
        return '<{} with source \'{}\'>'.format(self.__class__.__name__, self._filename)

    @property
    def gas_indices(self):
        return np.concatenate(([3,4,5], np.arange(8,36)))

    @property
    def gases(self):
        return self.concentrations[:,self.gas_indices]

    @property
    def year(self):
        return self.concentrations[:,0]

    @property
    def co2eq(self):
        return self.concentrations[:,1]

    @property
    def kyotoco2eq(self):
        return self.concentrations[:,2]

    @property
    def co2(self):
        return self.concentrations[:,3]

    @property
    def ch4(self):
        return self.concentrations[:,4]

    @property
    def n2o(self):
        return self.concentrations[:,5]

    @property
    def fgassum(self):
        return self.concentrations[:,6]

    @property
    def mhalosum(self):
        return self.concentrations[:,7]

    @property
    def cf4(self):
        return self.concentrations[:,8]

    @property
    def c2f6(self):
        return self.concentrations[:,9]

    @property
    def c6f14(self):
        return self.concentrations[:,10]

    @property
    def hfc23(self):
        return self.concentrations[:,11]

    @property
    def hfc32(self):
        return self.concentrations[:,12]

    @property
    def hfc43_10(self):
        return self.concentrations[:,13]

    @property
    def hfc125(self):
        return self.concentrations[:,14]

    @property
    def hfc134a(self):
        return self.concentrations[:,15]

    @property
    def hfc143a(self):
        return self.concentrations[:,16]

    @property
    def hfc227ea(self):
        return self.concentrations[:,17]

    @property
    def hfc245fa(self):
        return self.concentrations[:,18]

    @property
    def sf6(self):
        return self.concentrations[:,19]

    @property
    def cfc11(self):
        return self.concentrations[:,20]

    @property
    def cfc12(self):
        return self.concentrations[:,21]

    @property
    def cfc113(self):
        return self.concentrations[:,22]

    @property
    def cfc114(self):
        return self.concentrations[:,23]

    @property
    def cfc115(self):
        return self.concentrations[:,24]

    @property
    def carb_tet(self):
        return self.concentrations[:,25]

    @property
    def mcf(self):
        return self.concentrations[:,26]

    @property
    def hcfc22(self):
        return self.concentrations[:,27]

    @property
    def hcfc141b(self):
        return self.concentrations[:,28]

    @property
    def hcfc142b(self):
        return self.concentrations[:,29]

    @property
    def halon1211(self):
        return self.concentrations[:,30]

    @property
    def halon1202(self):
        return self.concentrations[:,31]

    @property
    def halon1301(self):
        return self.concentrations[:,32]

    @property
    def halon2402(self):
        return self.concentrations[:,33]

    @property
    def ch3br(self):
        return self.concentrations[:,34]

    @property
    def ch3cl(self):
        return self.concentrations[:,35]



class Forcing(object):
    def __init__(self, skiprows=59, **kwargs):
        self._filename = kwargs.get('filename')
        self.forcing = np.loadtxt(self._filename, skiprows=skiprows, delimiter=',')

    def __repr__(self):
        return '<{} with source \'{}\'>'.format(self.__class__.__name__, self._filename)

    @property
    def year(self):
        return self.forcing[:,0]

    @property
    def total(self):
        return self.forcing[:,1]

    @property
    def volcanic(self):
        return self.forcing[:,2]

    @property
    def solar(self):
        return self.forcing[:,3]

    @property
    def ghg(self):
        return self.forcing[:,5]

    @property
    def co2(self):
        return self.forcing[:,8]

    @property
    def ch4(self):
        return self.forcing[:,9]

    @property
    def n2o(self):
        return self.forcing[:,10]

    @property
    def fgas(self):
        return self.forcing[:,11]

    @property
    def halo(self):
        return self.forcing[:,12]

    @property
    def aero(self):
        return self.forcing[:,41]

    @property
    def cloud(self):
        return self.forcing[:,48]

    @property
    def strato3(self):
        return self.forcing[:,49]

    @property
    def tropo3(self):
        return self.forcing[:,50]

    @property
    def stwv(self):
        return self.forcing[:,51]

    @property
    def dust(self):
        return self.forcing[:,47]

    @property
    def landuse(self):
        return self.forcing[:,52]

    @property
    def bcsnow(self):
        return self.forcing[:,53]


