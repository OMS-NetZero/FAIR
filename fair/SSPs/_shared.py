import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ..constants import molwt

class Emissions:
    def __init__(self, loaded, startyear):
        self.year = loaded.columns[(startyear-1750)+7:].astype(int).to_numpy()
        self.co2 = (loaded[loaded["Variable"]=="Emissions|CO2"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.C/molwt.CO2 / 1000
        self.co2_fossil = (loaded[loaded["Variable"]=="Emissions|CO2|MAGICC Fossil and Industrial"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.C/molwt.CO2 / 1000
        self.co2_land = (loaded[loaded["Variable"]=="Emissions|CO2|MAGICC AFOLU"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.C/molwt.CO2 / 1000
        self.ch4 = (loaded[loaded["Variable"]=="Emissions|CH4"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.n2o = (loaded[loaded["Variable"]=="Emissions|N2O"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() / 1000
        self.sox = (loaded[loaded["Variable"]=="Emissions|Sulfur"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.S/molwt.SO2
        self.co = (loaded[loaded["Variable"]=="Emissions|CO"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.nmvoc = (loaded[loaded["Variable"]=="Emissions|VOC"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.nox = (loaded[loaded["Variable"]=="Emissions|NOx"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.N / molwt.NO2
        self.bc = (loaded[loaded["Variable"]=="Emissions|BC"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.oc = (loaded[loaded["Variable"]=="Emissions|OC"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.nh3 = (loaded[loaded["Variable"]=="Emissions|NH3"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.cf4 = (loaded[loaded["Variable"]=="Emissions|F-Gases|PFC|CF4"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.c2f6 = (loaded[loaded["Variable"]=="Emissions|F-Gases|PFC|C2F6"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.c6f14 = (loaded[loaded["Variable"]=="Emissions|F-Gases|PFC|C6F14"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc23 = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC23"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc32 = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC32"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc43_10 = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC4310mee"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc125 = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC125"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc134a = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC134a"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc143a = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC143a"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc227ea = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC227ea"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hfc245fa = (loaded[loaded["Variable"]=="Emissions|F-Gases|HFC|HFC245fa"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.sf6 = (loaded[loaded["Variable"]=="Emissions|F-Gases|SF6"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.cfc11 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CFC|CFC11"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.cfc12 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CFC|CFC12"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.cfc113 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CFC|CFC113"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.cfc114 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CFC|CFC114"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.cfc115 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CFC|CFC115"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.carb_tet = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CCl4"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.mcf = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CH3CCl3"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hcfc22 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|HCFC22"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hcfc141b = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|HCFC141b"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.hcfc142b = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|HCFC142b"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.halon1211 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|Halon1211"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.halon1202 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|Halon1202"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.halon1301 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|Halon1301"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.halon2402 = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|Halon2402"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.ch3br = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CH3Br"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.ch3cl = (loaded[loaded["Variable"]=="Emissions|Montreal Gases|CH3Cl"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()

        self.emissions = np.vstack((self.year, self.co2_fossil, self.co2_land, self.ch4, self.n2o, self.sox, self.co, self.nmvoc, self.nox, self.bc, self.oc, self.nh3, self.cf4, self.c2f6, self.c6f14, self.hfc23, self.hfc32, self.hfc43_10, self.hfc125, self.hfc134a, self.hfc143a, self.hfc227ea, self.hfc245fa, self.sf6, self.cfc11, self.cfc12, self.cfc113, self.cfc114, self.cfc115, self.carb_tet, self.mcf, self.hcfc22, self.hcfc141b, self.hcfc142b, self.halon1211, self.halon1202, self.halon1301, self.halon2402, self.ch3br, self.ch3cl)).T

        _nox_avi = (loaded[loaded["Variable"]=="Emissions|NOx|MAGICC Fossil and Industrial|Aircraft"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.N / molwt.NO2
        _ch4_fos = (loaded[loaded["Variable"]=="Emissions|CH4|MAGICC Fossil and Industrial"].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze()
        self.fossilCH4_frac = _ch4_fos / self.ch4
        self.aviNOx_frac = _nox_avi / self.nox
        

class Concentrations:
    def __init__(self, loaded, startyear):
        self.year = loaded.columns[(startyear-1700)+7:].astype(int).to_numpy()
        self.co2 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|CO2")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.ch4 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|CH4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.n2o = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|N2O")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cf4 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|PFC|CF4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.c2f6 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|PFC|C2F6")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.c6f14 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|PFC|C6F14")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc23 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC23")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc32 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC32")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc43_10 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC4310mee")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc125 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC125")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc134a = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC134a")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc143a = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC143a")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc227ea = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC227ea")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hfc245fa = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC245fa")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.sf6 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|F-Gases|SF6")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cfc11 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC11")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cfc12 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC12")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cfc113 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC113")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cfc114 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC114")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cfc115 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC115")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.carb_tet = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CCl4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.mcf = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CH3CCl3")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hcfc22 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|HCFC22")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hcfc141b = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|HCFC141b")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.hcfc142b = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|HCFC142b")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.halon1211 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|Halon1211")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.halon1202 = np.zeros(736)
        self.halon1301 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|Halon1301")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.halon2402 = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|Halon2402")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.ch3br = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CH3Br")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.ch3cl = (loaded[(loaded["Variable"]=="Atmospheric Concentrations|Montreal Gases|CH3Cl")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.concentrations = np.vstack((self.year, self.co2, self.ch4, self.n2o, self.cf4, self.c2f6, self.c6f14, self.hfc23, self.hfc32, self.hfc43_10, self.hfc125, self.hfc134a, self.hfc143a, self.hfc227ea, self.hfc245fa, self.sf6, self.cfc11, self.cfc12, self.cfc113, self.cfc114, self.cfc115, self.carb_tet, self.mcf, self.hcfc22, self.hcfc141b, self.hcfc142b, self.halon1211, self.halon1202, self.halon1301, self.halon2402, self.ch3br, self.ch3cl)).T
        self.gas_indices= np.arange(1, 32)


class Forcing:
    def __init__(self, loaded, startyear):
        self.year = loaded.columns[(startyear-1750)+7:].astype(int).to_numpy()
        self.total = (loaded[(loaded["Variable"]=="Effective Radiative Forcing")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.volcanic = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Natural|Volcanic")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.solar = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Natural|Solar")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.co2 = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|CO2")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.ch4 = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|CH4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.n2o = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|N2O")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.other_wmghgs = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|Other WMGHGs")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.aero = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Aerosols|Aerosols-radiation Interactions")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.cloud = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Aerosols|Aerosols-cloud Interactions")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.strato3 = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Stratospheric Ozone")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.tropo3 = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Tropospheric Ozone")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.stwv = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|CH4 Oxidation Stratospheric H2O")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.contrails = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|Contrails and Contrail-induced Cirrus")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.landuse = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|Albedo Change")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
        self.bcsnow = (loaded[(loaded["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|BC on Snow")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    
        
def load_emissions_data(filepath, scenario, startyear=1765):
    # do the loading, return something which has the attributes you want
    loaded = pd.read_csv(filepath)
    loaded = loaded[loaded["Scenario"]==scenario]
    out = Emissions(loaded, startyear)
    return out

def load_concentrations_data(filepath, scenario, startyear=1765):
    loaded = pd.read_csv(filepath)
    loaded = loaded[loaded["Scenario"]==scenario]
    out = Concentrations(loaded, startyear)
    return out

def load_forcing_data(filepath, scenario, startyear=1765):
    loaded = pd.read_csv(filepath)
    loaded = loaded[loaded["Scenario"]==scenario]
    out = Forcing(loaded, startyear)
    return out