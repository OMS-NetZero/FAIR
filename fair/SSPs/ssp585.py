# Convenience module for loading in SSP emissions datasets
#
# Usage:
#
# import ssp585
# ssp585.Emissions.co2

import numpy as np
import pandas as pd
import os

from ..constants import molwt

emissions_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-emissions-annual-means-5-1-0-ssp-only.csv')
concentrations_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-concentrations-annual-means-5-1-0-ssp-only.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-radiative-forcing-annual-means-5-1-0-ssp-only.csv')

scenario = "ssp585"
startyear = 1765
e = pd.read_csv(emissions_filename)
c = pd.read_csv(concentrations_filename)
f = pd.read_csv(forcing_filename)

_nox_avi = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|NOx|MAGICC Fossil and Industrial|Aircraft")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy() * molwt.N / molwt.NO2
_ch4_fos = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CH4|MAGICC Fossil and Industrial")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()

_nox = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|NOx")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy() * molwt.N / molwt.NO2
_ch4 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CH4")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()

fossilCH4_frac = _ch4_fos / _ch4
aviNOx_frac = _nox_avi / _nox

class Emissions:
    year = e.columns[(startyear-1750)+7:].astype(int).to_numpy()
    co2 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CO2")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.C/molwt.CO2 / 1000
    co2_fossil = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CO2|MAGICC Fossil and Industrial")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.C/molwt.CO2 / 1000
    co2_land = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CO2|MAGICC AFOLU")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy().squeeze() * molwt.C/molwt.CO2 / 1000
    ch4 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CH4")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    n2o = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|N2O")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy() / 1000
    sox = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Sulfur")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy() * molwt.S/molwt.SO2
    co = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|CO")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    nmvoc = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|VOC")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    nox = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|NOx")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy() * molwt.N / molwt.NO2
    bc = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|BC")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    oc = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|OC")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    nh3 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|NH3")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    cf4 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|PFC|CF4")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    c2f6 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|PFC|C2F6")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    c6f14 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|PFC|C6F14")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc23 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC23")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc32 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC32")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc43_10 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC4310mee")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc125 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC125")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc134a = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC134a")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc143a = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC143a")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc227ea = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC227ea")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hfc245fa = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|HFC|HFC245fa")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    sf6 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|F-Gases|SF6")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    cfc11 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CFC|CFC11")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    cfc12 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CFC|CFC12")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    cfc113 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CFC|CFC113")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    cfc114 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CFC|CFC114")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    cfc115 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CFC|CFC115")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    carb_tet = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CCl4")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    mcf = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CH3CCl3")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hcfc22 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|HCFC22")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hcfc141b = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|HCFC141b")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    hcfc142b = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|HCFC142b")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    halon1211 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|Halon1211")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    halon1202 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|Halon1202")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    halon1301 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|Halon1301")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    halon2402 = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|Halon2402")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    ch3br = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CH3Br")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()
    ch3cl = (e[(e["Scenario"]==scenario) & (e["Variable"]=="Emissions|Montreal Gases|CH3Cl")].loc[:,str(startyear):].interpolate(axis=1)).astype(float).to_numpy()

    emissions = np.vstack((year, co2_fossil, co2_land, ch4, n2o, sox, co, nmvoc, nox, bc, oc, nh3, cf4, c2f6, c6f14, hfc23, hfc32, hfc43_10, hfc125, hfc134a, hfc143a, hfc227ea, hfc245fa, sf6, cfc11, cfc12, cfc113, cfc114, cfc115, carb_tet, mcf, hcfc22, hcfc141b, hcfc142b, halon1211, halon1202,  halon1301, halon2402, ch3br, ch3cl)).T


class Concentrations:
    year = c.columns[(startyear-1700)+7:].astype(int).to_numpy()
    co2 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|CO2")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    ch4 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|CH4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    n2o = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|N2O")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cf4 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|PFC|CF4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    c2f6 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|PFC|C2F6")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    c6f14 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|PFC|C6F14")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc23 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC23")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc32 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC32")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc43_10 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC4310mee")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc125 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC125")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc134a = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC134a")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc143a = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC143a")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc227ea = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC227ea")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hfc245fa = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|HFC|HFC245fa")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    sf6 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|F-Gases|SF6")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cfc11 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC11")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cfc12 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC12")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cfc113 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC113")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cfc114 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC114")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cfc115 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CFC|CFC115")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    carb_tet = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CCl4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    mcf = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CH3CCl3")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hcfc22 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|HCFC22")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hcfc141b = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|HCFC141b")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    hcfc142b = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|HCFC142b")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    halon1211 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|Halon1211")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    halon1202 = np.zeros(736)
    halon1301 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|Halon1301")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    halon2402 = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|Halon2402")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    ch3br = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CH3Br")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    ch3cl = (c[(c["Scenario"]==scenario) & (c["Variable"]=="Atmospheric Concentrations|Montreal Gases|CH3Cl")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    concentrations = np.vstack((year, co2, ch4, n2o, cf4, c2f6, c6f14, hfc23, hfc32, hfc43_10, hfc125, hfc134a, hfc143a, hfc227ea, hfc245fa, sf6, cfc11, cfc12, cfc113, cfc114, cfc115, carb_tet, mcf, hcfc22, hcfc141b, hcfc142b, halon1211, halon1202, halon1301, halon2402, ch3br, ch3cl)).T
    gas_indices= np.arange(1,32)


class Forcing:
    year = f.columns[(startyear-1750)+7:].astype(int).to_numpy()
    total = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    volcanic = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Natural|Volcanic")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    solar = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Natural|Solar")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    co2 = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|CO2")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    ch4 = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|CH4")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    n2o = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|N2O")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    other_wmghgs = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|Other WMGHGs")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    aero = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Aerosols|Aerosols-radiation Interactions")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    cloud = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Aerosols|Aerosols-cloud Interactions")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    strato3 = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Stratospheric Ozone")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    tropo3 = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Tropospheric Ozone")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    stwv = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|CH4 Oxidation Stratospheric H2O")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    contrails = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|Contrails and Contrail-induced Cirrus")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    landuse = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|Albedo Change")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()
    bcsnow = (f[(f["Scenario"]==scenario) & (f["Variable"]=="Effective Radiative Forcing|Anthropogenic|Other|BC on Snow")].loc[:,str(startyear):]).astype(float).to_numpy().squeeze()