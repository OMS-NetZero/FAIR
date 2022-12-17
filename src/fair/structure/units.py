"""Define species units and conversions."""

# These are the emissions units that FaIR expects. One can specify emissions in
# a different unit; in the rcmip files they will get converted automatically.
# You are not required to use these species names, but sticking to these will
# make loading up defaults much less painful.

# publicly importable stuff:
# time_convert
# species_convert
# prefix_convert

import os

import numpy as np
import pandas as pd

from ..earth_params import seconds_per_year

HERE = os.path.dirname(os.path.realpath(__file__))
DEFAULT_PROPERTIES_FILE = os.path.join(
    HERE, "..", "defaults", "data", "ar6", "species_configs_properties.csv"
)

#: Desired emissions units for each specie.
desired_emissions_units = {
    "CO2 FFI": "Gt CO2/yr",
    "CO2 AFOLU": "Gt CO2/yr",
    "CO2": "Gt CO2/yr",
    "CH4": "Mt CH4/yr",
    "N2O": "Mt N2O/yr",
    "Sulfur": "Mt SO2/yr",
    "BC": "Mt BC/yr",
    "OC": "Mt OC/yr",
    "NH3": "Mt NH3/yr",
    "NOx": "Mt NO2/yr",
    "VOC": "Mt VOC/yr",
    "CO": "Mt CO/yr",
    "CFC-11": "kt CFC11/yr",
    "CFC-12": "kt CFC12/yr",
    "CFC-113": "kt CFC113/yr",
    "CFC-114": "kt CFC114/yr",
    "CFC-115": "kt CFC115/yr",
    "HCFC-22": "kt HCFC22/yr",
    "HCFC-141b": "kt HCFC141b/yr",
    "HCFC-142b": "kt HCFC142b/yr",
    "CCl4": "kt CCl4/yr",
    "CHCl3": "kt CHCl3/yr",
    "CH2Cl2": "kt CH2Cl2/yr",
    "CH3Cl": "kt CH3Cl/yr",
    "CH3CCl3": "kt CH3CCl3/yr",
    "CH3Br": "kt CH3Br/yr",
    "Halon-1202": "kt Halon1202/yr",
    "Halon-1211": "kt Halon1211/yr",
    "Halon-1301": "kt Halon1301/yr",
    "Halon-2402": "kt Halon2402/yr",
    "CF4": "kt CF4/yr",
    "C2F6": "kt C2F6/yr",
    "C3F8": "kt C3F8/yr",
    "c-C4F8": "kt cC4F8/yr",
    "C4F10": "kt C4F10/yr",
    "C5F12": "kt C5F12/yr",
    "C6F14": "kt C6F14/yr",
    "C7F16": "kt C7F16/yr",
    "C8F18": "kt C8F18/yr",
    "NF3": "kt NF3/yr",
    "SF6": "kt SF6/yr",
    "SO2F2": "kt SO2F2/yr",
    "HFC-125": "kt HFC125/yr",
    "HFC-134a": "kt HFC134a/yr",
    "HFC-143a": "kt HFC143a/yr",
    "HFC-152a": "kt HFC152a/yr",
    "HFC-227ea": "kt HFC227ea/yr",
    "HFC-23": "kt HFC23/yr",
    "HFC-236fa": "kt HFC236fa/yr",
    "HFC-245fa": "kt HFC245fa/yr",
    "HFC-32": "kt HFC32/yr",
    "HFC-365mfc": "kt HFC365mfc/yr",
    "HFC-4310mee": "kt HFC4310mee/yr",
    "NOx aviation": "Mt NO2/yr",
}

#: Desired concentration units for each default specie.
desired_concentration_units = {
    "CO2": "ppm",
    "CH4": "ppb",
    "N2O": "ppb",
    "CFC-11": "ppt",
    "CFC-12": "ppt",
    "CFC-113": "ppt",
    "CFC-114": "ppt",
    "CFC-115": "ppt",
    "HCFC-22": "ppt",
    "HCFC-141b": "ppt",
    "HCFC-142b": "ppt",
    "CCl4": "ppt",
    "CHCl3": "ppt",
    "CH2Cl2": "ppt",
    "CH3Cl": "ppt",
    "CH3CCl3": "ppt",
    "CH3Br": "ppt",
    "Halon-1202": "ppt",
    "Halon-1211": "ppt",
    "Halon-1301": "ppt",
    "Halon-2402": "ppt",
    "CF4": "ppt",
    "C2F6": "ppt",
    "C3F8": "ppt",
    "c-C4F8": "ppt",
    "C4F10": "ppt",
    "C5F12": "ppt",
    "C6F14": "ppt",
    "C7F16": "ppt",
    "C8F18": "ppt",
    "NF3": "ppt",
    "SF6": "ppt",
    "SO2F2": "ppt",
    "HFC-125": "ppt",
    "HFC-134a": "ppt",
    "HFC-143a": "ppt",
    "HFC-152a": "ppt",
    "HFC-227ea": "ppt",
    "HFC-23": "ppt",
    "HFC-236fa": "ppt",
    "HFC-245fa": "ppt",
    "HFC-32": "ppt",
    "HFC-365mfc": "ppt",
    "HFC-4310mee": "ppt",
    "Equivalent effective stratospheric chlorine": "ppt",
}


# Convert one compound to another based on molecular weight
compound_convert = {}
df = pd.read_csv(DEFAULT_PROPERTIES_FILE, index_col=0)
species_molwts = {
    specie.replace("-", ""): value
    for specie, value in dict(df["molecular_weight"]).items()
}
species_molwts["C"] = 12.011
species_molwts["NO2"] = 46.006
species_molwts["NO"] = 30.006
species_molwts["N"] = 14.007
species_molwts["N2"] = 28.014
species_molwts["SO2"] = 64.069
species_molwts["S"] = 32.07
for specie_from, mw_from in species_molwts.items():
    if ~np.isnan(mw_from):
        compound_convert[specie_from.replace("-", "")] = {}
        for specie_to, mw_to in species_molwts.items():
            if ~np.isnan(mw_to):
                compound_convert[specie_from.replace("-", "")][specie_to] = (
                    mw_to / mw_from
                )


# convert between time units
time_convert = {
    "yr": {
        "yr": 1,
        "day": 1 / (seconds_per_year / 3600 / 24),
        "s": 1 / seconds_per_year,
    },
    "day": {"yr": seconds_per_year / 3600 / 24, "day": 1, "s": 3600 * 24},
    "s": {"yr": seconds_per_year, "day": 1 / 3600 / 24, "s": 1},
}


# convert between prefixes
prefix_convert = {
    "Gt": {"Gt": 1, "Tg": 1000, "Mt": 1000, "kt": 1e6, "t": 1e9, "kg": 1e12, "g": 1e15},
    "Tg": {"Gt": 0.001, "Tg": 1, "Mt": 1, "kt": 1000, "t": 1e6, "kg": 1e9, "g": 1e12},
    "Mt": {"Gt": 0.001, "Tg": 1, "Mt": 1, "kt": 1000, "t": 1e6, "kg": 1e9, "g": 1e12},
    "kt": {"Gt": 1e-6, "Tg": 1e-3, "Mt": 1e-3, "kt": 1, "t": 1e3, "kg": 1e6, "g": 1e9},
    "t": {"Gt": 1e-9, "Tg": 1e-6, "Mt": 1e-6, "kt": 1e-3, "t": 1, "kg": 1e3, "g": 1e6},
    "kg": {
        "Gt": 1e-12,
        "Tg": 1e-9,
        "Mt": 1e-9,
        "kt": 1e-6,
        "t": 1e-3,
        "kg": 1,
        "g": 1e3,
    },
    "g": {
        "Gt": 1e-15,
        "Tg": 1e-12,
        "Mt": 1e-12,
        "kt": 1e-9,
        "t": 1e-6,
        "kg": 1e-3,
        "g": 1,
    },
}

mixing_ratio_convert = {
    "ppm": {"ppm": 1, "ppb": 1000, "ppt": 1e6},
    "ppb": {"ppm": 0.001, "ppb": 1, "ppt": 1000},
    "ppt": {"ppm": 1e-6, "ppb": 0.001, "ppt": 1},
}
