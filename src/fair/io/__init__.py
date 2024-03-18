"""Tools for getting data into and out of FaIR."""

import os

import pandas as pd

HERE = os.path.dirname(os.path.realpath(__file__))
DEFAULT_PROPERTIES_FILE = os.path.join(
    HERE, "..", "defaults", "data", "ar6", "species_configs_properties.csv"
)

_default_ghg_and_slcfs = [
    "CO2 FFI",
    "CO2 AFOLU",
    "CO2",
    "CH4",
    "N2O",
    "Sulfur",
    "BC",
    "OC",
    "NH3",
    "NOx",
    "VOC",
    "CO",
    "CFC-11",
    "CFC-12",
    "CFC-113",
    "CFC-114",
    "CFC-115",
    "HCFC-22",
    "HCFC-141b",
    "HCFC-142b",
    "CCl4",
    "CHCl3",
    "CH2Cl2",
    "CH3Cl",
    "CH3CCl3",
    "CH3Br",
    "Halon-1211",
    "Halon-1301",
    "Halon-2402",
    "CF4",
    "C2F6",
    "C3F8",
    "c-C4F8",
    "C4F10",
    "C5F12",
    "C6F14",
    "C7F16",
    "C8F18",
    "NF3",
    "SF6",
    "SO2F2",
    "HFC-125",
    "HFC-134a",
    "HFC-143a",
    "HFC-152a",
    "HFC-227ea",
    "HFC-23",
    "HFC-236fa",
    "HFC-245fa",
    "HFC-32",
    "HFC-365mfc",
    "HFC-4310mee",
    "NOx aviation",
]


def read_properties(filename=DEFAULT_PROPERTIES_FILE, species=None):
    """Get a properties file.

    Parameters
    ----------
    filename : str
        path to a csv file. Default is an AR6 WG1-like config for FaIR
        covering all of the species considered in CMIP6.
    species : list of str or None
        the species that are to be included in the FaIR run. All of these
        species should be present in the index (first column) of the csv. If
        None (default), return all of the species in the defaults.

    Returns
    -------
    species : list
        a list of species names that are included in the FaIR run.
    properties : dict
        species properties that control the FaIR run
    """
    df = pd.read_csv(filename, index_col=0)

    if species is None:
        species = list(df.index)

    properties = {}
    for specie in species:
        properties[specie] = {
            "type": df.loc[specie].type,
            "input_mode": df.loc[specie].input_mode,
            "greenhouse_gas": bool(df.loc[specie].greenhouse_gas),
            "aerosol_chemistry_from_emissions": bool(
                df.loc[specie].aerosol_chemistry_from_emissions
            ),
            "aerosol_chemistry_from_concentration": bool(
                df.loc[specie].aerosol_chemistry_from_concentration
            ),
        }
    return species, properties
