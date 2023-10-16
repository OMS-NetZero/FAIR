"""Tools for getting data into and out of FaIR."""

import copy
import os

import numpy as np
import pandas as pd
import pooch
from scipy.interpolate import interp1d

from .exceptions import FromCsvError
from .interface import fill
from .structure.units import (
    compound_convert,
    desired_concentration_units,
    desired_emissions_units,
    mixing_ratio_convert,
    prefix_convert,
    time_convert,
)

HERE = os.path.dirname(os.path.realpath(__file__))
DEFAULT_PROPERTIES_FILE = os.path.join(
    HERE, "defaults", "data", "ar6", "species_configs_properties.csv"
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


def fill_from_csv(self, filename):
    """Fill emissions, concentrations and/or forcing from a CSV file.

    Parameters
    ----------
    filename : str
        path to a csv file.
    """

    # I only require Scenario, Variable and Unit. FaIR may also use Specie.
    # Check present.
    df_input = pd.read_csv(filename)
    df_input.columns = df_input.columns.str.lower()
    df_input.rename(columns={"variable": "specie"}, inplace=True)
    required_columns = ["scenario", "specie", "unit"]
    if not (set(required_columns) < set(df_input.columns)):
        raise FromCsvError(
            f"Input file {filename} must contain 'scenario', 'specie' and 'unit' "
            f"column headers, and at least one timepoint or two timebounds."
        )

    # Now check scenarios
    if not (set(self.scenarios) <= set(df_input["scenario"].unique())):
        raise FromCsvError(
            f"Defined {set(self.scenarios) - set(df_input['scenario'].unique())} scenarios "
            f"were not found in {filename}."
        )

    time = []
    first_time_column = 0
    for column in df_input.columns:
        try:
            time.append(float(column))
        except:
            first_time_column = first_time_column + 1
            # String column afer numeric is an error, probably mis-typed input.
            # Either way, I don't know how to interpolate that, so I won't try.
            if len(time) > 0:
                raise FromCsvError(
                    f"Input file {filename} contains a string column {column} after "
                    f"{time[-1]}. "
                )
    if len(time) < 2:
        raise FromCsvError(
            f"Input file {filename} must contain at least two timepoints or "
            f"timebounds."
        )
    first_time = time[0]
    last_time = time[-1]
    time = np.array(time)

    # Check monotonicity of time axis
    if not np.all(time[1:] >= time[:-1]):
        raise FromCsvError(
            f"Input file {filename} does not have monotonically increasing time "
            f"indices."
        )

    # Use Scipy's interpolate rather than pandas, because the columns are not numeric
    # grab a 1D numpy array with NaNs removed
    for scenario in self.scenarios:
        for specie in self.species:
            if self.properties_df.loc[specie, "input_mode"] == "emissions":
                emis_in = (
                    df_input.loc[
                        (df_input["scenario"] == scenario)
                        & (df_input["specie"] == specie),
                        str(first_time):str(last_time),
                    ]
                    .dropna(axis=1)
                    .values.squeeze()
                )
                print(specie)
                print(emis_in)
                print(df_input)


                # throw error if data missing
                if emis_in.shape[0] == 0:
                    raise ValueError(
                        f"I can't find a value for scenario={scenario}, variable "
                        f"name {specie} in the RCMIP "
                        f"emissions database."
                    )

                # avoid NaNs from outside the interpolation range being mixed into
                # the results
                notnan = np.nonzero(~np.isnan(emis_in))

                # RCMIP are "annual averages"; for emissions this is basically
                # the emissions over the year, for concentrations and forcing
                # it would be midyear values. In every case, we can assume
                # midyear values and interpolate to our time grid.
                rcmip_index = np.arange(1750.5, 2501.5)
                interpolator = interp1d(
                    rcmip_index[notnan],
                    emis_in[notnan],
                    fill_value="extrapolate",
                    bounds_error=False,
                )
                emis = interpolator(self.timepoints)

                # We won't throw an error if the time is out of range for RCMIP,
                # but we will fill with NaN to allow a user to manually specify
                # pre- and post- emissions.
                emis[self.timepoints < 1750] = np.nan
                emis[self.timepoints > 2501] = np.nan

                # Parse and possibly convert unit in input file to what FaIR wants
                unit = df_emis.loc[
                    (df_emis["Scenario"] == scenario)
                    & (df_emis["Variable"].str.endswith("|" + specie_rcmip_name))
                    & (df_emis["Region"] == "World"),
                    "Unit",
                ].values[0]
                emis = emis * (
                    prefix_convert[unit.split()[0]][
                        desired_emissions_units[specie].split()[0]
                    ]
                    * compound_convert[unit.split()[1].split("/")[0]][
                        desired_emissions_units[specie].split()[1].split("/")[0]
                    ]
                    * time_convert[unit.split()[1].split("/")[1]][
                        desired_emissions_units[specie].split()[1].split("/")[1]
                    ]
                )  # * self.timestep

                # fill FaIR xarray
                fill(self.emissions, emis[:, None], specie=specie, scenario=scenario)

            if self.properties_df.loc[specie, "input_mode"] == "concentration":
                # Grab raw concentration from dataframe
                conc_in = (
                    df_conc.loc[
                        (df_conc["Scenario"] == scenario)
                        & (df_conc["Variable"].str.endswith("|" + specie_rcmip_name))
                        & (df_conc["Region"] == "World"),
                        "1700":"2500",
                    ]
                    .interpolate(axis=1)
                    .values.squeeze()
                )

                # throw error if data missing
                if conc_in.shape[0] == 0:
                    raise ValueError(
                        f"I can't find a value for scenario={scenario}, variable "
                        f"name ending with {specie_rcmip_name} in the RCMIP "
                        f"concentration database."
                    )

                # avoid NaNs from outside the interpolation range being mixed into
                # the results
                notnan = np.nonzero(~np.isnan(conc_in))

                # interpolate: this time to timebounds
                rcmip_index = np.arange(1700.5, 2501.5)
                interpolator = interp1d(
                    rcmip_index[notnan],
                    conc_in[notnan],
                    fill_value="extrapolate",
                    bounds_error=False,
                )
                conc = interpolator(self.timebounds)

                # strip out pre- and post-
                conc[self.timebounds < 1700] = np.nan
                conc[self.timebounds > 2501] = np.nan

                # Parse and possibly convert unit in input file to what FaIR wants
                unit = df_conc.loc[
                    (df_conc["Scenario"] == scenario)
                    & (df_conc["Variable"].str.endswith("|" + specie_rcmip_name))
                    & (df_conc["Region"] == "World"),
                    "Unit",
                ].values[0]
                conc = conc * (
                    mixing_ratio_convert[unit][desired_concentration_units[specie]]
                )

                # fill FaIR xarray
                fill(
                    self.concentration,
                    conc[:, None],
                    specie=specie,
                    scenario=scenario,
                )

            if self.properties_df.loc[specie, "input_mode"] == "forcing":
                # Grab raw concentration from dataframe
                forc_in = (
                    df_forc.loc[
                        (df_forc["Scenario"] == scenario)
                        & (df_forc["Variable"].str.endswith("|" + specie_rcmip_name))
                        & (df_forc["Region"] == "World"),
                        "1750":"2500",
                    ]
                    .interpolate(axis=1)
                    .values.squeeze()
                )

                # throw error if data missing
                if forc_in.shape[0] == 0:
                    raise ValueError(
                        f"I can't find a value for scenario={scenario}, variable "
                        f"name ending with {specie_rcmip_name} in the RCMIP "
                        f"radiative forcing database."
                    )

                # avoid NaNs from outside the interpolation range being mixed into
                # the results
                notnan = np.nonzero(~np.isnan(forc_in))

                # interpolate: this time to timebounds
                rcmip_index = np.arange(1750.5, 2501.5)
                interpolator = interp1d(
                    rcmip_index[notnan],
                    forc_in[notnan],
                    fill_value="extrapolate",
                    bounds_error=False,
                )
                forc = interpolator(self.timebounds)

                # strip out pre- and post-
                forc[self.timebounds < 1750] = np.nan
                forc[self.timebounds > 2501] = np.nan

                # Forcing so far is always W m-2, but perhaps this will change.

                # fill FaIR xarray
                fill(self.forcing, forc[:, None], specie=specie, scenario=scenario)


def fill_from_rcmip(self):
    """Fill emissions, concentrations and/or forcing from RCMIP scenarios."""
    # lookup converting FaIR default names to RCMIP names
    species_to_rcmip = {specie: specie.replace("-", "") for specie in self.species}
    species_to_rcmip["CO2 FFI"] = "CO2|MAGICC Fossil and Industrial"
    species_to_rcmip["CO2 AFOLU"] = "CO2|MAGICC AFOLU"
    species_to_rcmip["NOx aviation"] = "NOx|MAGICC Fossil and Industrial|Aircraft"
    species_to_rcmip[
        "Aerosol-radiation interactions"
    ] = "Aerosols-radiation interactions"
    species_to_rcmip["Aerosol-cloud interactions"] = "Aerosols-radiation interactions"
    species_to_rcmip["Contrails"] = "Contrails and Contrail-induced Cirrus"
    species_to_rcmip["Light absorbing particles on snow and ice"] = "BC on Snow"
    species_to_rcmip["Stratospheric water vapour"] = "CH4 Oxidation Stratospheric H2O"
    species_to_rcmip["Land use"] = "Albedo Change"

    species_to_rcmip_copy = copy.deepcopy(species_to_rcmip)

    for specie in species_to_rcmip_copy:
        if specie not in self.species:
            del species_to_rcmip[specie]

    rcmip_emissions_file = pooch.retrieve(
        url="doi:10.5281/zenodo.4589756/rcmip-emissions-annual-means-v5-1-0.csv",
        known_hash="md5:4044106f55ca65b094670e7577eaf9b3",
    )

    rcmip_concentration_file = pooch.retrieve(
        url=(
            "doi:10.5281/zenodo.4589756/" "rcmip-concentrations-annual-means-v5-1-0.csv"
        ),
        known_hash="md5:0d82c3c3cdd4dd632b2bb9449a5c315f",
    )

    rcmip_forcing_file = pooch.retrieve(
        url=(
            "doi:10.5281/zenodo.4589756/"
            "rcmip-radiative-forcing-annual-means-v5-1-0.csv"
        ),
        known_hash="md5:87ef6cd4e12ae0b331f516ea7f82ccba",
    )

    df_emis = pd.read_csv(rcmip_emissions_file)
    df_conc = pd.read_csv(rcmip_concentration_file)
    df_forc = pd.read_csv(rcmip_forcing_file)

    for scenario in self.scenarios:
        for specie, specie_rcmip_name in species_to_rcmip.items():
            if self.properties_df.loc[specie, "input_mode"] == "emissions":
                # Grab raw emissions from dataframe
                emis_in = (
                    df_emis.loc[
                        (df_emis["Scenario"] == scenario)
                        & (df_emis["Variable"].str.endswith("|" + specie_rcmip_name))
                        & (df_emis["Region"] == "World"),
                        "1750":"2500",
                    ]
                    .interpolate(axis=1)
                    .values.squeeze()
                )

                # throw error if data missing
                if emis_in.shape[0] == 0:
                    raise ValueError(
                        f"I can't find a value for scenario={scenario}, variable "
                        f"name ending with {specie_rcmip_name} in the RCMIP "
                        f"emissions database."
                    )

                # avoid NaNs from outside the interpolation range being mixed into
                # the results
                notnan = np.nonzero(~np.isnan(emis_in))

                # RCMIP are "annual averages"; for emissions this is basically
                # the emissions over the year, for concentrations and forcing
                # it would be midyear values. In every case, we can assume
                # midyear values and interpolate to our time grid.
                rcmip_index = np.arange(1750.5, 2501.5)
                interpolator = interp1d(
                    rcmip_index[notnan],
                    emis_in[notnan],
                    fill_value="extrapolate",
                    bounds_error=False,
                )
                emis = interpolator(self.timepoints)

                # We won't throw an error if the time is out of range for RCMIP,
                # but we will fill with NaN to allow a user to manually specify
                # pre- and post- emissions.
                emis[self.timepoints < 1750] = np.nan
                emis[self.timepoints > 2501] = np.nan

                # Parse and possibly convert unit in input file to what FaIR wants
                unit = df_emis.loc[
                    (df_emis["Scenario"] == scenario)
                    & (df_emis["Variable"].str.endswith("|" + specie_rcmip_name))
                    & (df_emis["Region"] == "World"),
                    "Unit",
                ].values[0]
                emis = emis * (
                    prefix_convert[unit.split()[0]][
                        desired_emissions_units[specie].split()[0]
                    ]
                    * compound_convert[unit.split()[1].split("/")[0]][
                        desired_emissions_units[specie].split()[1].split("/")[0]
                    ]
                    * time_convert[unit.split()[1].split("/")[1]][
                        desired_emissions_units[specie].split()[1].split("/")[1]
                    ]
                )  # * self.timestep

                # fill FaIR xarray
                fill(self.emissions, emis[:, None], specie=specie, scenario=scenario)

            if self.properties_df.loc[specie, "input_mode"] == "concentration":
                # Grab raw concentration from dataframe
                conc_in = (
                    df_conc.loc[
                        (df_conc["Scenario"] == scenario)
                        & (df_conc["Variable"].str.endswith("|" + specie_rcmip_name))
                        & (df_conc["Region"] == "World"),
                        "1700":"2500",
                    ]
                    .interpolate(axis=1)
                    .values.squeeze()
                )

                # throw error if data missing
                if conc_in.shape[0] == 0:
                    raise ValueError(
                        f"I can't find a value for scenario={scenario}, variable "
                        f"name ending with {specie_rcmip_name} in the RCMIP "
                        f"concentration database."
                    )

                # avoid NaNs from outside the interpolation range being mixed into
                # the results
                notnan = np.nonzero(~np.isnan(conc_in))

                # interpolate: this time to timebounds
                rcmip_index = np.arange(1700.5, 2501.5)
                interpolator = interp1d(
                    rcmip_index[notnan],
                    conc_in[notnan],
                    fill_value="extrapolate",
                    bounds_error=False,
                )
                conc = interpolator(self.timebounds)

                # strip out pre- and post-
                conc[self.timebounds < 1700] = np.nan
                conc[self.timebounds > 2501] = np.nan

                # Parse and possibly convert unit in input file to what FaIR wants
                unit = df_conc.loc[
                    (df_conc["Scenario"] == scenario)
                    & (df_conc["Variable"].str.endswith("|" + specie_rcmip_name))
                    & (df_conc["Region"] == "World"),
                    "Unit",
                ].values[0]
                conc = conc * (
                    mixing_ratio_convert[unit][desired_concentration_units[specie]]
                )

                # fill FaIR xarray
                fill(
                    self.concentration,
                    conc[:, None],
                    specie=specie,
                    scenario=scenario,
                )

            if self.properties_df.loc[specie, "input_mode"] == "forcing":
                # Grab raw concentration from dataframe
                forc_in = (
                    df_forc.loc[
                        (df_forc["Scenario"] == scenario)
                        & (df_forc["Variable"].str.endswith("|" + specie_rcmip_name))
                        & (df_forc["Region"] == "World"),
                        "1750":"2500",
                    ]
                    .interpolate(axis=1)
                    .values.squeeze()
                )

                # throw error if data missing
                if forc_in.shape[0] == 0:
                    raise ValueError(
                        f"I can't find a value for scenario={scenario}, variable "
                        f"name ending with {specie_rcmip_name} in the RCMIP "
                        f"radiative forcing database."
                    )

                # avoid NaNs from outside the interpolation range being mixed into
                # the results
                notnan = np.nonzero(~np.isnan(forc_in))

                # interpolate: this time to timebounds
                rcmip_index = np.arange(1750.5, 2501.5)
                interpolator = interp1d(
                    rcmip_index[notnan],
                    forc_in[notnan],
                    fill_value="extrapolate",
                    bounds_error=False,
                )
                forc = interpolator(self.timebounds)

                # strip out pre- and post-
                forc[self.timebounds < 1750] = np.nan
                forc[self.timebounds > 2501] = np.nan

                # Forcing so far is always W m-2, but perhaps this will change.

                # fill FaIR xarray
                fill(self.forcing, forc[:, None], specie=specie, scenario=scenario)
