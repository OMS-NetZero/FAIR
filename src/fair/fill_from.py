"""Methods for filling data from scenario files."""

import copy

import numpy as np
import pandas as pd
import pooch
from scipy.interpolate import interp1d

from .interface import fill
from .structure.units import (
    compound_convert,
    desired_concentration_units,
    desired_emissions_units,
    mixing_ratio_convert,
    prefix_convert,
    time_convert,
)


def fill_from_rcmip(
    self,
    emissions_file=None,
    concentration_file=None,
    forcing_file=None,
):
    """Fill emissions, concentrations and/or forcing from CSV scenarios.

    This method is part of the `FAIR` class. It uses self.scenarios to look up
    the scenario to extract from the given files. If None, the emissions/concentration/
    forcing is filled in from the RCMIP database.

    Parameters
    ----------
    emissions_file : str
        filename of emissions to fill.
    concentration_file : str
        filename of concentrations to fill.
    forcing_file : str
        filename of effective radiative forcing to fill.

    Raises
    ------
    ValueError:
        if the scenario isn't found in the emissions database.
    """
    # lookup converting FaIR default names to RCMIP names
    species_to_rcmip = {specie: specie.replace("-", "") for specie in self.species}
    species_to_rcmip["CO2 FFI"] = "CO2|MAGICC Fossil and Industrial"
    species_to_rcmip["CO2 AFOLU"] = "CO2|MAGICC AFOLU"
    species_to_rcmip["NOx aviation"] = "NOx|MAGICC Fossil and Industrial|Aircraft"
    species_to_rcmip["Aerosol-radiation interactions"] = (
        "Aerosols-radiation interactions"
    )
    species_to_rcmip["Aerosol-cloud interactions"] = "Aerosols-radiation interactions"
    species_to_rcmip["Contrails"] = "Contrails and Contrail-induced Cirrus"
    species_to_rcmip["Light absorbing particles on snow and ice"] = "BC on Snow"
    species_to_rcmip["Stratospheric water vapour"] = "CH4 Oxidation Stratospheric H2O"
    species_to_rcmip["Land use"] = "Albedo Change"

    species_to_rcmip_copy = copy.deepcopy(species_to_rcmip)

    for specie in species_to_rcmip_copy:
        if specie not in self.species:
            del species_to_rcmip[specie]

    if emissions_file is None:
        emissions_file = pooch.retrieve(
            url=(
                "https://zenodo.org/records/4589756/files/"
                "rcmip-emissions-annual-means-v5-1-0.csv"
            ),
            known_hash="md5:4044106f55ca65b094670e7577eaf9b3",
        )

    if concentration_file is None:
        concentration_file = pooch.retrieve(
            url=(
                "https://zenodo.org/records/4589756/files/"
                "rcmip-concentrations-annual-means-v5-1-0.csv"
            ),
            known_hash="md5:0d82c3c3cdd4dd632b2bb9449a5c315f",
        )

    if forcing_file is None:
        forcing_file = pooch.retrieve(
            url=(
                "https://zenodo.org/records/4589756/files/"
                "rcmip-radiative-forcing-annual-means-v5-1-0.csv"
            ),
            known_hash="md5:87ef6cd4e12ae0b331f516ea7f82ccba",
        )

    df_emis = pd.read_csv(emissions_file)
    df_conc = pd.read_csv(concentration_file)
    df_forc = pd.read_csv(forcing_file)

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
