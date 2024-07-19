"""Methods for filling data from scenario files."""

import copy
import logging

import numpy as np
import pandas as pd
import pooch
from scipy.interpolate import interp1d

from ..exceptions import (
    DuplicateScenarioError,
    MetaAfterValueError,
    MissingColumnError,
    MissingDataError,
    NonMonotonicError,
    UnitParseError,
)
from ..interface import fill
from ..structure.units import (
    compound_convert,
    desired_concentration_units,
    desired_emissions_units,
    mixing_ratio_convert,
    prefix_convert,
    time_convert,
)

logger = logging.getLogger(__name__)


def _check_csv(df, runmode):
    # check our three metadata columns are present
    required_columns = ["scenario", "region", "variable", "unit"]
    for required_column in required_columns:
        if required_column not in df.columns:
            raise MissingColumnError(
                f"{required_column} is not in the {runmode} file. Please ensure you "
                f"have {required_columns} defined."
            )

    # check that dates come after all metadata and are in chrological order
    # there's no critical reason for erroring here; we can always do sorting on dates
    # and delete columns that shouldn't be there, but there's a chance the user has
    # made a mistake if not.
    times = []
    first_time = False
    for col in df.columns:
        try:
            float(col)
            times.append(col)
            first_time = True
        except ValueError:
            if not first_time:
                pass
            else:
                raise MetaAfterValueError(
                    f"There is a string-value column label {col} in the {runmode} file "
                    "that comes after the first time value. Please check your input "
                    "file and ensure time values are uninterrupted."
                )

    # check strictly monotonically increasing
    if np.any(np.diff(np.array(times, dtype=float)) <= 0):
        raise NonMonotonicError(
            f"Time values in the {runmode} file must be strictly monotonically "
            "increasing."
        )

    return times


def _bounds_warning(firstlast, runmode, filetime, problemtime):
    # Don't raise error if time is out of range because we can still fill in emissions
    # by directly modifying attributes, but user might have made a mistake so warn
    earlierlater = {"first": "later", "last": "earlier"}
    return logger.warning(
        f"The {firstlast} time in the {runmode} file ({filetime}) is "
        f"{earlierlater[firstlast]} than the {firstlast} time in the problem "
        f"definition ({problemtime})."
    )


def _parse_unit(unit, specie, is_ghg):
    try:
        prefix = unit.split()[0]
        compound = unit.split()[1].split("/")[0]
        time = unit.split()[1].split("/")[1]
        logger.debug(f"prefix={prefix}, compound={compound}, time={time}")
    except IndexError:
        raise UnitParseError(
            "Units must be given in the format MASS SPECIE/TIME (with a whitespace "
            "between MASS and SPECIE)"
        )

    # prefix and time need to be from pre-defined list
    if prefix not in prefix_convert:
        raise UnitParseError(
            f"Unit mass given ({prefix}) is not in the list of recognised values, "
            f"which are {list(prefix_convert.keys())}."
        )
    if time not in time_convert:
        raise UnitParseError(
            f"Unit time given ({time}) is not in the list of recognised values, which "
            f"are {list(time_convert.keys())}."
        )

    # compound may be novel if it is user defined, but we can't convert it if so. In
    # which case add to our desired unit lists to prevent later errors.
    if compound not in compound_convert:
        logger.warning(
            f"{compound} is not in fair's default list of species masses for "
            f"{specie}, so I can't convert it. For my non-native species, greenhouse "
            "gas emissions are reported in kt/yr, short-lived forcer emissions in "
            "Mt/yr, greenhouse gas concentrations in ppt, and forcings in W/m2."
        )
        if is_ghg:
            desired_emissions_units[specie] = f"kt {compound}/yr"
            desired_concentration_units[specie] = "ppt"
        else:
            desired_emissions_units[specie] = f"Mt {compound}/yr"
        compound_convert[compound] = {compound: 1}

    return prefix, compound, time


def _emissions_unit_convert(emissions, unit, specie, is_ghg):
    # parse the unit
    prefix, compound, time = _parse_unit(unit, specie, is_ghg)

    emissions = emissions * (
        prefix_convert[prefix][desired_emissions_units[specie].split()[0]]
        * compound_convert[compound][
            desired_emissions_units[specie].split()[1].split("/")[0]
        ]
        * time_convert[time][desired_emissions_units[specie].split()[1].split("/")[1]]
    )  # * self.timestep
    return emissions


def _concentration_unit_convert(concentration, unit, specie):
    if unit not in mixing_ratio_convert:
        raise UnitParseError(
            f"Unit mixing ratio given ({unit}) is not in the list of recognised "
            f"values, which are {list(mixing_ratio_convert.keys())}."
        )
    if specie not in desired_concentration_units:
        logger.warning(
            f"{specie} is not in the default list of greenhouse gases known to fair, "
            "so I'm going to convert concentrations to ppt and report back-calculated "
            "emissions in kt/yr."
        )
        desired_concentration_units[specie] = "ppt"
    concentration = concentration * (
        mixing_ratio_convert[unit][desired_concentration_units[specie]]
    )
    return concentration


def fill_from_csv(
    self, emissions_file=None, concentration_file=None, forcing_file=None
):
    """Fill emissions, concentration and/or forcing from a CSV file.

    This method is part of the `FAIR` class. It uses self.scenarios to look up
    the scenario to extract from the given files.

    Unlike `fill_from_rcmip`, it expects to see species given in the names that they
    have been defined in self.scenarios.

    Additionally, unlike `fill_from_rcmip`, emissions are expected on timepoints.

    CSV files should be provided in horizontal form. Expected columns are Scenario,
    Variable and Unit. Multiple scenarios per file can be supplied. Unit conversion is
    done automatically if the units are properly defined.

    Parameters
    ----------
    emissions_file : str
        filename of emissions to fill.
    concentration_file : str
        filename of concentrations to fill.
    forcing_file : str
        filename of effective radiative forcing to fill.
    """
    mode_options = {
        "emissions": {
            "file": emissions_file,
            "time": self.timepoints,
            "var": self.emissions,
        },
        "concentration": {
            "file": concentration_file,
            "time": self.timebounds,
            "var": self.concentration,
        },
        "forcing": {"file": forcing_file, "time": self.timebounds, "var": self.forcing},
    }
    for mode in mode_options:
        if mode_options[mode]["file"] is not None:
            df = pd.read_csv(mode_options[mode]["file"])
            df.columns = df.columns.str.lower()
            times = _check_csv(df, runmode=mode)  # list of strings
            if float(times[0]) > mode_options[mode]["time"][0]:
                _bounds_warning("first", mode, times[0], mode_options[mode]["time"][0])
            if float(times[-1]) < self.timepoints[-1]:
                _bounds_warning("last", mode, times[-1], mode_options[mode]["time"][-1])
            times_array = np.array(times, dtype=float)

            for scenario in self.scenarios:
                for specie in self.species:
                    if self.properties_df.loc[specie, "input_mode"] == mode:
                        # Grab raw emissions from dataframe
                        data_in = df.loc[
                            (df["scenario"] == scenario)
                            & (df["variable"] == specie)
                            & (df["region"].str.lower() == "world"),
                            times[0] : times[-1],
                        ].values

                        # warn if data missing; it might be an error by the user, but
                        # it's not fatal; we can fill in later
                        if data_in.shape[0] == 0:
                            logger.warning(
                                f"I can't find a value for scenario='{scenario}', "
                                f"variable='{specie}', region='World' in "
                                f"{mode_options[mode]['file']} file."
                            )
                            continue
                        # duplicates are ambigious however, and are an error
                        elif data_in.shape[0] > 1:
                            raise DuplicateScenarioError(
                                f"In {mode_options[mode]['file']} there are duplicate "
                                f"rows for variable='{specie}, scenario='{scenario}'."
                            )
                        # now cast to 1D
                        data_in = data_in.squeeze()

                        # interpolate from the supplied file to our desired timepoints
                        interpolator = interp1d(
                            times_array, data_in, bounds_error=False
                        )
                        data = interpolator(mode_options[mode]["time"])

                        # Parse and possibly convert unit in input to what FaIR wants
                        unit = df.loc[
                            (df["scenario"] == scenario)
                            & (df["variable"] == specie)
                            & (df["region"].str.lower() == "world"),
                            "unit",
                        ].values[0]
                        is_ghg = self.properties_df.loc[specie, "greenhouse_gas"]
                        if mode == "emissions":
                            data = _emissions_unit_convert(data, unit, specie, is_ghg)
                        elif mode == "concentration":
                            data = _concentration_unit_convert(data, unit, specie)

                        # fill FaIR xarray
                        fill(
                            getattr(self, mode),
                            data[:, None],
                            specie=specie,
                            scenario=scenario,
                        )


# TO DO: make part of fill_from_csv
def fill_from_rcmip(self):
    """Fill emissions, concentrations and/or forcing from RCMIP scenarios.

    This method is part of the `FAIR` class. It uses self.scenarios to look up
    the scenario to extract from the given files. If None, the emissions/concentration/
    forcing is filled in from the RCMIP database.
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

    emissions_file = pooch.retrieve(
        url=(
            "https://zenodo.org/records/4589756/files/"
            "rcmip-emissions-annual-means-v5-1-0.csv"
        ),
        known_hash="md5:4044106f55ca65b094670e7577eaf9b3",
    )

    concentration_file = pooch.retrieve(
        url=(
            "https://zenodo.org/records/4589756/files/"
            "rcmip-concentrations-annual-means-v5-1-0.csv"
        ),
        known_hash="md5:0d82c3c3cdd4dd632b2bb9449a5c315f",
    )

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
                    raise MissingDataError(
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
                    raise MissingDataError(
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
                    raise MissingDataError(
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
