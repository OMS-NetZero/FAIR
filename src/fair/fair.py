"""Finite-amplitude Impulse Response (FaIR) simple climate model."""

import os
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from tqdm.auto import tqdm

from .constants import SPECIES_AXIS, TIME_AXIS
from .earth_params import (
    earth_radius,
    mass_atmosphere,
    molecular_weight_air,
    seconds_per_year,
)
from .energy_balance_model import (
    calculate_toa_imbalance_postrun,
    multi_ebm,
    step_temperature,
)
from .forcing.aerosol.erfaci import logsum
from .forcing.aerosol.erfari import calculate_erfari_forcing
from .forcing.ghg import etminan2016, leach2021ghg, meinshausen2020, myhre1998
from .forcing.minor import calculate_linear_forcing
from .forcing.ozone import thornhill2021
from .gas_cycle import calculate_alpha
from .gas_cycle.ch4_lifetime import calculate_alpha_ch4
from .gas_cycle.eesc import calculate_eesc
from .gas_cycle.forward import step_concentration
from .gas_cycle.inverse import unstep_concentration
from .interface import fill
from .structure.species import multiple_allowed, species_types, valid_input_modes
from .structure.species_configs import SPECIES_CONFIGS_EXCL_GASBOX

HERE = os.path.dirname(os.path.realpath(__file__))
DEFAULT_SPECIES_CONFIG_FILE = os.path.join(
    HERE, "defaults", "data", "ar6", "species_configs_properties.csv"
)


class FAIR:
    r"""Initialise FaIR.

    Parameters
    ----------
    n_gasboxes : int
        the number of atmospheric greenhouse gas boxes to run the model with
    n_layers : int
        the number of ocean layers in the energy balance or impulse
        response model to run with
    iirf_max : float
        limit for time-integral of greenhouse gas impulse response function.
    br_cl_ods_potential : float
        factor describing the ratio of efficiency that each bromine atom
        has as an ozone depleting substance relative to each chlorine atom.
    ghg_method : str
        method to use for calculating greenhouse gas forcing from CO\ :sub:`2`,
        CH\ :sub:`4` and N\ :sub:`2`\ O. Valid options are {"leach2021",
        "meinshausen2020", "etminan2016", "myhre1998"}
    ch4_method : str
        method to use for calculating methane lifetime change. Valid options are
        {"leach2021", "thornhill2021"}.
    temperature_prescribed : bool
        Run FaIR with temperatures prescribed.

    Raises
    ------
    ValueError
        if ``ghg_method`` or ``ch4_method`` given are not valid options.
    """

    def __init__(
        self,
        n_gasboxes=4,
        n_layers=3,
        iirf_max=100,
        br_cl_ods_potential=45,
        ghg_method="meinshausen2020",
        ch4_method="leach2021",
        temperature_prescribed=False,
    ):
        """Initialise FaIR."""
        self._ghg_method = ghg_method
        self._ch4_method = ch4_method
        self.gasboxes = range(n_gasboxes)
        self.layers = range(n_layers)
        self.iirf_max = iirf_max
        self.br_cl_ods_potential = br_cl_ods_potential
        self._n_gasboxes = n_gasboxes
        self._n_layers = n_layers
        self.temperature_prescribed = temperature_prescribed

    # attach fill_from methods
    from .io.fill_from import fill_from_csv, fill_from_rcmip
    from .io.param_sets import override_defaults

    # must be a less cumbsersome way to code this
    @property
    def ch4_method(self):
        """Return methane lifetime method."""
        return self._ch4_method.lower()

    @ch4_method.setter
    def ch4_method(self, value):
        if value.lower() in ["thornhill2021", "leach2021"]:
            self._ch4_method = value.lower()
        else:
            raise ValueError(
                f"ch4_method should be ``thornhill2021`` or ``leach2021``; you "
                f"provided {value.lower()}."
            )

    @property
    def ghg_method(self):
        """Return greenhouse gas forcing method."""
        return self._ghg_method.lower()

    @ghg_method.setter
    def ghg_method(self, value):
        if value.lower() in [
            "leach2021",
            "meinshausen2020",
            "etminan2016",
            "myhre1998",
        ]:
            self._ghg_method = value.lower()
        else:
            raise ValueError(
                f"``ghg_method`` should be one of 'leach2021', 'meinshausen2020', "
                f"'etminan2016' or 'myhre1998'; you provided '{value.lower()}'."
            )

    def define_time(self, start, end, step):
        """Define timebounds vector to run FaIR.

        Parameters
        ----------
        start : float
            first timebound of the model (year)
        end : float
            last timebound of the model (year)
        step : float
            timestep (year)
        """
        self.timebounds = np.arange(start, end + step / 2, step)
        self.timepoints = 0.5 * (self.timebounds[1:] + self.timebounds[:-1])
        self.timestep = step
        self._n_timebounds = len(self.timebounds)
        self._n_timepoints = len(self.timepoints)

    def define_scenarios(self, scenarios):
        """Define scenarios to analyse in FaIR.

        Parameters
        ----------
        scenarios : list
            scenario names to run
        """
        self.scenarios = scenarios
        self._n_scenarios = len(scenarios)

    def define_configs(self, configs):
        """Define configs to analyse in FaIR.

        Parameters
        ----------
        configs : list
            config names to run
        """
        self.configs = configs
        self._n_configs = len(configs)

    def define_species(self, species, properties):
        """Define species to run in FaIR.

        Parameters
        ----------
        species : list
            names of species to include in FaIR
        properties : dict
            mapping of each specie to particular run properties. This is a
            nested dict, where each dict key contains a dict of 5 keys as follows:

            ``type`` : str
                the type of specie that is being provided. Valid inputs are
                "co2 ffi", "co2 afolu", "co2", "ch4", "n2o", "cfc-11",
                "other halogen", "f-gas", "sulfur", "black carbon",
                "organic carbon", "other slcf", "nox aviation", "eesc", "ozone",
                "ari", "aci", "contrails", "lapsi", "h2o stratospheric", "land use",
                "volcanic", "solar", "unspecified",
            ``input_mode`` : {'emissions', 'concentration', 'forcing', 'calculated'}
                describes how the specie is input into the model.
            ``greenhouse_gas`` : bool
                is the specie a greenhouse gas?
            ``aerosol_chemistry_from_emissions`` : bool
                does the specie's emissions affect aerosols and/or chemistry?
            ``aerosol_chemistry_from_concentration`` : bool
                does the specie's concentration affect aerosols and/or chemistry?

        Raises
        ------
        ValueError
            if a specie in species does not have a matching key in properties.
        ValueError
            if an invalid species type is specified.
        ValueError
            if an invalid input_type (driving mode) is provided for a particular
            type.
        ValueError
            if duplicate species types are provided for types that must be
            unique.
        """
        self.species = species
        self._n_species = len(species)

        # everything we want to run with defined?
        for specie in species:
            if specie not in properties:
                raise ValueError(
                    f"{specie} does not have a corresponding key in `properties`."
                )

            # everything a valid species type?
            if properties[specie]["type"] not in species_types:
                raise ValueError(
                    f"{properties[specie]['type']} is not a valid species type. Valid "
                    f"types are: {[t for t in species_types]}"
                )

            # input_modes valid?
            if (
                properties[specie]["input_mode"]
                not in valid_input_modes[properties[specie]["type"]]
            ):
                raise ValueError(
                    f"{properties[specie]['input_mode']} is not a valid input mode for "
                    f"{properties[specie]['type']}. Valid input modes are: "
                    f"{[m for m in valid_input_modes[properties[specie]['type']]]}"
                )

        # on the way in, we don't mind if properties is over-specified, but
        # by the time we call allocate(), species and properties must align, so
        # we trim self.properties to match species.
        self.properties = properties
        self.properties_df = pd.DataFrame(self.properties).T.reindex(self.species)

        # 4. check that unique species actually are
        for specie_type in self.properties_df["type"].unique():
            n_repeats = sum(self.properties_df["type"] == specie_type)
            if n_repeats > 1 and not multiple_allowed[specie_type]:
                raise ValueError(
                    f"{specie_type} is defined {n_repeats} times in the problem, but "
                    f"must be unique."
                )

    def allocate(self):
        """Create ``xarray`` DataArrays of data input and output."""
        # check dimensions declared
        required_attributes_and_uncalled_method = {
            "_n_timepoints": "define_time()",
            "_n_scenarios": "define_scenarios()",
            "_n_configs": "define_configs()",
            "_n_species": "define_species()",
        }
        for attr, method in required_attributes_and_uncalled_method.items():
            if not hasattr(self, attr):
                raise AttributeError(
                    f"``FAIR`` object has no attribute '{attr}'. Did you forget to "
                    f"call ``{method}``?"
                )

        # driver/output variables
        self.emissions = xr.DataArray(
            np.ones(
                (
                    self._n_timepoints,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            )
            * np.nan,
            coords=(self.timepoints, self.scenarios, self.configs, self.species),
            dims=("timepoints", "scenario", "config", "specie"),
        )
        self.concentration = xr.DataArray(
            np.ones(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            )
            * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs, self.species),
            dims=("timebounds", "scenario", "config", "specie"),
        )
        self.forcing = xr.DataArray(
            np.ones(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            )
            * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs, self.species),
            dims=("timebounds", "scenario", "config", "specie"),
        )
        self.temperature = xr.DataArray(
            np.ones(
                (self._n_timebounds, self._n_scenarios, self._n_configs, self._n_layers)
            )
            * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs, self.layers),
            dims=("timebounds", "scenario", "config", "layer"),
        )

        # output variables
        self.airborne_emissions = xr.DataArray(
            np.zeros(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            ),
            coords=(self.timebounds, self.scenarios, self.configs, self.species),
            dims=("timebounds", "scenario", "config", "specie"),
        )
        self.alpha_lifetime = xr.DataArray(
            np.ones(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            )
            * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs, self.species),
            dims=("timebounds", "scenario", "config", "specie"),
        )
        self.cumulative_emissions = xr.DataArray(
            np.ones(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            )
            * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs, self.species),
            dims=("timebounds", "scenario", "config", "specie"),
        )
        self.airborne_fraction = xr.DataArray(
            np.ones(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_species,
                )
            )
            * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs, self.species),
            dims=("timebounds", "scenario", "config", "specie"),
        )
        # init with NaNs better than zeros, but makes code cleaner later
        self.ocean_heat_content_change = xr.DataArray(
            np.zeros((self._n_timebounds, self._n_scenarios, self._n_configs)),
            coords=(self.timebounds, self.scenarios, self.configs),
            dims=("timebounds", "scenario", "config"),
        )
        self.toa_imbalance = xr.DataArray(
            np.ones((self._n_timebounds, self._n_scenarios, self._n_configs)) * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs),
            dims=("timebounds", "scenario", "config"),
        )
        self.stochastic_forcing = xr.DataArray(
            np.ones((self._n_timebounds, self._n_scenarios, self._n_configs)) * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs),
            dims=("timebounds", "scenario", "config"),
        )
        self.forcing_sum = xr.DataArray(
            np.ones((self._n_timebounds, self._n_scenarios, self._n_configs)) * np.nan,
            coords=(self.timebounds, self.scenarios, self.configs),
            dims=("timebounds", "scenario", "config"),
        )
        # This is about the only exception to the dimension ordering structure;
        # but testing shows that 5D arrays are too memory consuming and we only
        # want the partitions on the first and last timestep to use in restarts.
        self.gas_partitions = xr.DataArray(
            np.zeros(
                (self._n_scenarios, self._n_configs, self._n_species, self._n_gasboxes)
            ),
            coords=(self.scenarios, self.configs, self.species, self.gasboxes),
            dims=("scenario", "config", "specie", "gasbox"),
        )

        # climate configs
        self.climate_configs = xr.Dataset(
            {
                "ocean_heat_transfer": (
                    ["config", "layer"],
                    np.ones((self._n_configs, self._n_layers)) * np.nan,
                ),
                "ocean_heat_capacity": (
                    ["config", "layer"],
                    np.ones((self._n_configs, self._n_layers)) * np.nan,
                ),
                "deep_ocean_efficacy": ("config", np.ones(self._n_configs) * np.nan),
                "stochastic_run": ("config", np.zeros(self._n_configs, dtype=bool)),
                "sigma_eta": ("config", np.ones(self._n_configs) * 0.5),
                "sigma_xi": ("config", np.ones(self._n_configs) * 0.5),
                "gamma_autocorrelation": ("config", np.ones(self._n_configs) * 2),
                "seed": ("config", np.zeros(self._n_configs, dtype=np.uint32)),
                "use_seed": ("config", np.zeros(self._n_configs, dtype=bool)),
                "forcing_4co2": ("config", np.ones(self._n_configs) * 8),
            },
            coords={"config": self.configs, "layer": self.layers},
        )

        # species configs
        self.species_configs = xr.Dataset(
            {
                # general parameters applicable to all species
                # NB: at present forcing_scale has NO EFFECT on species provided
                # as prescribed forcing.
                "tropospheric_adjustment": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "forcing_efficacy": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)),
                ),
                "forcing_temperature_feedback": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "forcing_scale": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)),
                ),
                # greenhouse gas parameters
                "partition_fraction": (
                    ["config", "specie", "gasbox"],
                    np.ones((self._n_configs, self._n_species, self._n_gasboxes))
                    * np.nan,
                ),
                "unperturbed_lifetime": (
                    ["config", "specie", "gasbox"],
                    np.ones((self._n_configs, self._n_species, self._n_gasboxes))
                    * np.nan,
                ),
                "molecular_weight": ("specie", np.ones(self._n_species) * np.nan),
                "baseline_concentration": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "iirf_0": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "iirf_airborne": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "iirf_uptake": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "iirf_temperature": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "baseline_emissions": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "g0": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "g1": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "concentration_per_emission": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "forcing_reference_concentration": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                # general parameters relating emissions, concentration or forcing of one
                # species to forcing of another.
                # these are all linear factors
                "greenhouse_gas_radiative_efficiency": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "contrails_radiative_efficiency": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "erfari_radiative_efficiency": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "h2o_stratospheric_factor": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "lapsi_radiative_efficiency": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "land_use_cumulative_emissions_to_forcing": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                "ozone_radiative_efficiency": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                # specific parameters for aerosol-cloud interactions
                "aci_scale": (
                    ["config"],
                    np.ones((self._n_configs)) * np.nan,
                ),
                "aci_shape": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                # specific parameters for ozone-depleting GHGs
                "cl_atoms": ("specie", np.zeros(self._n_species)),
                "br_atoms": ("specie", np.zeros(self._n_species)),
                "fractional_release": (
                    ["config", "specie"],
                    np.zeros((self._n_configs, self._n_species)),
                ),
                # specific parameters for methane lifetime
                "ch4_lifetime_chemical_sensitivity": (
                    ["config", "specie"],
                    np.ones((self._n_configs, self._n_species)) * np.nan,
                ),
                "lifetime_temperature_sensitivity": (
                    ["config"],
                    np.ones((self._n_configs)) * np.nan,
                ),
            },
            coords={
                "config": self.configs,
                "specie": self.species,
                "gasbox": self.gasboxes,
            },
        )

    def fill_species_configs(self, filename=DEFAULT_SPECIES_CONFIG_FILE):
        """Fill the species_configs with values from a CSV file.

        Parameters
        ----------
        filename : str
            Path of the CSV file to read the species configs from. If omitted, the
            default configs file will be read.
        """
        df = pd.read_csv(filename, index_col=0)
        for specie in self.species:
            for config in SPECIES_CONFIGS_EXCL_GASBOX:
                fill(
                    self.species_configs[config],
                    df.loc[specie, config],
                    specie=specie,
                )
            for gasbox in range(self._n_gasboxes):
                fill(
                    self.species_configs["partition_fraction"],
                    df.loc[specie, f"partition_fraction{gasbox}"],
                    specie=specie,
                    gasbox=gasbox,
                )
                fill(
                    self.species_configs["unperturbed_lifetime"],
                    df.loc[specie, f"unperturbed_lifetime{gasbox}"],
                    specie=specie,
                    gasbox=gasbox,
                )
        if len(df.loc[df["type"] == "aci"]) > 0:
            fill(
                self.species_configs["aci_scale"],
                df.loc[df["type"] == "aci"].aci_scale,
            )
        if len(df.loc[df["type"] == "ch4"]) > 0:
            fill(
                self.species_configs["lifetime_temperature_sensitivity"],
                df.loc[df["type"] == "ch4"].lifetime_temperature_sensitivity,
            )
        self.calculate_concentration_per_emission()

    # greenhouse gas convenience functions
    def calculate_iirf0(self, iirf_horizon=100):
        r"""Convert greenhouse gas lifetime to time-integrated airborne fraction.

        iirf_0 is the 100-year time-integrated airborne fraction to a pulse
        emission. We know that the gas's atmospheric airborne fraction :math:`a(t)` for
        a gas with lifetime :math:`\tau` after time :math:`t` is therefore

        .. math::
            a(t) = \exp(-t/tau)

        and integrating this for :math:`H` years after a pulse emissions gives us:

        .. math::
            r_0(t) = \int_0^{H} \exp(-t/\tau) dt = \tau (1 - \exp (-H / \tau)).

        :math:`H` = 100 years is the default time horizon in FaIR but this can be set to
        any value.

        Parameters
        ----------
        iirf_horizon : float
            time horizon for time-integrated airborne fraction (yr).
        """
        gasbox_axis = self.species_configs["partition_fraction"].get_axis_num("gasbox")
        self.species_configs["iirf_0"] = np.sum(
            self.species_configs["unperturbed_lifetime"]
            * (1 - np.exp(-iirf_horizon / self.species_configs["unperturbed_lifetime"]))
            * self.species_configs["partition_fraction"],
            gasbox_axis,
        )

    def calculate_g(self, iirf_horizon=100):
        """Calculate lifetime scaling parameters.

        Parameters
        ----------
        iirf_horizon : float
            time horizon for time-integrated airborne fraction (yr).
        """
        gasbox_axis = self.species_configs["partition_fraction"].get_axis_num("gasbox")
        self.species_configs["g1"] = np.sum(
            self.species_configs["partition_fraction"]
            * self.species_configs["unperturbed_lifetime"]
            * (
                1
                - (1 + iirf_horizon / self.species_configs["unperturbed_lifetime"])
                * np.exp(-iirf_horizon / self.species_configs["unperturbed_lifetime"])
            ),
            axis=gasbox_axis,
        )
        self.species_configs["g0"] = np.exp(
            -1
            * np.sum(
                (self.species_configs["partition_fraction"])
                * self.species_configs["unperturbed_lifetime"]
                * (
                    1
                    - np.exp(
                        -iirf_horizon / self.species_configs["unperturbed_lifetime"]
                    )
                ),
                axis=gasbox_axis,
            )
            / self.species_configs["g1"]
        )

    def calculate_concentration_per_emission(
        self, mass_atmosphere=mass_atmosphere, molecular_weight_air=molecular_weight_air
    ):
        """Calculate change in atmospheric concentration for each unit emission."""
        self.species_configs["concentration_per_emission"] = 1 / (
            mass_atmosphere
            / 1e18
            * self.species_configs["molecular_weight"]
            / molecular_weight_air
        )

    # climate response
    def _make_ebms(self):
        # First check for NaNs
        for var in [
            "ocean_heat_capacity",
            "ocean_heat_transfer",
            "deep_ocean_efficacy",
            "gamma_autocorrelation",
        ]:
            if np.isnan(self.climate_configs[var]).sum() > 0:
                raise ValueError(
                    f"There are NaN values in FAIR.climate_configs['{var}']"
                )
        if self.climate_configs["stochastic_run"].sum() > 0:
            for var in ["sigma_eta", "sigma_xi", "seed"]:
                if np.isnan(self.climate_configs[var]).sum() > 0:
                    raise ValueError(
                        f"There are NaN values in climate_configs['{var}'], which is "
                        f"not allowed for FAIR.climate_configs['stochastic_run']=True"
                    )

        self.ebms = multi_ebm(
            self.configs,
            ocean_heat_capacity=self.climate_configs["ocean_heat_capacity"],
            ocean_heat_transfer=self.climate_configs["ocean_heat_transfer"],
            deep_ocean_efficacy=self.climate_configs["deep_ocean_efficacy"],
            stochastic_run=self.climate_configs["stochastic_run"],
            sigma_eta=self.climate_configs["sigma_eta"],
            sigma_xi=self.climate_configs["sigma_xi"],
            gamma_autocorrelation=self.climate_configs["gamma_autocorrelation"],
            seed=self.climate_configs["seed"],
            use_seed=self.climate_configs["use_seed"],
            forcing_4co2=self.climate_configs["forcing_4co2"],
            timestep=self.timestep,
            timebounds=self.timebounds,
        )

    def _check_properties(self):
        def _raise_if_nan(specie, input_mode):
            raise ValueError(
                f"{specie} contains NaN values in its {input_mode} array, which you "
                f"are trying to drive the simulation with."
            )

        self._routine_flags = {
            "ghg": False,
            "ari": False,
            "aci": False,
            "eesc": False,
            "contrails": False,
            "ozone": False,
            "land use": False,
            "lapsi": False,
            "h2o stratospheric": False,
            "temperature": True,
        }
        # check if emissions, concentration, forcing have been defined and
        # that we have non-nan data in every case
        for specie in self.species:
            if self.properties[specie]["input_mode"] == "emissions":
                n_nan = np.isnan(self.emissions.loc[dict(specie=specie)]).sum()
                if n_nan > 0:
                    _raise_if_nan(specie, "emissions")
            elif self.properties[specie]["input_mode"] == "concentration":
                n_nan = np.isnan(self.concentration.loc[dict(specie=specie)]).sum()
                if n_nan > 0:
                    _raise_if_nan(specie, "concentration")
            elif self.properties[specie]["input_mode"] == "forcing":
                n_nan = np.isnan(self.forcing.loc[dict(specie=specie)]).sum()
                if n_nan > 0:
                    _raise_if_nan(specie, "forcing")

        # same for if we are prescribing temperature; we must have non-nan
        # values in the surface level
        if self.temperature_prescribed:
            n_nan = np.isnan(self.temperature.loc[dict(layer=0)]).sum()
            if n_nan > 0:
                raise ValueError(
                    "You are running with prescribed temperatures, but the "
                    "FAIR.temperature xarray contains NaN values in the surface layer."
                )

        # special dependency cases
        if "co2" in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "calculated"][
                "type"
            ]
        ):
            if "co2 ffi" not in list(
                self.properties_df.loc[self.properties_df["input_mode"] == "emissions"][
                    "type"
                ]
            ) or "co2 afolu" not in list(
                self.properties_df.loc[self.properties_df["input_mode"] == "emissions"][
                    "type"
                ]
            ):
                raise ValueError(
                    "co2 in calculated mode requires co2 ffi and co2 afolu in "
                    "emissions mode."
                )

        if "land use" in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "calculated"][
                "type"
            ]
        ):
            if "co2 afolu" not in list(
                self.properties_df.loc[self.properties_df["input_mode"] == "emissions"][
                    "type"
                ]
            ):
                raise ValueError(
                    "land use in calculated mode requires co2 afolu in emissions "
                    "mode."
                )

        if "eesc" not in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "concentration"][
                "type"
            ]
        ) and (
            "eesc"
            in list(
                self.properties_df.loc[
                    self.properties_df["input_mode"] == "calculated"
                ]["type"]
            )
            and (
                "cfc-11"
                not in list(
                    self.properties_df.loc[
                        self.properties_df["input_mode"] == "emissions"
                    ]["type"]
                )
                and "cfc-11"
                not in list(
                    self.properties_df.loc[
                        self.properties_df["input_mode"] == "concentration"
                    ]["type"]
                )
            )
        ):
            if self.ch4_method == "thornhill2021":
                raise ValueError(
                    "For ch4_method = thornhill2021, EESC needs to be input as "
                    "concentrations, or to be calculated from emissions of "
                    "halogenated species which requires at least cfc-11 to be "
                    "provided in emissions or concentration driven mode."
                )

        co2_to_forcing = False
        ch4_to_forcing = False
        n2o_to_forcing = False

        if (
            "co2"
            in list(
                self.properties_df.loc[
                    self.properties_df["input_mode"] == "calculated"
                ]["type"]
            )
            or "co2"
            in list(
                self.properties_df.loc[self.properties_df["input_mode"] == "emissions"][
                    "type"
                ]
            )
            or "co2"
            in list(
                self.properties_df.loc[
                    self.properties_df["input_mode"] == "concentration"
                ]["type"]
            )
        ):
            co2_to_forcing = True
        if "ch4" in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "emissions"][
                "type"
            ]
        ) or "ch4" in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "concentration"][
                "type"
            ]
        ):
            ch4_to_forcing = True
        if "n2o" in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "emissions"][
                "type"
            ]
        ) or "n2o" in list(
            self.properties_df.loc[self.properties_df["input_mode"] == "concentration"][
                "type"
            ]
        ):
            n2o_to_forcing = True
        if self.ghg_method in ["meinshausen2020", "etminan2016"]:
            if 0 < co2_to_forcing + ch4_to_forcing + n2o_to_forcing < 3:
                raise ValueError(
                    "For ghg_method in meinshausen2020, etminan2016, either all of "
                    "co2, ch4 and n2o must be provided in a form that can be "
                    "converted to concentrations, or none"
                )
        elif self.ghg_method == "myhre1998":
            if 0 < ch4_to_forcing + n2o_to_forcing < 2:
                raise ValueError(
                    "for ghg_method=myhre1998, either both of ch4 and n2o must be "
                    "provided, or neither."
                )

        for flag in [
            "ari",
            "aci",
            "ozone",
            "contrails",
            "lapsi",
            "land use",
            "h2o stratospheric",
            "eesc",
        ]:
            if flag in list(
                self.properties_df.loc[
                    self.properties_df["input_mode"] == "calculated"
                ]["type"]
            ):
                self._routine_flags[flag] = True

        # if at least one GHG is emissions, concentration or calculated from
        # precursor emissions, we want to run the forcing calculation
        if (
            (
                self.properties_df.loc[self.properties_df["greenhouse_gas"]].input_mode
                == "concentration"
            ).sum()
            + (
                self.properties_df.loc[self.properties_df["greenhouse_gas"]].input_mode
                == "emissions"
            ).sum()
            + (
                self.properties_df.loc[self.properties_df["greenhouse_gas"]].input_mode
                == "calculated"
            ).sum()
        ):
            self._routine_flags["ghg"] = True

        if self.temperature_prescribed:
            self._routine_flags["temperature"] = False

    def _make_indices(self):
        # the following are all n_species-length boolean arrays

        # these define which species do what in FaIR
        self._ghg_indices = np.asarray(
            self.properties_df.loc[:, "greenhouse_gas"].values, dtype=bool
        )
        self._co2_ffi_indices = np.asarray(
            self.properties_df["type"] == "co2 ffi", dtype=bool
        )
        self._co2_afolu_indices = np.asarray(
            self.properties_df["type"] == "co2 afolu", dtype=bool
        )
        self._co2_indices = np.asarray(self.properties_df["type"] == "co2", dtype=bool)
        self._ch4_indices = np.asarray(self.properties_df["type"] == "ch4", dtype=bool)
        self._n2o_indices = np.asarray(self.properties_df["type"] == "n2o", dtype=bool)
        self._cfc11_indices = np.asarray(
            self.properties_df["type"] == "cfc-11", dtype=bool
        )
        self._sulfur_indices = np.asarray(
            self.properties_df["type"] == "sulfur", dtype=bool
        )
        self._bc_indices = np.asarray(
            self.properties_df["type"] == "black carbon", dtype=bool
        )
        self._oc_indices = np.asarray(
            self.properties_df["type"] == "organic carbon", dtype=bool
        )
        self._aviation_nox_indices = np.asarray(
            self.properties_df["type"] == "aviation nox", dtype=bool
        )
        self._ari_indices = np.asarray(self.properties_df["type"] == "ari", dtype=bool)
        self._aci_indices = np.asarray(self.properties_df["type"] == "aci", dtype=bool)
        self._ozone_indices = np.asarray(
            self.properties_df["type"] == "ozone", dtype=bool
        )
        self._contrails_indices = np.asarray(
            self.properties_df["type"] == "contrails", dtype=bool
        )
        self._lapsi_indices = np.asarray(
            self.properties_df["type"] == "lapsi", dtype=bool
        )
        self._landuse_indices = np.asarray(
            self.properties_df["type"] == "land use", dtype=bool
        )
        self._h2ostrat_indices = np.asarray(
            self.properties_df["type"] == "h2o stratospheric", dtype=bool
        )
        self._eesc_indices = np.asarray(
            self.properties_df["type"] == "eesc", dtype=bool
        )
        self._minor_ghg_indices = (
            self._ghg_indices
            ^ self._co2_indices
            ^ self._ch4_indices
            ^ self._n2o_indices
        )
        self._halogen_indices = self._cfc11_indices | np.asarray(
            self.properties_df["type"] == "other halogen", dtype=bool
        )
        self._aerosol_chemistry_from_emissions_indices = np.asarray(
            self.properties_df.loc[:, "aerosol_chemistry_from_emissions"].values,
            dtype=bool,
        )
        self._aerosol_chemistry_from_concentration_indices = np.asarray(
            self.properties_df.loc[:, "aerosol_chemistry_from_concentration"].values,
            dtype=bool,
        )

        # and these ones are more specific, tripping certain behaviours or functions
        self._ghg_forward_indices = np.asarray(
            (
                (
                    (self.properties_df.loc[:, "input_mode"] == "emissions")
                    | (self.properties_df.loc[:, "input_mode"] == "calculated")
                )
                & (self.properties_df.loc[:, "greenhouse_gas"])
            ).values,
            dtype=bool,
        )
        self._ghg_inverse_indices = np.asarray(
            (
                (self.properties_df.loc[:, "input_mode"] == "concentration")
                & (self.properties_df.loc[:, "greenhouse_gas"])
            ).values,
            dtype=bool,
        )

    def run(self, progress=True, suppress_warnings=True):
        """Run the FaIR model.

        Parameters
        ----------
        progress : bool
            Display progress bar.
        suppress_warnings : bool
            Hide warnings relating to covariance in energy balance matrix.
        """
        self._check_properties()
        self._make_indices()
        if self._routine_flags["temperature"]:
            with warnings.catch_warnings():
                if suppress_warnings:
                    warnings.filterwarnings(
                        "ignore",
                        category=RuntimeWarning,
                        module="scipy.stats._multivariate",
                    )
                self._make_ebms()

        # part of pre-run: TODO move to a new method
        if (
            self._co2_indices.sum()
            + self._co2_ffi_indices.sum()
            + self._co2_afolu_indices.sum()
            == 3
        ):
            self.emissions[..., self._co2_indices] = (
                self.emissions[..., self._co2_ffi_indices].data
                + self.emissions[..., self._co2_afolu_indices].data
            )
        self.cumulative_emissions[1:, ...] = (
            self.emissions.cumsum(dim="timepoints", skipna=False) * self.timestep
            + self.cumulative_emissions[0, ...]
        ).data

        # create numpy arrays
        alpha_lifetime_array = self.alpha_lifetime.data
        airborne_emissions_array = self.airborne_emissions.data
        baseline_concentration_array = self.species_configs[
            "baseline_concentration"
        ].data
        baseline_emissions_array = self.species_configs["baseline_emissions"].data
        br_atoms_array = self.species_configs["br_atoms"].data
        ch4_lifetime_chemical_sensitivity_array = self.species_configs[
            "ch4_lifetime_chemical_sensitivity"
        ].data
        lifetime_temperature_sensitivity_array = self.species_configs[
            "lifetime_temperature_sensitivity"
        ].data
        cl_atoms_array = self.species_configs["cl_atoms"].data
        concentration_array = self.concentration.data
        concentration_per_emission_array = self.species_configs[
            "concentration_per_emission"
        ].data
        contrails_radiative_efficiency_array = self.species_configs[
            "contrails_radiative_efficiency"
        ].data
        cummins_state_array = (
            np.ones(
                (
                    self._n_timebounds,
                    self._n_scenarios,
                    self._n_configs,
                    self._n_layers + 1,
                )
            )
            * np.nan
        )
        cumulative_emissions_array = self.cumulative_emissions.data
        deep_ocean_efficacy_array = self.climate_configs["deep_ocean_efficacy"].data
        emissions_array = self.emissions.data
        erfari_radiative_efficiency_array = self.species_configs[
            "erfari_radiative_efficiency"
        ].data
        erfaci_scale_array = self.species_configs["aci_scale"].data
        erfaci_shape_array = self.species_configs["aci_shape"].data
        forcing_array = self.forcing.data
        forcing_scale_array = self.species_configs["forcing_scale"].data * (
            1 + self.species_configs["tropospheric_adjustment"].data
        )
        forcing_efficacy_array = self.species_configs["forcing_efficacy"].data
        forcing_efficacy_sum_array = (
            np.ones((self._n_timebounds, self._n_scenarios, self._n_configs)) * np.nan
        )
        forcing_reference_concentration_array = self.species_configs[
            "forcing_reference_concentration"
        ].data
        forcing_sum_array = self.forcing_sum.data
        forcing_temperature_feedback_array = self.species_configs[
            "forcing_temperature_feedback"
        ].data
        fractional_release_array = self.species_configs["fractional_release"].data
        g0_array = self.species_configs["g0"].data
        g1_array = self.species_configs["g1"].data
        gas_partitions_array = self.gas_partitions.data
        greenhouse_gas_radiative_efficiency_array = self.species_configs[
            "greenhouse_gas_radiative_efficiency"
        ].data
        h2o_stratospheric_factor_array = self.species_configs[
            "h2o_stratospheric_factor"
        ].data
        iirf_0_array = self.species_configs["iirf_0"].data
        iirf_airborne_array = self.species_configs["iirf_airborne"].data
        iirf_temperature_array = self.species_configs["iirf_temperature"].data
        iirf_uptake_array = self.species_configs["iirf_uptake"].data
        land_use_cumulative_emissions_to_forcing_array = self.species_configs[
            "land_use_cumulative_emissions_to_forcing"
        ].data
        lapsi_radiative_efficiency_array = self.species_configs[
            "lapsi_radiative_efficiency"
        ].data
        ocean_heat_transfer_array = self.climate_configs["ocean_heat_transfer"].data
        ozone_radiative_efficiency_array = self.species_configs[
            "ozone_radiative_efficiency"
        ].data
        partition_fraction_array = self.species_configs["partition_fraction"].data
        unperturbed_lifetime_array = self.species_configs["unperturbed_lifetime"].data

        if self._routine_flags["temperature"]:
            eb_matrix_d_array = self.ebms["eb_matrix_d"].data
            forcing_vector_d_array = self.ebms["forcing_vector_d"].data
            stochastic_d_array = self.ebms["stochastic_d"].data

        # forcing should be initialised so this should not be nan. We could check, or
        # allow silent fail as some species don't take forcings and would correctly be
        # nan.
        forcing_sum_array[0:1, ...] = np.nansum(
            forcing_array[0:1, ...], axis=SPECIES_AXIS
        )

        # this is the most important state vector
        cummins_state_array[0, ..., 0] = forcing_sum_array[0, ...]
        cummins_state_array[..., 1:] = self.temperature.data

        # non-linear forcing relationships need an offset. To save calculating
        # them every timestep, we'll pre-determine the forcing to use as the
        # baseline values.
        # GHGs forcing under Meinshausen2020
        # This check, and others, need to come earlier.
        if self._routine_flags["ghg"] and self.ghg_method == "meinshausen2020":
            if (
                np.sum(
                    np.isnan(
                        forcing_reference_concentration_array[:, self._ghg_indices]
                    )
                )
                > 0
            ):
                raise ValueError(
                    "There are NaNs in "
                    "FAIR.species_configs['forcing_reference_concentration'] which "
                    "means that I can't calculate greenhouse gas forcing."
                )

            # Allow for a user-specified offset (provided through self). This will allow
            # different baseline and pre-industrial concentrations, for example if we
            # want to include natural emissions in CH4. In this case
            # baseline_concentration is zero, but the offset should be w.r.t initial
            # (usually pre-industrial) concentration.
            if not hasattr(self, "ghg_forcing_offset"):
                self.ghg_forcing_offset = meinshausen2020(
                    baseline_concentration_array[None, None, ...],
                    forcing_reference_concentration_array[None, None, ...],
                    forcing_scale_array[None, None, ...],
                    greenhouse_gas_radiative_efficiency_array[None, None, ...],
                    self._co2_indices,
                    self._ch4_indices,
                    self._n2o_indices,
                    self._minor_ghg_indices,
                )

        # Do we also need to check Leach2021 and ozone forcing?

        # it's all been leading up to this : FaIR MAIN LOOP
        for i_timepoint in tqdm(
            range(self._n_timepoints),
            disable=1 - progress,
            desc=f"Running {self._n_scenarios*self._n_configs} projections in parallel",
            unit="timesteps",
        ):
            if self._routine_flags["ghg"]:
                # 1. alpha scaling
                alpha_lifetime_array[
                    i_timepoint : i_timepoint + 1, ..., self._ghg_indices
                ] = calculate_alpha(  # this timepoint
                    airborne_emissions_array[
                        i_timepoint : i_timepoint + 1, ..., self._ghg_indices
                    ],  # last timebound
                    cumulative_emissions_array[
                        i_timepoint : i_timepoint + 1, ..., self._ghg_indices
                    ],  # last timebound
                    g0_array[None, None, ..., self._ghg_indices],
                    g1_array[None, None, ..., self._ghg_indices],
                    iirf_0_array[None, None, ..., self._ghg_indices],
                    iirf_airborne_array[None, None, ..., self._ghg_indices],
                    iirf_temperature_array[None, None, ..., self._ghg_indices],
                    iirf_uptake_array[None, None, ..., self._ghg_indices],
                    cummins_state_array[i_timepoint : i_timepoint + 1, ..., 1:2],
                    self.iirf_max,
                )

                # 2. multi-species methane lifetime if desired; update GHG concentration
                # for CH4
                # needs previous timebound but this is no different to the generic
                if self.ch4_method == "thornhill2021":
                    alpha_lifetime_array[
                        i_timepoint : i_timepoint + 1, ..., self._ch4_indices
                    ] = calculate_alpha_ch4(
                        emissions_array[i_timepoint : i_timepoint + 1, ...],
                        concentration_array[i_timepoint : i_timepoint + 1, ...],
                        cummins_state_array[i_timepoint : i_timepoint + 1, ..., 1:2],
                        baseline_emissions_array[None, None, ...],
                        baseline_concentration_array[None, None, ...],
                        ch4_lifetime_chemical_sensitivity_array[None, None, ...],
                        lifetime_temperature_sensitivity_array[None, None, :, None],
                        self._aerosol_chemistry_from_emissions_indices,
                        self._aerosol_chemistry_from_concentration_indices,
                    )

                # 3. greenhouse emissions to concentrations; include methane from IIRF
                (
                    concentration_array[
                        i_timepoint + 1 : i_timepoint + 2,
                        ...,
                        self._ghg_forward_indices,
                    ],
                    gas_partitions_array[..., self._ghg_forward_indices, :],
                    airborne_emissions_array[
                        i_timepoint + 1 : i_timepoint + 2,
                        ...,
                        self._ghg_forward_indices,
                    ],
                ) = step_concentration(
                    emissions_array[
                        i_timepoint : i_timepoint + 1,
                        ...,
                        self._ghg_forward_indices,
                        None,
                    ],  # this timepoint
                    gas_partitions_array[
                        ..., self._ghg_forward_indices, :
                    ],  # last timebound
                    airborne_emissions_array[
                        i_timepoint + 1 : i_timepoint + 2,
                        ...,
                        self._ghg_forward_indices,
                        None,
                    ],  # last timebound
                    alpha_lifetime_array[
                        i_timepoint : i_timepoint + 1,
                        ...,
                        self._ghg_forward_indices,
                        None,
                    ],
                    baseline_concentration_array[
                        None, None, ..., self._ghg_forward_indices
                    ],
                    baseline_emissions_array[
                        None, None, ..., self._ghg_forward_indices, None
                    ],
                    concentration_per_emission_array[
                        None, None, ..., self._ghg_forward_indices
                    ],
                    unperturbed_lifetime_array[
                        None, None, ..., self._ghg_forward_indices, :
                    ],
                    #        oxidation_matrix,
                    partition_fraction_array[
                        None, None, ..., self._ghg_forward_indices, :
                    ],
                    self.timestep,
                )

                # 4. greenhouse gas concentrations to emissions
                (
                    emissions_array[
                        i_timepoint : i_timepoint + 1, ..., self._ghg_inverse_indices
                    ],
                    gas_partitions_array[..., self._ghg_inverse_indices, :],
                    airborne_emissions_array[
                        i_timepoint + 1 : i_timepoint + 2,
                        ...,
                        self._ghg_inverse_indices,
                    ],
                ) = unstep_concentration(
                    concentration_array[
                        i_timepoint + 1 : i_timepoint + 2,
                        ...,
                        self._ghg_inverse_indices,
                    ],  # this timepoint
                    gas_partitions_array[
                        None, ..., self._ghg_inverse_indices, :
                    ],  # last timebound
                    airborne_emissions_array[
                        i_timepoint : i_timepoint + 1,
                        ...,
                        self._ghg_inverse_indices,
                        None,
                    ],  # last timebound
                    alpha_lifetime_array[
                        i_timepoint : i_timepoint + 1,
                        ...,
                        self._ghg_inverse_indices,
                        None,
                    ],
                    baseline_concentration_array[
                        None, None, ..., self._ghg_inverse_indices
                    ],
                    baseline_emissions_array[
                        None, None, ..., self._ghg_inverse_indices
                    ],
                    concentration_per_emission_array[
                        None, None, ..., self._ghg_inverse_indices
                    ],
                    unperturbed_lifetime_array[
                        None, None, ..., self._ghg_inverse_indices, :
                    ],
                    #        oxidation_matrix,
                    partition_fraction_array[
                        None, None, ..., self._ghg_inverse_indices, :
                    ],
                    self.timestep,
                )
                cumulative_emissions_array[
                    i_timepoint + 1, ..., self._ghg_inverse_indices
                ] = (
                    cumulative_emissions_array[
                        i_timepoint, ..., self._ghg_inverse_indices
                    ]
                    + emissions_array[i_timepoint, ..., self._ghg_inverse_indices]
                    * self.timestep
                )

                # 5. greenhouse gas concentrations to forcing
                if self.ghg_method == "leach2021":
                    forcing_array[
                        i_timepoint + 1 : i_timepoint + 2, ..., self._ghg_indices
                    ] = leach2021ghg(
                        concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                        baseline_concentration_array[None, None, ...]
                        * np.ones(
                            (1, self._n_scenarios, self._n_configs, self._n_species)
                        ),
                        forcing_scale_array[None, None, ...],
                        greenhouse_gas_radiative_efficiency_array[None, None, ...],
                        self._co2_indices,
                        self._ch4_indices,
                        self._n2o_indices,
                        self._minor_ghg_indices,
                    )[
                        0:1, ..., self._ghg_indices
                    ]
                if self.ghg_method == "meinshausen2020":
                    forcing_array[
                        i_timepoint + 1 : i_timepoint + 2, ..., self._ghg_indices
                    ] = meinshausen2020(
                        concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                        forcing_reference_concentration_array[None, None, ...]
                        * np.ones(
                            (1, self._n_scenarios, self._n_configs, self._n_species)
                        ),
                        forcing_scale_array[None, None, ...],
                        greenhouse_gas_radiative_efficiency_array[None, None, ...],
                        self._co2_indices,
                        self._ch4_indices,
                        self._n2o_indices,
                        self._minor_ghg_indices,
                    )[
                        0:1, ..., self._ghg_indices
                    ]
                    forcing_array[
                        i_timepoint + 1 : i_timepoint + 2, ..., self._ghg_indices
                    ] = (
                        forcing_array[
                            i_timepoint + 1 : i_timepoint + 2, ..., self._ghg_indices
                        ]
                        - self.ghg_forcing_offset[..., self._ghg_indices]
                    )
                elif self.ghg_method == "etminan2016":
                    forcing_array[
                        i_timepoint + 1 : i_timepoint + 2, ..., self._ghg_indices
                    ] = etminan2016(
                        concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                        baseline_concentration_array[None, None, ...]
                        * np.ones(
                            (1, self._n_scenarios, self._n_configs, self._n_species)
                        ),
                        forcing_scale_array[None, None, ...],
                        greenhouse_gas_radiative_efficiency_array[None, None, ...],
                        self._co2_indices,
                        self._ch4_indices,
                        self._n2o_indices,
                        self._minor_ghg_indices,
                    )[
                        0:1, ..., self._ghg_indices
                    ]
                elif self.ghg_method == "myhre1998":
                    forcing_array[
                        i_timepoint + 1 : i_timepoint + 2, ..., self._ghg_indices
                    ] = myhre1998(
                        concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                        baseline_concentration_array[None, None, ...]
                        * np.ones(
                            (1, self._n_scenarios, self._n_configs, self._n_species)
                        ),
                        forcing_scale_array[None, None, ...],
                        greenhouse_gas_radiative_efficiency_array[None, None, ...],
                        self._co2_indices,
                        self._ch4_indices,
                        self._n2o_indices,
                        self._minor_ghg_indices,
                    )[
                        0:1, ..., self._ghg_indices
                    ]

            # 6. aerosol direct forcing
            if self._routine_flags["ari"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._ari_indices
                ] = calculate_erfari_forcing(
                    emissions_array[i_timepoint : i_timepoint + 1, ...],
                    concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                    baseline_emissions_array[None, None, ...],
                    baseline_concentration_array[None, None, ...],
                    forcing_scale_array[None, None, ...],
                    erfari_radiative_efficiency_array[None, None, ...],
                    self._aerosol_chemistry_from_emissions_indices,
                    self._aerosol_chemistry_from_concentration_indices,
                )

            # 7. aerosol indirect forcing
            if self._routine_flags["aci"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._aci_indices
                ] = logsum(
                    emissions_array[i_timepoint : i_timepoint + 1, ...],
                    concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                    baseline_emissions_array[None, None, ...],
                    baseline_concentration_array[None, None, ...],
                    forcing_scale_array[None, None, ..., self._aci_indices],
                    erfaci_scale_array[None, None, :],
                    erfaci_shape_array[None, None, ...],
                    self._aerosol_chemistry_from_emissions_indices,
                    self._aerosol_chemistry_from_concentration_indices,
                )

            # 8. calculate EESC this timestep for ozone forcing (and use it for
            # methane lifetime in the following timestep)
            if self._routine_flags["eesc"]:
                concentration_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._eesc_indices
                ] = calculate_eesc(
                    concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                    fractional_release_array[None, None, ...],
                    cl_atoms_array[None, None, ...],
                    br_atoms_array[None, None, ...],
                    self._cfc11_indices,
                    self._halogen_indices,
                    self.br_cl_ods_potential,
                )

            # 9. ozone emissions & concentrations to forcing
            if self._routine_flags["ozone"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._ozone_indices
                ] = thornhill2021(
                    emissions_array[i_timepoint : i_timepoint + 1, ...],
                    concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                    baseline_emissions_array[None, None, ...],
                    baseline_concentration_array[None, None, ...],
                    forcing_scale_array[None, None, ..., self._ozone_indices],
                    ozone_radiative_efficiency_array[None, None, ...],
                    self._aerosol_chemistry_from_emissions_indices,
                    self._aerosol_chemistry_from_concentration_indices,
                )

            # 10. contrails forcing from NOx emissions
            if self._routine_flags["contrails"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._contrails_indices
                ] = calculate_linear_forcing(
                    emissions_array[i_timepoint : i_timepoint + 1, ...],
                    0,
                    forcing_scale_array[None, None, ..., self._contrails_indices],
                    contrails_radiative_efficiency_array[None, None, ...],
                )

            # 11. LAPSI forcing from BC and OC emissions
            if self._routine_flags["lapsi"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._lapsi_indices
                ] = calculate_linear_forcing(
                    emissions_array[i_timepoint : i_timepoint + 1, ...],
                    baseline_emissions_array[None, None, ...],
                    forcing_scale_array[None, None, ..., self._lapsi_indices],
                    lapsi_radiative_efficiency_array[None, None, ...],
                )

            # 12. concentration to stratospheric water vapour forcing
            # 12. concentration to stratospheric water vapour forcing
            if self._routine_flags["h2o stratospheric"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._h2ostrat_indices
                ] = calculate_linear_forcing(
                    concentration_array[i_timepoint + 1 : i_timepoint + 2, ...],
                    baseline_concentration_array[None, None, ...],
                    forcing_scale_array[None, None, ..., self._h2ostrat_indices],
                    h2o_stratospheric_factor_array[None, None, ...],
                )

            # 13. CO2 cumulative emissions to land use change forcing
            if self._routine_flags["land use"]:
                forcing_array[
                    i_timepoint + 1 : i_timepoint + 2, ..., self._landuse_indices
                ] = calculate_linear_forcing(
                    cumulative_emissions_array[i_timepoint + 1 : i_timepoint + 2, ...],
                    0,
                    forcing_scale_array[None, None, ..., self._landuse_indices],
                    land_use_cumulative_emissions_to_forcing_array[None, None, ...],
                )

            # 14. apply temperature-forcing feedback here.
            forcing_array[i_timepoint + 1 : i_timepoint + 2, ...] = (
                forcing_array[i_timepoint + 1 : i_timepoint + 2, ...]
                + cummins_state_array[i_timepoint : i_timepoint + 1, ..., 1:2]
                * forcing_temperature_feedback_array[None, None, ...]
            )

            # 15. sum forcings
            forcing_sum_array[i_timepoint + 1 : i_timepoint + 2, ...] = np.nansum(
                forcing_array[i_timepoint + 1 : i_timepoint + 2, ...], axis=SPECIES_AXIS
            )
            forcing_efficacy_sum_array[i_timepoint + 1 : i_timepoint + 2, ...] = (
                np.nansum(
                    forcing_array[i_timepoint + 1 : i_timepoint + 2, ...]
                    * forcing_efficacy_array[None, None, ...],
                    axis=SPECIES_AXIS,
                )
            )

            # 16. forcing to temperature
            if self._routine_flags["temperature"]:
                cummins_state_array[i_timepoint + 1 : i_timepoint + 2, ...] = (
                    step_temperature(
                        cummins_state_array[i_timepoint : i_timepoint + 1, ...],
                        eb_matrix_d_array[None, None, ...],
                        forcing_vector_d_array[None, None, ...],
                        stochastic_d_array[
                            i_timepoint + 1 : i_timepoint + 2, None, ...
                        ],
                        forcing_efficacy_sum_array[
                            i_timepoint + 1 : i_timepoint + 2, ..., None
                        ],
                    )
                )

        # 17. TOA imbalance
        # forcing is not efficacy adjusted here, is this correct?
        toa_imbalance_array = calculate_toa_imbalance_postrun(
            cummins_state_array,
            forcing_sum_array,  # [..., None],
            ocean_heat_transfer_array,
            deep_ocean_efficacy_array,
        )

        # 18. Ocean heat content change
        # if a restart value is present, include it and subtract off the TOA
        # imbalance from the restart as this would be double counting
        # in non-restart runs, both OHC and N are zero in first timebound
        ocean_heat_content_change_array = self.ocean_heat_content_change[0:1, ...] + (
            (
                np.cumsum(toa_imbalance_array, axis=TIME_AXIS)
                - toa_imbalance_array[0:1, ...]
            )
            * self.timestep
            * earth_radius**2
            * 4
            * np.pi
            * seconds_per_year
        )

        # 19. calculate airborne fraction - we have NaNs and zeros we know about, and we
        # don't mind
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            airborne_fraction_array = (
                airborne_emissions_array / cumulative_emissions_array
            )

        # 20. (Re)allocate to xarray
        self.temperature.data = cummins_state_array[..., 1:]
        self.concentration.data = concentration_array
        self.emissions.data = emissions_array
        self.forcing.data = forcing_array
        self.forcing_sum.data = forcing_sum_array
        self.cumulative_emissions.data = cumulative_emissions_array
        self.airborne_emissions.data = airborne_emissions_array
        self.airborne_fraction.data = airborne_fraction_array
        self.gas_partitions.data = gas_partitions_array
        self.ocean_heat_content_change.data = ocean_heat_content_change_array
        self.toa_imbalance.data = toa_imbalance_array
        self.stochastic_forcing.data = cummins_state_array[..., 0]

    def to_netcdf(self, filename):
        """Write out FaIR scenario data to a netCDF file.

        Parameters
        ----------
        filename : str
            file path of the file to write.
        """
        ds = xr.Dataset(
            data_vars=dict(
                emissions=(
                    ["timepoint", "scenario", "config", "specie"],
                    self.emissions.data,
                ),
                concentration=(
                    ["timebound", "scenario", "config", "specie"],
                    self.concentration.data,
                ),
                forcing=(
                    ["timebound", "scenario", "config", "specie"],
                    self.forcing.data,
                ),
                forcing_sum=(
                    ["timebound", "scenario", "config"],
                    self.forcing_sum.data,
                ),
                temperature=(
                    ["timebound", "scenario", "config", "layer"],
                    self.temperature.data,
                ),
                airborne_emissions=(
                    ["timebound", "scenario", "config", "specie"],
                    self.airborne_emissions.data,
                ),
                airborne_fraction=(
                    ["timebound", "scenario", "config", "specie"],
                    self.airborne_fraction.data,
                ),
                cumulative_emissions=(
                    ["timebound", "scenario", "config", "specie"],
                    self.cumulative_emissions.data,
                ),
                ocean_heat_content_change=(
                    ["timebound", "scenario", "config"],
                    self.ocean_heat_content_change.data,
                ),
                stochastic_forcing=(
                    ["timebound", "scenario", "config"],
                    self.stochastic_forcing.data,
                ),
                toa_imbalance=(
                    ["timebound", "scenario", "config"],
                    self.toa_imbalance.data,
                ),
            ),
            coords=dict(
                timepoint=self.timepoints,
                timebound=self.timebounds,
                scenario=self.scenarios,
                config=self.configs,
                specie=self.species,
                layer=self.layers,
            ),
        )
        ds.to_netcdf(filename)
