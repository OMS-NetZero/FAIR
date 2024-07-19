"""Methods for filling species and climate config parameters in."""

import logging

import pandas as pd

from ..interface import fill

logger = logging.getLogger(__name__)


energy_balance_parameters = [
    "gamma_autocorrelation",
    "ocean_heat_capacity",
    "ocean_heat_transfer",
    "deep_ocean_efficacy",
    "sigma_eta",
    "sigma_xi",
    "forcing_4co2",
    "seed",
    "use_seed",
    "stochastic_run",
]


def override_defaults(self, filename):
    """Fill climate and species for each config from a CSV file.

    This method is part of the `FAIR` class. It uses self.configs to look up
    the config to extract from the given files.

    This method would be used to read in output from the `fair-calibrate` package
    directly into `fair` for the calibrated, constrained parameter sets that are
    produced.

    The column headers are 'config' for the first column (under which appear the config
    names), then the name of the `climate_config` or `species_config` in `fair` (for
    example, 'deep_ocean_efficacy'). When a `climate_config` or `species_config` takes
    a second dimension (often `layer` or `specie`), this is passed inside a square
    bracket in the column header (e.g. 'ocean_heat_transfer[0]'). For `layer` indices,
    a zero-based convention is used.

    Parameters
    ----------
    filename : str
        file location of the configs file
    """
    df_configs = pd.read_csv(filename, index_col=0)
    for config in self.configs:
        logger.debug("Checking for missing config label")
        # warn if config is not present in file; it might be an error by the user, but
        # its not fatal; it can still be filled in before calling run()
        if config not in df_configs.index:
            logger.warning(
                f"I can't find a config with label '{config}' in the supplied file "
                f"{filename}."
            )
            continue
        for col in df_configs.columns:
            logger.debug("Checking whether this is an array")
            # should really do some error checking here; hopefully python's default
            # errors will be obvious enough if a user made a mistake
            if len(col.split("[")) > 1:
                param_name = col.split("[")[0]
                param_index = col.split("[")[1][:-1]
            else:
                param_name = col
                param_index = None

            if param_name in energy_balance_parameters:
                logger.debug(f"Found climate_config parameter {param_name}")
                # error checking required?
                if param_index is not None:
                    logger.debug(f"Filling layer {param_index}")
                    fill(
                        self.climate_configs[param_name],
                        df_configs.loc[config, col],
                        layer=int(param_index),
                        config=config,
                    )
                else:
                    fill(
                        self.climate_configs[param_name],
                        df_configs.loc[config, col],
                        config=config,
                    )

            else:
                logger.debug(f"Found species_config parameter {param_name}")
                # error checking required?
                if param_index is not None:
                    if param_index not in self.species:
                        logger.warning(
                            f"{param_index} is not a specie defined in this `fair` "
                            f"instance for column name {col} in {filename}."
                        )
                        continue
                    logger.debug(f"Filling specie {param_index}")
                    fill(
                        self.species_configs[param_name],
                        df_configs.loc[config, col],
                        specie=param_index,
                        config=config,
                    )
                else:
                    fill(
                        self.species_configs[param_name],
                        df_configs.loc[config, col],
                        config=config,
                    )
