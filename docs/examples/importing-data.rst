Importing emissions files
=========================

This notebook example shows how to run FaIR with your own emissions
scenarios. This feature was introduced in fair v2.2.0.

The data is in the ``examples/data/importing-data`` directory of the
fair repository.

The structure of this example closely follows the ``basic-run-example``.

.. code:: ipython3

    import pandas as pd
    from fair import FAIR
    from fair.io import read_properties
    from fair.interface import fill, initialise
    
    import matplotlib.pyplot as pl

.. code:: ipython3

    f = FAIR()

Set up the run
--------------

.. code:: ipython3

    f.define_time(2000, 2050, 1)  # annual timestep, running from 2000 to 2050
    f.define_scenarios(["renewable", "fossil"])  # define two emissions scenarios
    f.define_configs(["one", "two", "three"])  # three climate ensemble members

In this example the ``species_configs_properties`` are read in from an
external file.

``species_configs_properties`` contains the list of species that you
want to run with, how you want to run them, and the default (though
modifiable) parameter values that you give them.

.. code:: ipython3

    species, properties = read_properties("data/importing-data/species_configs_properties.csv")

.. code:: ipython3

    species

.. code:: ipython3

    f.define_species(species, properties)

.. code:: ipython3

    f.allocate()

Read in our driving data
------------------------

First, we’ll inspect it using ``pandas``.

Note that not every ``specie`` defined is included here, since some are
calculated from other species (CO2, aerosol-radiation and aerosol-cloud
forcing).

Remember also that emissions are on ``timepoints`` - so will be
calculated on half years.

The ``fill_from_csv`` function will do our interpolation for us, so it’s
fine to provide 10-year data for an annual problem.

.. code:: ipython3

    pd.read_csv('data/importing-data/demo-emissions.csv')

.. code:: ipython3

    pd.read_csv('data/importing-data/demo-concentration.csv')

.. code:: ipython3

    pd.read_csv('data/importing-data/demo-forcing.csv')

.. code:: ipython3

    f.fill_from_csv(
        emissions_file='data/importing-data/demo-emissions.csv',
        concentration_file='data/importing-data/demo-concentration.csv',
        forcing_file='data/importing-data/demo-forcing.csv'
    )

Now we fill in the climate and species configs
----------------------------------------------

First take the defaults from the same file as the species/properties
definition:

.. code:: ipython3

    f.fill_species_configs('data/importing-data/species_configs_properties.csv')

Then, for each config set (climate ensemble member) we want to
**override** to the default values in the ``species_config_default``
file. Note that no climate configs are given by default because we want
users to think about what they are doing.

.. code:: ipython3

    df_configs = pd.read_csv('data/importing-data/calibrated_constrained_parameters.csv', index_col=0)

.. code:: ipython3

    energy_balance_parameters = [
        'gamma_autocorrelation',
        'ocean_heat_capacity',
        'ocean_heat_transfer',
        'deep_ocean_efficacy',
        'sigma_eta',
        'sigma_xi',
        'forcing_4co2',
        'seed',
        'use_seed',
        'stochastic_run'
    ]

.. code:: ipython3

    f.climate_configs

.. code:: ipython3

    for config in f.configs:
        for col in df_configs.columns:
            if len(col.split("[")) > 1:
                param_name = col.split("[")[0]
                param_index = (col.split("[")[1][:-1])
            else:
                param_name = col
                param_index = None
    
            if param_name in energy_balance_parameters:
                if param_index is not None:
                    fill(f.climate_configs[param_name], df_configs.loc[config, col], layer=int(param_index), config=config)
                else:
                    fill(f.climate_configs[param_name], df_configs.loc[config, col], config=config)
    
            else:
                if param_index is not None:
                    fill(f.species_configs[param_name], df_configs.loc[config, col], specie=param_index, config=config)
                else:
                    fill(f.species_configs[param_name], df_configs.loc[config, col], config=config)

.. code:: ipython3

    fill(f.climate_configs['stochastic_run'], True)

Initial conditions
------------------

What do we assume at the first time bound (2000.0)?

.. code:: ipython3

    initialise(f.concentration, 278, specie='CO2')
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)

Run
---

.. code:: ipython3

    f.run()

Analyse results
---------------

.. code:: ipython3

    f.temperature

.. code:: ipython3

    for config in f.configs:
        for scenario in f.scenarios:
            pl.plot(f.timebounds, f.temperature.sel(layer=0, scenario=scenario, config=config), label=f"{scenario}:{config}");
    pl.legend()

.. code:: ipython3

    for config in f.configs:
        for scenario in f.scenarios:
            pl.plot(f.timebounds, f.concentration.sel(scenario=scenario, config=config, specie="CO2"), label=f"{scenario}:{config}");
    pl.legend()

.. code:: ipython3

    for config in f.configs:
        for scenario in f.scenarios:
            pl.plot(f.timebounds, f.forcing_sum.sel(scenario=scenario, config=config), label=f"{scenario}:{config}");
    pl.legend()

