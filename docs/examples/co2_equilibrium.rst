CO2 equilibrium runs
====================

Recreate abrupt 4xCO2 runs from all 66 CMIP6 model calibrations and
create Gregory plots.

This demonstrates the flexibility of FaIR in that it can be applied to
mimic specific ESM experiments.

.. code:: ipython3

    import matplotlib.pyplot as pl
    import numpy as np
    import pandas as pd
    
    from fair import FAIR
    from fair.interface import fill, initialise

1. Create FaIR instance
-----------------------

.. code:: ipython3

    f = FAIR()

2. Define time horizon
----------------------

I want 1000 years, even though 4xCO2 is only 150 year experiment.

.. code:: ipython3

    f.define_time(0, 1000, 1)
    f.timebounds

3. Define scenarios
-------------------

This is easy: there’s only one

.. code:: ipython3

    f.define_scenarios(['abrupt-4xCO2'])

4. Define configs
-----------------

.. code:: ipython3

    df = pd.read_csv("../tests/test_data/4xCO2_cummins_ebm3.csv")
    models = df['model'].unique()
    configs = []
    
    for imodel, model in enumerate(models):
        for run in df.loc[df['model']==model, 'run']:
            configs.append(f"{model}_{run}")
    f.define_configs(configs)

5. Define species
-----------------

Note we set the ``input_mode`` to forcing, as we are running with
prescribed forcing from the 4xCO2 Gregory.

.. code:: ipython3

    species = ['CO2']

.. code:: ipython3

    properties = {
        'CO2': {
            'type': 'co2',
            'input_mode': 'forcing',
            'greenhouse_gas': True,
            'aerosol_chemistry_from_emissions': False,
            'aerosol_chemistry_from_concentration': False,
        },
    }

.. code:: ipython3

    f.define_species(species, properties)

6. Modifying run options
------------------------

Not applicable

7. Create input and output data
-------------------------------

.. code:: ipython3

    f.allocate()

8. fill in everything
---------------------

.. code:: ipython3

    initialise(f.temperature, 0)

.. code:: ipython3

    df = pd.read_csv("../tests/test_data/4xCO2_cummins_ebm3.csv")
    models = df['model'].unique()
    
    seed = 0
    
    for config in configs:
        model, run = config.split('_')
        condition = (df['model']==model) & (df['run']==run)
        fill(f.climate_configs['ocean_heat_capacity'], df.loc[condition, 'C1':'C3'].values.squeeze(), config=config)
        fill(f.climate_configs['ocean_heat_transfer'], df.loc[condition, 'kappa1':'kappa3'].values.squeeze(), config=config)
        fill(f.climate_configs['deep_ocean_efficacy'], df.loc[condition, 'epsilon'].values[0], config=config)
        fill(f.climate_configs['gamma_autocorrelation'], df.loc[condition, 'gamma'].values[0], config=config)
        fill(f.climate_configs['sigma_eta'], df.loc[condition, 'sigma_eta'].values[0], config=config)
        fill(f.climate_configs['sigma_xi'], df.loc[condition, 'sigma_xi'].values[0], config=config)
        fill(f.climate_configs['stochastic_run'], True, config=config)
        fill(f.climate_configs['use_seed'], True, config=config)
        fill(f.climate_configs['seed'], seed, config=config)
        
        # We want to fill in a constant 4xCO2 forcing (for each model) across the run.
        fill(f.forcing, df.loc[condition, 'F_4xCO2'].values[0], config=config, specie='CO2')
        
        seed = seed + 10101

.. code:: ipython3

    df

.. code:: ipython3

    f.fill_species_configs()

.. code:: ipython3

    fill(f.species_configs['tropospheric_adjustment'], 0, specie='CO2')

9. run FaIR
-----------

.. code:: ipython3

    f.run()

10. Show results
----------------

Although we can get convincing internal variability for T and N
individually, it appears that the stochastic variability is correlated.

.. code:: ipython3

    fig, ax = pl.subplots()
    ax.plot(f.timebounds, f.temperature.loc[dict(layer=0, scenario='abrupt-4xCO2')]);
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 13)
    ax.set_ylabel('Global mean warming above pre-industrial, °C')
    ax.set_xlabel('Year')
    ax.set_title('CMIP6 abrupt-4xCO$_2$ emulations, FaIR v2.1')
    fig.tight_layout()

.. code:: ipython3

    pl.plot(f.timebounds, f.toa_imbalance.loc[dict(scenario='abrupt-4xCO2')]);

.. code:: ipython3

    pl.plot(f.timebounds, f.forcing_sum.loc[dict(scenario='abrupt-4xCO2')]);

.. code:: ipython3

    pl.plot(f.timebounds[800:], f.toa_imbalance.loc[dict(scenario='abrupt-4xCO2')][800:,...])
    pl.axhline(0, color='k')

.. code:: ipython3

    fig, ax = pl.subplots(11, 6, figsize=(16, 30))
    
    for i, config in enumerate(configs):
        ax[i//6,i%6].scatter(f.temperature.loc[dict(layer=0, scenario='abrupt-4xCO2', config=config)], f.toa_imbalance.loc[dict(scenario='abrupt-4xCO2', config=config)])
        ax[i//6,i%6].set_xlim(0,13)
        ax[i//6,i%6].set_ylim(-1, 10)
        ax[i//6,i%6].axhline(0, color='k')
        ax[i//6,i%6].set_title(config, fontsize=6)

