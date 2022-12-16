Run the n-layer energy balance model
====================================

This notebook shows examples of extending the 3-layer energy balance
model to general n.

For the two and three layer cases we’ll take the MLE estimates from
Cummins et al. (2020) for HadGEM2-ES, and we’ll use the GISS forcing.
Where n > 3 the data is fake.

.. code:: ipython3

    import matplotlib.pyplot as pl
    import numpy as np
    import pandas as pd
    import pooch
    
    from fair.energy_balance_model import EnergyBalanceModel

.. code:: ipython3

    df_forcing = pd.read_csv('../tests/test_data/RFMIP_ERF_tier2_GISS-E2-1-G.csv')

.. code:: ipython3

    ebm3 = EnergyBalanceModel(
        ocean_heat_capacity=[3.62, 9.47, 98.66],
        ocean_heat_transfer=[0.54, 2.39, 0.63],
        deep_ocean_efficacy=1.59,
        gamma_autocorrelation=1.73,
        sigma_xi=0.32,
        sigma_eta=0.43,
        forcing_4co2=6.35,
        stochastic_run=True,
        seed=16
    )

.. code:: ipython3

    ebm3.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm3.run()

.. code:: ipython3

    ebm3.temperature

.. code:: ipython3

    time = np.arange(1850.5, 2101)

.. code:: ipython3

    pl.plot(time, ebm3.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm3.temperature[:,1], label='second ocean layer')
    pl.plot(time, ebm3.temperature[:,2], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    ebm3.emergent_parameters()
    ebm3.ecs, ebm3.tcr

.. code:: ipython3

    ebm2 = EnergyBalanceModel(
        ocean_heat_capacity=[7.73, 89.29],
        ocean_heat_transfer=[0.63, 0.52],
        deep_ocean_efficacy=1.52,
        gamma_autocorrelation=1.58,
        sigma_xi=0.64,
        sigma_eta=0.43,
        stochastic_run=True,
        forcing_4co2=6.86,
        seed=16
    )

.. code:: ipython3

    ebm2.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm2.emergent_parameters()
    ebm2.ecs, ebm2.tcr

.. code:: ipython3

    ebm2.run()

.. code:: ipython3

    pl.plot(time, ebm2.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm2.temperature[:,1], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    # this is not based on a tuning to any existing CMIP6 model, but I try to get the TCR close to the 
    # HadGEM2 2- and 3-layer cases.
    ebm4 = EnergyBalanceModel(
        ocean_heat_capacity=[1.3, 9, 20, 80],
        ocean_heat_transfer=[0.54, 3, 3, 0.63],
        deep_ocean_efficacy=1.2,
        gamma_autocorrelation=1.73,
        sigma_xi=0.32,
        sigma_eta=0.43,
        forcing_4co2=6.35,
        stochastic_run=True,
        seed=16
    )

.. code:: ipython3

    ebm4.emergent_parameters()
    ebm4.ecs, ebm4.tcr

.. code:: ipython3

    ebm4.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm4.run()

.. code:: ipython3

    pl.plot(time, ebm4.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm4.temperature[:,1], label='second ocean layer')
    pl.plot(time, ebm4.temperature[:,2], label='third ocean layer')
    pl.plot(time, ebm4.temperature[:,3], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    # let's go totally crazy
    ebm10 = EnergyBalanceModel(
        ocean_heat_capacity=[0.6, 1.3, 2, 5, 7, 10, 45, 70, 80, 130],
        ocean_heat_transfer=[0.54, 4, 5, 5, 5, 5, 5, 5, 5, 0.63],
        deep_ocean_efficacy=1.2,
        gamma_autocorrelation=1.73,
        sigma_xi=0.32,
        sigma_eta=0.43,
        forcing_4co2=6.35,
        stochastic_run=True,
        seed=16
    )

.. code:: ipython3

    ebm10.emergent_parameters()
    ebm10.ecs, ebm10.tcr

.. code:: ipython3

    ebm10.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm10.run()

.. code:: ipython3

    pl.plot(time, ebm10.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm10.temperature[:,1], label='second ocean layer')
    pl.plot(time, ebm10.temperature[:,2], label='third ocean layer')
    pl.plot(time, ebm10.temperature[:,3], label='fourth ocean layer')
    pl.plot(time, ebm10.temperature[:,4], label='fifth ocean layer')
    pl.plot(time, ebm10.temperature[:,5], label='sixth ocean layer')
    pl.plot(time, ebm10.temperature[:,6], label='seventh ocean layer')
    pl.plot(time, ebm10.temperature[:,7], label='eighth ocean layer')
    pl.plot(time, ebm10.temperature[:,8], label='ninth ocean layer')
    pl.plot(time, ebm10.temperature[:,9], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    pl.plot(time, ebm2.temperature[:,0], label='two layer model')
    pl.plot(time, ebm3.temperature[:,0], label='three layer model')
    pl.plot(time, ebm4.temperature[:,0], label='four layer model')
    pl.plot(time, ebm10.temperature[:,0], label='ten layer model')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    pl.plot(time, ebm2.toa_imbalance, label='two layer model')
    pl.plot(time, ebm3.toa_imbalance, label='three layer model')
    pl.plot(time, ebm4.toa_imbalance, label='four layer model')
    pl.plot(time, ebm10.toa_imbalance, label='ten layer model')
    pl.ylabel('W/m2 relative to 1850')
    pl.title('SSP2-4.5 TOA radiation change')
    pl.legend()

Repeat everything with stochastic forcing switched off
------------------------------------------------------

.. code:: ipython3

    ebm3 = EnergyBalanceModel(
        ocean_heat_capacity=[3.62, 9.47, 98.66],
        ocean_heat_transfer=[0.54, 2.39, 0.63],
        deep_ocean_efficacy=1.59,
        gamma_autocorrelation=1.73,
        sigma_xi=0.32,
        sigma_eta=0.43,
        forcing_4co2=6.35,
        stochastic_run=False,
        seed=16
    )

.. code:: ipython3

    ebm3.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm3.run()

.. code:: ipython3

    ebm3.temperature

.. code:: ipython3

    time = np.arange(1850.5, 2101)

.. code:: ipython3

    pl.plot(time, ebm3.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm3.temperature[:,1], label='second ocean layer')
    pl.plot(time, ebm3.temperature[:,2], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    ebm3.emergent_parameters()
    ebm3.ecs, ebm3.tcr

.. code:: ipython3

    ebm2 = EnergyBalanceModel(
        ocean_heat_capacity=[7.73, 89.29],
        ocean_heat_transfer=[0.63, 0.52],
        deep_ocean_efficacy=1.52,
        gamma_autocorrelation=1.58,
        sigma_xi=0.64,
        sigma_eta=0.43,
        stochastic_run=False,
        forcing_4co2=6.86,
        seed=16
    )

.. code:: ipython3

    ebm2.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm2.emergent_parameters()
    ebm2.ecs, ebm2.tcr

.. code:: ipython3

    ebm2.run()

.. code:: ipython3

    pl.plot(time, ebm2.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm2.temperature[:,1], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    # this is not based on a tuning to any existing CMIP6 model, but I try to get the TCR close to the 
    # HadGEM2 2- and 3-layer cases.
    ebm4 = EnergyBalanceModel(
        ocean_heat_capacity=[1.3, 9, 20, 80],
        ocean_heat_transfer=[0.54, 3, 3, 0.63],
        deep_ocean_efficacy=1.2,
        gamma_autocorrelation=1.73,
        sigma_xi=0.32,
        sigma_eta=0.43,
        forcing_4co2=6.35,
        stochastic_run=False,
        seed=16
    )

.. code:: ipython3

    ebm4.emergent_parameters()
    ebm4.ecs, ebm4.tcr

.. code:: ipython3

    ebm4.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm4.run()

.. code:: ipython3

    pl.plot(time, ebm4.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm4.temperature[:,1], label='second ocean layer')
    pl.plot(time, ebm4.temperature[:,2], label='third ocean layer')
    pl.plot(time, ebm4.temperature[:,3], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    # let's go totally crazy
    ebm10 = EnergyBalanceModel(
        ocean_heat_capacity=[0.6, 1.3, 2, 5, 7, 10, 45, 70, 80, 130],
        ocean_heat_transfer=[0.54, 4, 5, 5, 5, 5, 5, 5, 5, 0.63],
        deep_ocean_efficacy=1.2,
        gamma_autocorrelation=1.73,
        sigma_xi=0.32,
        sigma_eta=0.43,
        forcing_4co2=6.35,
        stochastic_run=False,
        seed=16
    )

.. code:: ipython3

    ebm10.emergent_parameters()
    ebm10.ecs, ebm10.tcr

.. code:: ipython3

    ebm10.add_forcing(forcing = df_forcing['GISS-E2-1-G TOT'].values, timestep=1)

.. code:: ipython3

    ebm10.run()

.. code:: ipython3

    pl.plot(time, ebm10.temperature[:,0], label='surface / top ocean layer')
    pl.plot(time, ebm10.temperature[:,1], label='second ocean layer')
    pl.plot(time, ebm10.temperature[:,2], label='third ocean layer')
    pl.plot(time, ebm10.temperature[:,3], label='fourth ocean layer')
    pl.plot(time, ebm10.temperature[:,4], label='fifth ocean layer')
    pl.plot(time, ebm10.temperature[:,5], label='sixth ocean layer')
    pl.plot(time, ebm10.temperature[:,6], label='seventh ocean layer')
    pl.plot(time, ebm10.temperature[:,7], label='eighth ocean layer')
    pl.plot(time, ebm10.temperature[:,8], label='ninth ocean layer')
    pl.plot(time, ebm10.temperature[:,9], label='deep ocean layer')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    pl.plot(time, ebm2.temperature[:,0], label='two layer model')
    pl.plot(time, ebm3.temperature[:,0], label='three layer model')
    pl.plot(time, ebm4.temperature[:,0], label='four layer model')
    pl.plot(time, ebm10.temperature[:,0], label='ten layer model')
    pl.ylabel('K relative to 1850')
    pl.title('SSP2-4.5 temperature change')
    pl.legend()

.. code:: ipython3

    pl.plot(time, ebm2.toa_imbalance, label='two layer model')
    pl.plot(time, ebm3.toa_imbalance, label='three layer model')
    pl.plot(time, ebm4.toa_imbalance, label='four layer model')
    pl.plot(time, ebm10.toa_imbalance, label='ten layer model')
    pl.ylabel('W/m2 relative to 1850')
    pl.title('SSP2-4.5 TOA radiation change')
    pl.legend()

