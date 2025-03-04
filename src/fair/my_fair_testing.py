import matplotlib.pyplot as pl
import pandas as pd
import numpy as np

from fair import FAIR
from fair.io import read_properties
from fair.interface import fill
from fair.io import read_properties
from fair.interface import initialise

# WARNING: This may not actually be fully emulating UKESM. If you succeed in coupling it
# it's not a very good emulation, you may want to check the code in here to improve the 
# UKESM emulation:
# https://github.com/chrisroadmap/ukesm-fair2.1/blob/main/notebooks/15_try-coupling.ipynb

# ====================
# 1. Initialise FaIR
# ====================

f = FAIR(carbon_method='romero-prieto2024')

# ====================
# 2. Define time horizon
# ====================

f.define_time(1850, 2100, 1)

# ====================
# 3. Define scenarios
# ====================

# scenarios = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp434', 'ssp460', 'ssp534-over', 'ssp585']
scenarios = ['ssp119', 'ssp585']
f.define_scenarios(scenarios)

# ====================
# 4. Define configs ??? how to get the ukesm one?
# ====================

# f.define_configs(['high', 'central', 'low'])

df = pd.read_csv("/home/eearp/code/FAIR/tests/test_data/4xCO2_cummins_ebm3.csv")
models = df['model'].unique()
configs = []

for imodel, model in enumerate(models):
    if model == "UKESM1-0-LL":  # My hack to only get UKESM
        for run in df.loc[df['model']==model, 'run']:
            configs.append(f"{model}_{run}")
f.define_configs(configs)
# This just adds 'UKESM1-0-LL_r1i1p1f2'

# ==================================
# 5. Define species and properties
# ==================================

# species, properties = read_properties('/home/eearp/code/FAIR/examples/data/importing-data/species_configs_properties.csv')
species, properties = read_properties()
f.define_species(species, properties)

# ==================================
# 6. Modify run options (optional)
# ==================================

# f.ghg_method='Myhre1998'

# ==================================
# 7. Initialise arrays
# ==================================

f.allocate()

# ==================================
# 8. Fill in data
# ==================================

# fill(f.emissions, 40, scenario='abrupt', specie='CO2 FFI')

# f.fill_from_csv(
#     emissions_file='/home/eearp/code/FAIR/examples/data/basic_run_example/emissions.csv',
#     concentration_file='/home/eearp/code/FAIR/examples/data/basic_run_example/concentration.csv',
#     forcing_file='/home/eearp/code/FAIR/examples/data/basic_run_example/forcing.csv'
# )

# ==================================
# 8a. Fill in data - species configs
# ==================================

# f.fill_species_configs() # happy with defaults

# f.fill_species_configs('/home/eearp/code/FAIR/examples/data/importing-data/species_configs_properties.csv')

f.fill_species_configs('/home/eearp/code/FAIR/src/fair/defaults/data/ukesm/species_configs_properties_1.4.0_ukesm_methane.csv')

# ==================================
# 8b. Fill in data - emissions
# ==================================

f.fill_from_rcmip() # happy with defaults, but changing volcanic

df_volcanic = pd.read_csv('/home/eearp/code/FAIR/tests/test_data/volcanic_ERF_monthly_175001-201912.csv', index_col='year')

# # overwrite volcanic
# volcanic_forcing = np.zeros(351)
# volcanic_forcing[:271] = df_volcanic[1749:].groupby(np.ceil(df_volcanic[1749:].index) // 1).mean().squeeze().values
# fill(f.forcing, volcanic_forcing[:, None, None], specie="Volcanic")  # sometimes need to expand the array

# initialising 
initialise(f.concentration, f.species_configs['baseline_concentration'])
initialise(f.forcing, 0)
initialise(f.temperature, 0)
initialise(f.cumulative_emissions, 0)
initialise(f.airborne_emissions, 0)

# ==================================
# 8c. Fill in data - climate configs
# ==================================

df = pd.read_csv("/home/eearp/code/FAIR/tests/test_data/4xCO2_cummins_ebm3.csv")
models = df['model'].unique()

seed = 1355763

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
    
    seed = seed + 399

# f.override_defaults('/home/eearp/code/FAIR/examples/data/basic_run_example/configs_ensemble.csv')

# ==================================
# 9 Run FaIR
# ==================================

f.run()

# ==================================
# 10 Pretty plots!
# ==================================

# pl.plot(f.timebounds, f.temperature.loc[dict(scenario='ramp', layer=0)], label=f.configs)
# pl.title('Ramp scenario: temperature')
# pl.xlabel('year')
# pl.ylabel('Temperature anomaly (K)')
# pl.legend()
# pl.show()

pl.plot(f.timebounds, f.concentration.loc[dict(scenario='ssp119', specie='CO2')], label=f.configs);
pl.title('ssp119: CO2 concentration')
pl.xlabel('year')
pl.ylabel('ppm')
pl.show()