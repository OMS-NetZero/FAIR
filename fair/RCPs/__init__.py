import os

import numpy as np

from ._shared import load_emissions_data, load_concentrations_data, load_forcing_data

aviNOx_filename = os.path.join(os.path.dirname(__file__), 'data/aviNOx_fraction.csv')
fossilCH4_filename = os.path.join(os.path.dirname(__file__), 'data/fossilCH4_fraction.csv')

class rcp26:
    emissions_filename = os.path.join(os.path.dirname(__file__), 'data/RCP3PD_EMISSIONS.csv')
    concentrations_filename = os.path.join(os.path.dirname(__file__), 'data/RCP3PD_MIDYEAR_CONCENTRATIONS.csv')
    forcing_filename = os.path.join(os.path.dirname(__file__), 'data/RCP3PD_MIDYEAR_RADFORCING.csv')
    
    Emissions = load_emissions_data(emissions_filename, skiprows=37)
    Concentrations = load_concentrations_data(concentrations_filename, skiprows=39)
    Forcing = load_forcing_data(forcing_filename, skiprows=59)
    aviNOx_frac = np.loadtxt(aviNOx_filename, skiprows=5, usecols=(1,), delimiter=',')
    fossilCH4_frac = np.loadtxt(fossilCH4_filename, skiprows=5, usecols=(1,), delimiter=',')

rcp3pd = rcp26
    
class rcp45:
    emissions_filename = os.path.join(os.path.dirname(__file__), 'data/RCP45_EMISSIONS.csv')
    concentrations_filename = os.path.join(os.path.dirname(__file__), 'data/RCP45_MIDYEAR_CONCENTRATIONS.csv')
    forcing_filename = os.path.join(os.path.dirname(__file__), 'data/RCP45_MIDYEAR_RADFORCING.csv')
    
    Emissions = load_emissions_data(emissions_filename, skiprows=37)
    Concentrations = load_concentrations_data(concentrations_filename, skiprows=38)
    Forcing = load_forcing_data(forcing_filename, skiprows=59)
    aviNOx_frac = np.loadtxt(aviNOx_filename, skiprows=5, usecols=(1,), delimiter=',')
    fossilCH4_frac = np.loadtxt(fossilCH4_filename, skiprows=5, usecols=(1,), delimiter=',')
    
class rcp60:
    emissions_filename = os.path.join(os.path.dirname(__file__), 'data/RCP6_EMISSIONS.csv')
    concentrations_filename = os.path.join(os.path.dirname(__file__), 'data/RCP6_MIDYEAR_CONCENTRATIONS.csv')
    forcing_filename = os.path.join(os.path.dirname(__file__), 'data/RCP6_MIDYEAR_RADFORCING.csv')
    
    Emissions = load_emissions_data(emissions_filename, skiprows=37)
    Concentrations = load_concentrations_data(concentrations_filename, skiprows=39)
    Forcing = load_forcing_data(forcing_filename, skiprows=59)
    aviNOx_frac = np.loadtxt(aviNOx_filename, skiprows=5, usecols=(1,), delimiter=',')
    fossilCH4_frac = np.loadtxt(fossilCH4_filename, skiprows=5, usecols=(1,), delimiter=',')

rcp6 = rcp60
    
class rcp85:
    emissions_filename = os.path.join(os.path.dirname(__file__), 'data/RCP85_EMISSIONS.csv')
    concentrations_filename = os.path.join(os.path.dirname(__file__), 'data/RCP85_MIDYEAR_CONCENTRATIONS.csv')
    forcing_filename = os.path.join(os.path.dirname(__file__), 'data/RCP85_MIDYEAR_RADFORCING.csv')
    
    Emissions = load_emissions_data(emissions_filename, skiprows=37)
    Concentrations = load_concentrations_data(concentrations_filename, skiprows=39)
    Forcing = load_forcing_data(forcing_filename, skiprows=59)
    aviNOx_frac = np.loadtxt(aviNOx_filename, skiprows=5, usecols=(1,), delimiter=',')
    fossilCH4_frac = np.loadtxt(fossilCH4_filename, skiprows=5, usecols=(1,), delimiter=',')