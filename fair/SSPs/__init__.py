import os
from ._shared import load_emissions_data, load_concentrations_data, load_forcing_data

emissions_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-emissions-annual-means-5-1-0-ssp-only.csv')
concentrations_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-concentrations-annual-means-5-1-0-ssp-only.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-radiative-forcing-annual-means-5-1-0-ssp-only.csv')


# I know we can refactor this but can't think immediately how
class ssp119:
    Emissions = load_emissions_data(emissions_filename, "ssp119")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp119")
    Forcing = load_forcing_data(forcing_filename, "ssp119")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac
    
    
class ssp126:
    Emissions = load_emissions_data(emissions_filename, "ssp126")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp126")
    Forcing = load_forcing_data(forcing_filename, "ssp126")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac
    
    
class ssp245:
    Emissions = load_emissions_data(emissions_filename, "ssp245")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp245")
    Forcing = load_forcing_data(forcing_filename, "ssp245")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac

    
class ssp370:
    Emissions = load_emissions_data(emissions_filename, "ssp370")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp370")
    Forcing = load_forcing_data(forcing_filename, "ssp370")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac
    

class ssp370_lowNTCF:
    Emissions = load_emissions_data(emissions_filename, "ssp370-lowNTCF-aerchemmip")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp370-lowNTCF-aerchemmip")
    Forcing = load_forcing_data(forcing_filename, "ssp370-lowNTCF-aerchemmip")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac

    
class ssp434:
    Emissions = load_emissions_data(emissions_filename, "ssp434")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp434")
    Forcing = load_forcing_data(forcing_filename, "ssp434")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac
    

class ssp460:
    Emissions = load_emissions_data(emissions_filename, "ssp460")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp460")
    Forcing = load_forcing_data(forcing_filename, "ssp460")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac

    
class ssp534:
    Emissions = load_emissions_data(emissions_filename, "ssp534-over")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp534-over")
    Forcing = load_forcing_data(forcing_filename, "ssp534-over")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac

    
class ssp585:
    Emissions = load_emissions_data(emissions_filename, "ssp585")
    Concentrations = load_concentrations_data(concentrations_filename, "ssp585")
    Forcing = load_forcing_data(forcing_filename, "ssp585")
    fossilCH4_frac = Emissions.fossilCH4_frac
    aviNOx_frac = Emissions.aviNOx_frac