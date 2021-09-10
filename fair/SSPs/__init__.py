import os
from ._shared import load_emissions_data, load_concentrations_data, load_forcing_data

emissions_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-emissions-annual-means-5-1-0-ssp-only.csv')
concentrations_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-concentrations-annual-means-5-1-0-ssp-only.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'data/rcmip-radiative-forcing-annual-means-5-1-0-ssp-only.csv')

class FaIRHelper:
    pass
    
    
def get_ssp(
        emissions_filename,
        concentration_filename,
        forcing_filename,
        scenario
    ):
    emissions = load_emissions_data(emissions_filename, scenario)
    concentrations = load_concentrations_data(concentrations_filename, scenario)
    forcing = load_forcing_data(forcing_filename, scenario)

    out = FaIRHelper()
    out.Emissions = emissions
    out.Concentrations = concentrations
    out.Forcing = forcing
    out.aviNOx_frac = emissions.aviNOx_frac
    out.fossilCH4_frac = emissions.fossilCH4_frac
    return out


ssp119 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp119")
ssp126 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp126")
ssp245 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp245")
ssp370 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp370")
ssp370_lowNTCF = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp370-lowNTCF-aerchemmip")
ssp434 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp434")
ssp460 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp460")
ssp534 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp534-over")
ssp585 = get_ssp(emissions_filename, concentrations_filename, forcing_filename, "ssp585")