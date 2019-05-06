# Convenience module for loading in RCP emissions datasets
#
# Usage:
#
# import rcp45
# rcp45.Emissions.co2

import os

from . import types


emissions_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP45_EMISSIONS.csv')
concentrations_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP45_MIDYEAR_CONCENTRATIONS.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP45_MIDYEAR_RADFORCING.csv')

Emissions = types.Emissions(filename=emissions_filename)

Concentrations = types.Concentrations(filename=concentrations_filename)

Forcing = types.Forcing(filename=forcing_filename)
