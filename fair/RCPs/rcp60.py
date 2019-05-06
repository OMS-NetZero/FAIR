# Convenience module for loading in RCP emissions datasets
#
# Usage:
#
# import rcp6
# rcp6.Emissions.co2

import os

from . import types


emissions_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP6_EMISSIONS.csv')
concentrations_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP6_MIDYEAR_CONCENTRATIONS.csv')
forcing_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP6_MIDYEAR_RADFORCING.csv')

Emissions = types.Emissions(filename=emissions_filename)

Concentrations = types.Concentrations(filename=concentrations_filename, skiprows=39)

Forcing = types.Forcing(filename=forcing_filename)
