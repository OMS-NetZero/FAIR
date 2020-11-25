import pyam
import os

data_filename = os.path.join(
    os.path.dirname(__file__), 'data/RCP6.csv')

df = pyam.IamDataFrame(data=data_filename)

class Emissions:
    emissions = df.filter(variable = 'Emissions*')

class Concentrations:
    concentrations = df.filter(variable = 'Atmospheric Concentrations*')

class Forcing:
    forcing = df.filter(variable = 'Effective Radiative Forcing*')

