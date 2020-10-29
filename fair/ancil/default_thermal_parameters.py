import numpy as np
import os
import pandas as pd

thermal_params_filename = os.path.join(os.path.dirname(__file__), 'default_thermal_parameters.csv')
thermal_params_df = pd.read_csv(thermal_params_filename, skiprows=1, index_col='Thermal Box').T


def Thermal_Params():
    return thermal_params_df
