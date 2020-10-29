import numpy as np
import os
import pandas as pd

gas_params_filename = os.path.join(os.path.dirname(__file__), 'default_gas_parameters.csv')
gas_params_df = pd.read_csv(gas_params_filename, skiprows=1, index_col='Gas').T


def Gas_Params(gas_list = False):
    if gas_list != False:
        res = gas_params_df[gas_list]
    else:
        res = gas_params_df
    return res
