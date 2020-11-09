import os

import pandas as pd

units_filename = os.path.join(os.path.dirname(__file__), "units.csv")
units_df = pd.read_csv(units_filename, skiprows=1, index_col="Variable")
units_dict = units_df.to_dict()["Unit"]


class Units(dict):
    # TODO: update this so it only loads the csv once
    # rather than once on every instantiation
    def __init__(self):
        self.update(units_dict)
    
    def __getitem__(self, key):
        if key in self:
            val = dict.__getitem__(self, key)
        elif '|' in key:
            sub_key = key[0:key.index('|')]
            val = self.__getitem__(sub_key)
        else:
            val = dict.__getitem__(self, key)
        return val
