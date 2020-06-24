from __future__ import division

import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ..constants import molwt

def scmdf_to_emissions(scmdf, include_cfcs=True, startyear=1765, endyear=2100):

    """
    Opens a ScmDataFrame and extracts the data. Interpolates linearly
    between non-consecutive years in the SCEN file. Fills in chlorinated gases
    from a specified SSP scenario.

    Note this is a temporary fix for FaIR 1.6.

    Inputs:
        scmdf: ScmDataFrame

    Keywords:
        include_cfcs: bool
            MAGICC files do not come loaded with CFCs (indices 24-39).
            - if True, use the values from RCMIP for SSPs (all scenarios are
                the same).
            - Use False to ignore and create a 24-species emission file.
        startyear: First year of output file.
        endyear: Last year of output file.

    Returns:
        nt x 40 numpy emissions array
    """

    # We expect that aeneris and silicone are going to give us a nicely
    # formatted ScmDataFrame with all 23 species present and correct at
    # timesteps 2015, 2020 and ten-yearly to 2100.
    # We also implicitly assume that data up until 2014 will follow SSP
    # historical.
    # This adapter will not be tested on anything else!

    n_timepoints = scmdf.shape[1]
    n_cols = 40
    nt = endyear - startyear + 1

    data_out = np.ones((nt, n_cols)) * np.nan
    data_out[:,0] = np.arange(startyear, endyear+1)

    # fill in 1765 to 2014 from SSP emissions
    ssp_df = pd.read_csv(os.path.join(os.path.dirname(__file__), '../SSPs/data/rcmip-emissions-annual-means-4-0-0-ssp-only.csv'))

    years = scmdf.columns
    first_scenyear = years[0]
    last_scenyear = years[-1]
    first_row = int(first_scenyear-startyear)
    last_row = int(last_scenyear-startyear)
    
    species = [  # in fair 1.6, order is important
        '|CO2|Energy and Industrial Processes',
        '|CO2|AFOLU',
        '|CH4',
        '|N2O',
        '|Sulfur',
        '|CO',
        '|VOC',
        '|NOx',
        '|BC',
        '|OC',
        '|NH3',
        '|CF4',
        '|C2F6',
        '|C6F14',
        '|HFC23',
        '|HFC32',
        '|HFC43-10',
        '|HFC125',
        '|HFC134a',
        '|HFC143a',
        '|HFC227ea',
        '|HFC245ca',
        '|SF6',
    ]

    emissions_file_species = species.copy()
    emissions_file_species[0] = '|CO2|MAGICC Fossil and Industrial'
    emissions_file_species[1] = '|CO2|MAGICC AFOLU'
    emissions_file_species[16] = '|HFC4310mee'
    emissions_file_species[21] = '|HFC245fa'
    emissions_file_species.extend([
        '|CFC11',
        '|CFC12',
        '|CFC113',
        '|CFC114',
        '|CFC115',
        '|CCl4',
        '|CH3CCl3',
        '|HCFC22',
        '|HCFC141b',
        '|HCFC142b',
        '|Halon1211',
        '|Halon1202',
        '|Halon1301',
        '|Halon2402',
        '|CH3Br',
        '|CH3Cl',
    ])

    # Assume that units coming out of aneris don't change. One day I'll do unit parsing
    unit_convert = np.ones(40)
    unit_convert[1] = molwt.C/molwt.CO2/1000
    unit_convert[2] = molwt.C/molwt.CO2/1000
    unit_convert[4] = molwt.N2/molwt.N2O/1000
    unit_convert[5] = molwt.S/molwt.SO2
    unit_convert[8] = molwt.N/molwt.NO2

    years_future = [2015] + list(range(2020,2501,10))
    for i, specie in enumerate(emissions_file_species):
        data_out[:first_row,i+1] = ssp_df.loc[(ssp_df['Scenario']=='ssp245')&(ssp_df['Variable'].str.endswith(specie)),str(startyear):'2014']*unit_convert[i+1]
        if i<23:
            f = interp1d(years, scmdf[scmdf.index.get_level_values('variable').str.endswith(species[i])])
            data_out[first_row:(last_row+1), i+1] = f(np.arange(first_scenyear, last_scenyear+1))*unit_convert[i+1]
        else:
            f = interp1d(years_future, ssp_df.loc[(ssp_df['Scenario']=='ssp245')&(ssp_df['Variable'].str.endswith(specie)),'2015':'2500'].dropna(axis=1))
            data_out[first_row:(last_row+1), i+1] = f(np.arange(first_scenyear, last_scenyear+1))*unit_convert[i+1]

    return data_out
