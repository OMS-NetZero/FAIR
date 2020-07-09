from __future__ import division

import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ..constants import molwt

try:
    from scmdata import ScmDataFrame
    has_scmdata = True
except ImportError:
    has_scmdata = False


EMISSIONS_SPECIES_UNITS_CONTEXT = (  # in fair 1.6, order is important
    # @chrisroadmap can you check please
    ('|CO2|MAGICC Fossil and Industrial', 'GtC / yr', None),
    ('|CO2|MAGICC AFOLU', 'GtC / yr', None),
    ('|CH4', 'MtCH4 / yr', None),
    ('|N2O', 'MtN / yr', None),
    ('|Sulfur', 'MtS / yr', None),
    ('|CO', 'MtCO / yr', None),
    ('|VOC', 'MtNMVOC / yr', None),
    ('|NOx', 'MtN / yr', "NOx_conversions"),
    ('|BC', 'MtBC / yr', None),
    ('|OC', 'MtOC / yr', None),
    ('|NH3', 'MtN / yr', None),
    ('|CF4', 'ktCF4 / yr', None),
    ('|C2F6', 'ktC2F6 / yr', None),
    ('|C6F14', 'ktC6F14 / yr', None),
    ('|HFC23', 'ktHFC23 / yr', None),
    ('|HFC32', 'ktHFC32 / yr', None),
    ('|HFC4310mee', 'ktHFC4310mee / yr', None),
    ('|HFC125', 'ktHFC125 / yr', None),
    ('|HFC134a', 'ktHFC134a / yr', None),
    ('|HFC143a', 'ktHFC143a / yr', None),
    ('|HFC227ea', 'ktHFC227ea / yr', None),
    ('|HFC245fa', 'ktHFC245fa / yr', None),
    ('|SF6', 'ktSF6 / yr', None),
    ('|CFC11', 'ktCFC11 / yr', None),
    ('|CFC12', 'ktCFC12 / yr', None),
    ('|CFC113', 'ktCFC113 / yr', None),
    ('|CFC114', 'ktCFC114 / yr', None),
    ('|CFC115', 'ktCFC115 / yr', None),
    ('|CCl4', 'ktCCl4 / yr', None),
    ('|CH3CCl3', 'ktCH3CCl3 / yr', None),
    ('|HCFC22', 'ktHCFC22 / yr', None),
    ('|HCFC141b', 'ktHCFC141b / yr', None),
    ('|HCFC142b', 'ktHCFC142b / yr', None),
    ('|Halon1211', 'ktHalon1211 / yr', None),
    ('|Halon1202', 'ktHalon1202 / yr', None),
    ('|Halon1301', 'ktHalon1301 / yr', None),
    ('|Halon2402', 'ktHalon2402 / yr', None),
    ('|CH3Br', 'ktCH3Br / yr', None),
    ('|CH3Cl', 'ktCH3Cl / yr', None),
)


def scmdf_to_emissions(scmdf, include_cfcs=True, startyear=1765, endyear=2100):
    """
    Opens an ScmDataFrame and extracts the data. Interpolates linearly
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
            - Use False to ignore and create a 23-species emission file.
        startyear: First year of output file.
        endyear: Last year of output file.

    Returns:
        nt x 40 numpy emissions array (nt x 23 if ``include_cfcs`` is ``False``)
    """

    # We expect that aeneris and silicone are going to give us a nicely
    # formatted ScmDataFrame with all 23 species present and correct at
    # timesteps 2015, 2020 and ten-yearly to 2100.
    # We also implicitly assume that data up until 2014 will follow SSP
    # historical.
    # This adapter will not be tested on anything else!

    n_cols = 40
    nt = endyear - startyear + 1

    data_out = np.ones((nt, n_cols)) * np.nan
    data_out[:,0] = np.arange(startyear, endyear+1)

    if not has_scmdata:
        raise ImportError("This is not going to work without having scmdata installed")

    if not isinstance(scmdf, ScmDataFrame):
        raise TypeError("scmdf must be an scmdata.ScmDataFrame instance")

    if not include_cfcs:
        raise NotImplementedError("include_cfcs equal to False")

    if scmdf[["model", "scenario"]].drop_duplicates().shape[0] != 1:
        raise AssertionError("Should only have one model-scenario pair")

    # fill in 1765 to 2014 from SSP emissions
    ssp_df = ScmDataFrame(os.path.join(os.path.dirname(__file__), '../SSPs/data/rcmip-emissions-annual-means-4-0-0-ssp-only.csv'))

    years = scmdf["year"].values
    first_scenyear = years[0]
    last_scenyear = years[-1]
    first_scen_row = int(first_scenyear-startyear)
    last_scen_row = int(last_scenyear-startyear)

    for i, (specie, unit, context) in enumerate(EMISSIONS_SPECIES_UNITS_CONTEXT):
        data_out[:first_scen_row, i+1] = ssp_df.filter(
            variable="*{}".format(specie),
            region="World",
            scenario="ssp245",
            year=range(startyear, 2015)
        ).convert_unit(unit, context=context).values.squeeze()

        if i < 23:
            if not any([specie in v for v in scmdf.get_unique_meta("variable")]):
                raise AssertionError("{} not available in scmdf".format(specie))

            f = interp1d(
                years,
                scmdf.filter(
                    variable="*{}".format(specie),
                    region="World"
                ).convert_unit(unit, context=context).values.squeeze()
            )
            data_out[first_scen_row:(last_scen_row+1), i+1] = f(
                np.arange(first_scenyear, last_scenyear+1)
            )

        else:
            if not any([specie in v for v in ssp_df.get_unique_meta("variable")]):
                raise AssertionError("{} not available in ssp_df".format(specie))

            filler_data = ssp_df.filter(
                scenario="ssp245",
                variable="*{}".format(specie),
                year=range(2015, 2500 + 1),
            )

            f = interp1d(
                filler_data["year"].values,
                filler_data.convert_unit(unit, context=context).values.squeeze()
            )
            data_out[first_scen_row:(last_scen_row+1), i+1] = f(
                np.arange(first_scenyear, last_scenyear+1)
            )

    return data_out
