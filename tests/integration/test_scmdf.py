import os.path

import numpy.testing as npt
import pytest
from scmdata import ScmDataFrame

from fair.tools.scmdf import scmdf_to_emissions


@pytest.fixture(scope="session")
def test_scenarios():
    scenarios = ScmDataFrame(os.path.join(os.path.dirname(__file__), "tiny_example_emissions_harmonized_infilled.csv"))

    scenarios = scenarios.filter(variable="*Infilled*")
    return scenarios


@pytest.mark.parametrize("startyear,endyear", (
    (1765, 2100),
    (1765, 2300),
    (1850, 2300),
    (1850, 2100),
))
def test_scmdf_to_emissions_all_ssps(test_scenarios, startyear, endyear):
    assert test_scenarios["year"].min() == 2015

    for scen_scmdf in test_scenarios.groupby("scenario"):
        for scen_model_scmdf in scen_scmdf.groupby("model"):
            res = scmdf_to_emissions(
                scen_model_scmdf, startyear=startyear, endyear=endyear
            )

            npt.assert_allclose(res[:, 0], range(startyear, endyear + 1))

            for yr in [
                scen_model_scmdf["year"].min(),
                1850,
                1900,
                2015,
                2020,
                2050,
                2100,
                scen_model_scmdf["year"].max()
            ]:
                yr = int(yr)
                row_year = yr - startyear

                for var, idx in (
                    ("|CO2|Energy and Industrial Processes", 1),
                    ("|CO2|AFOLU", 2),
                    ("|CH4", 3),
                    ("|SF6", -1),
                ):
                    raw_val = scen_model_scmdf.filter(
                        variable=var,
                        year=yr,
                        region="World",
                    ).values.squeeze()
                    npt.assert_allclose(res[row_year, idx], raw_val)
