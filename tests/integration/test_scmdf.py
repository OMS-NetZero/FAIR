import os.path

import numpy as np
import numpy.testing as npt
import pytest
from scmdata import ScmRun

from fair.tools.scmdf import scmdf_to_emissions, _get_fair_col_unit_context


scenarios_to_test = ["ssp119", "ssp245", "ssp585"]
scenarios_to_test = ["ssp119"]
SCENARIOS = ScmRun(
    os.path.join(
        os.path.dirname(__file__), "rcmip_scen_ssp_world_emissions.csv"
    ),
    lowercase_cols=True
).filter(scenario=scenarios_to_test)

SSP245_EMMS = ScmRun(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "fair",
        "SSPs",
        "data",
        "rcmip-emissions-annual-means-4-0-0-ssp-only.csv"
    ),
    lowercase_cols=True
).filter(scenario="ssp245")


MODEL_SCEN_DFS = []
for scen_scmdf in SCENARIOS.groupby("scenario"):
    for scen_model_scmdf in scen_scmdf.groupby("model"):
        MODEL_SCEN_DFS.append(scen_model_scmdf)


@pytest.fixture(params=MODEL_SCEN_DFS)
def scen_model_scmdfs(request):
    yield request.param


@pytest.mark.parametrize("startyear,endyear", (
    (1765, 2100),
    (1765, 2150),
    (1850, 2150),
    (1850, 2100),
))
def test_scmdf_to_emissions_all_ssps(scen_model_scmdfs, startyear, endyear):
    res = scmdf_to_emissions(
        scen_model_scmdfs, startyear=startyear, endyear=endyear
    )
    assert not np.isnan(res).any()

    npt.assert_allclose(res[:, 0], range(startyear, endyear + 1))

    for yr in [
        startyear,
        1900,
        1950,
        2014,
        2015,
        2020,
        2050,
        2100,
        scen_model_scmdfs["year"].max()
    ]:
        yr = int(yr)
        row_year = yr - startyear

        for var, idx in (
            ("|CO2|MAGICC Fossil and Industrial", 1),
            ("|CO2|MAGICC AFOLU", 2),
            ("|N2O", 4),
            ("|Sulfur", 5),
            ("|NOx", 8),
            ("|CH4", 3),
            ("|SF6", 23),
            ("|CH3Cl", 39),
        ):
            _, fair_unit, fair_context = _get_fair_col_unit_context(var)

            if idx > 23 or yr < 2015:
                # default filler
                checker = SSP245_EMMS
            else:
                checker = scen_model_scmdfs

            raw_val = checker.filter(
                variable="*{}".format(var),
                year=yr,
                region="World",
            ).convert_unit(fair_unit, context=fair_context).values.squeeze()

            assert raw_val.shape != (0, 0)
            assert not np.isnan(res[row_year, idx])
            npt.assert_allclose(res[row_year, idx], raw_val)
