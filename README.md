[![image](https://github.com/OMS-NetZero/FAIR/actions/workflows/checks.yml/badge.svg)](https://github.com/OMS-NetZero/FAIR/actions)
[![image](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/OMS-NetZero/FAIR/master?filepath=examples/basic_run_example.ipynb)
[![Documentation Status](https://readthedocs.org/projects/fair/badge/?version=latest)](http://fair.readthedocs.io/en/latest/?badge=latest)
[![image](https://zenodo.org/badge/DOI/10.5281/zenodo.1340643.svg)](https://doi.org/10.5281/zenodo.1340643)
[![image](https://codecov.io/gh/OMS-NetZero/FAIR/branch/master/graph/badge.svg)](https://codecov.io/gh/OMS-NetZero/FAIR)
[![image](https://img.shields.io/pypi/v/fair)](https://pypi.org/project/fair/) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/fair/badges/version.svg)](https://anaconda.org/conda-forge/fair)

# FaIR

FaIR (the Finite-amplitude Impulse-Response) climate model is a simple
climate model, or *emulator*, useful for producing global mean
temperature projections from a wide range of emissions or prescribed
forcing scenarios.

## Requirements

-   python 3.8, 3.9, 3.10, 3.11 or 3.12

## Installation

### From anaconda (recommended)

**NEW!** from v2.1.4, `fair` is available on `conda-forge`:

    conda install -c conda-forge fair

Older versions of `fair` (1.6.2+, 2.1.0-4) can be installed from the `chrisroadmap` channel:

    conda install -c chrisroadmap fair==X.Y.Z

### From the Python Package Index

    pip install fair

### From source

Refer to [the
documentation](https://fair.readthedocs.io/en/latest/install.html)

## Usage

FaIR can be driven by emissions of greenhouse gases (GHGs) and
short-lived forcers (SLCFs), concentrations of GHGs, or effective
radiative forcing (ERF), with different input methods for different
species possible in the same run. If run concentration-driven, emissions
are back-calculated. Custom GHGs and SLCFs can be defined, and all
components are optional allowing experiments such as pulse-response
analyses to single forcers or gathering up non-CO~2~ species as an
aggregate forcing.

## Examples

The examples directory contains Jupyter notebooks with some
simple examples showing how to run FaIR and the standalone energy
balance model.

If you want to try this out online, [go
here](https://mybinder.org/v2/gh/OMS-NetZero/FAIR/master?filepath=examples/basic_run_example.ipynb).

## Important: A note about calibrating and constraining

FaIR is naive. It will run whatever climate scenario and climate
configuration you give it. If you violate the laws of physics, FaIR
won\'t stop you. For simple climate models as for complex, garbage in
leads to garbage out. More subtle to spot are those analyses with simple
climate models where the present day warming (or historical) is wrong or
the climate is warming too slowly or too quickly. At least, plot a
historical temperature reconstruction over your results and see if it
looks right.

We have produced IPCC AR6 Working Group 1 consistent probabilistic
ensembles to run with. The calibration data can be obtained
[here](https://doi.org/10.5281/zenodo.7112539). These parameter sets are
calibrated to CMIP6 models, run in a large Monte Carlo ensemble, and
constrained based on observed and assessed climate metrics. For an
example of how to use this calibration data set with SSP emissions, see
[this
example](https://docs.fairmodel.net/en/latest/examples/calibrated_constrained_ensemble.html).
If you\'re writing a paper using FaIR, you should use these. A paper describing this method has been submitted, but for now please cite the Zenodo DOI.

## Citation

If you use FaIR in your work, please cite the following references
depending on the version:

-   **v2.0+:** Leach, N. J., Jenkins, S., Nicholls, Z., Smith, C. J.,
    Lynch, J., Cain, M., Walsh, T., Wu, B., Tsutsui, J., and Allen, M.
    R.: FaIRv2.0.0: a generalized impulse response model for climate
    uncertainty and future scenario exploration, Geosci. Model Dev., 14,
    3007--3036, <https://doi.org/10.5194/gmd-14-3007-2021>, 2021
-   **v1.1-v1.6**: Smith, C. J., Forster, P. M., Allen, M., Leach, N.,
    Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.3: A
    simple emissions-based impulse response and carbon cycle model,
    Geosci. Model Dev.,
    <https://doi.org/10.5194/gmd-11-2273-2018>, 2018.
-   **v1.0** (or the concept of the state-dependent impulse-response
    function for CO<sub>2</sub>): Millar, R. J., Nicholls, Z. R., Friedlingstein,
    P., and Allen, M. R.: A modified impulse-response representation of
    the global near-surface air temperature and atmospheric
    concentration response to carbon dioxide emissions, Atmos. Chem.
    Phys., 17, 7213-7228,
    <https://doi.org/10.5194/acp-17-7213-2017>, 2017.
