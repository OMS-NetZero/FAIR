|Build Status| |Binder| |Docs Status| |Zenodo| |Codecov| |pypi|

FaIR
====

Finite Amplitude Impulse-Response simple climate model

Requirements
------------

- python 3.6+

Installation
------------

We recommend installing `FaIR` within a virtual environment, such as `virtualenv` or `conda`.

From `pypi`
~~~~~~~~~~~
``pip install fair``

From source
~~~~~~~~~~~
Clone the repository from https://github.com/OMS-NetZero/FAIR.

If you would like to contribute to the development of FaIR, please fork the repository
to your own GitHub account, read the contributing guide, and when ready make a pull
request.


Usage
-----

FaIR can be driven by emissions of greenhouse gases (GHGs) and short-lived forcers
(SLCFs), concentrations of GHGs, or effective radiative forcing (ERF), with different
input methods for different species possible in the same run. If run concentration-
driven, emissions are back-calculated. Custom GHGs and SLCFs can
be defined, and all components are optional allowing experiments such as pulse-response
analyses to single forcers or gathering up non-CO2 species as an aggregate forcing.




For further information, see the example ipython notebooks contained in
the GitHub repo at https://github.com/OMS-NetZero/FAIR and the documentation.


Calibration, constraint, and ensembles
--------------------------------------

FaIR is naive. It will run whatever climate scenario and climate configuration you give
it. If you create a runaway greenhouse effect, a snowball Earth, or violate the laws of
physics, FaIR won't judge you. For simple climate models as for complex, garbage in
leads to garbage out. We provide some emissions scenarios (the CMIP6
SSPs, as produced by rcmip.org), and reasonably sensible default configurations
(as of late 2022) for you to play with. These should not be taken to be best estimates
or up-to-date projections of the historical or future climate if you use them. We've
seen plenty of examples of old versions, or out-of-date calibrations, of FaIR being used in the
literature leading to complaints that FaIR is too cold, or too hot, and/or has eaten Baby Bear's
porridge. Calibrations quickly become out of date as new information comes in, as any
good Bayesian will tell you [note 1].

We have produced IPCC
AR6 Working Group 1 consistent probabilistic ensembles elsewhere. These parameter sets
are calibrated to CMIP6 models, run in a large Monte Carlo ensemble, and constrained
based on six observed or assessed climate metrics. If you're writing a paper using FaIR,
you should use these. Give us time and we'll write this into a paper and release the
code; if you want them now, get in touch and we'll give them to you.

[1] we're moving towards releasing annual updates based on re-calibrations to historical emissions
and temperature time series as they come in.

Citation
--------

If you use FaIR in your work, please cite the following reference depending on the version:

- **v2.0+:** Leach, N. J., Jenkins, S., Nicholls, Z., Smith, C. J., Lynch, J., Cain, M., Walsh, T., Wu, B., Tsutsui, J., and Allen, M. R.: FaIRv2.0.0: a generalized impulse response model for climate uncertainty and future scenario exploration, Geosci. Model Dev., 14, 3007–3036, https://doi.org/10.5194/gmd-14-3007-2021, 2021
- **v1.1-v1.6**: Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.3: A simple emissions-based impulse response and carbon cycle model, Geosci. Model Dev., https://doi.org/10.5194/gmd-11-2273-2018, 2018.
- **v1.0** (or the concept of the state-dependent impulse-response function for CO2): Millar, R. J., Nicholls, Z. R., Friedlingstein, P., and Allen, M. R.: A modified impulse-response representation of the global near-surface air temperature and atmospheric concentration response to carbon dioxide emissions, Atmos. Chem. Phys., 17, 7213-7228, https://doi.org/10.5194/acp-17-7213-2017, 2017.

.. |Build Status| image:: https://github.com/OMS-NetZero/FAIR/actions/workflows/continuous_integration.yml/badge.svg
.. |Binder| image:: https://mybinder.org/badge.svg
   :target: https://mybinder.org/v2/gh/OMS-NetZero/FAIR/master?filepath=notebooks/Example-Usage.ipynb
.. |Docs Status| image:: https://readthedocs.org/projects/fair/badge/?version=latest
   :target: http://fair.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1247898.svg
   :target: https://doi.org/10.5281/zenodo.1247898
.. |Codecov| image:: https://codecov.io/gh/OMS-NetZero/FAIR/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/OMS-NetZero/FAIR
.. |pypi| image:: https://img.shields.io/pypi/v/fair
   :target: https://pypi.org/project/fair/
