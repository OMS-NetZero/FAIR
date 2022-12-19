.. image:: https://github.com/OMS-NetZero/FAIR/actions/workflows/checks.yml/badge.svg
   :target: https://github.com/OMS-NetZero/FAIR/actions

.. image:: https://mybinder.org/badge.svg
   :target: https://mybinder.org/v2/gh/OMS-NetZero/FAIR/master?filepath=examples/basic_run_example.ipynb

.. image:: https://readthedocs.org/projects/fair/badge/?version=latest
   :target: http://fair.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1247898.svg
   :target: https://doi.org/10.5281/zenodo.1247898

.. image:: https://codecov.io/gh/OMS-NetZero/FAIR/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/OMS-NetZero/FAIR

.. image:: https://img.shields.io/pypi/v/fair
   :target: https://pypi.org/project/fair/


FaIR
====

FaIR (the Finite-amplitude Impulse-Response) climate model is a simple climate model, or *emulator*, useful for producing global mean temperature projections from a wide range of emissions or prescribed forcing scenarios.

Requirements
------------

- python 3.6+


Installation
------------

From the Python Package Index::

    pip install fair

For other options refer to `the documentation <https://fair.readthedocs.io/en/latest/install.html>`_

Usage
-----

FaIR can be driven by emissions of greenhouse gases (GHGs) and short-lived forcers (SLCFs), concentrations of GHGs, or effective radiative forcing (ERF), with different input methods for different species possible in the same run. If run concentration-driven, emissions are back-calculated. Custom GHGs and SLCFs can be defined, and all components are optional allowing experiments such as pulse-response analyses to single forcers or gathering up non-CO\ :sub:`2` species as an aggregate forcing.

Examples
--------

The `examples <examples/>`_ folder contains Jupyter notebooks with some simple examples showing how to run FaIR and the standalone energy balance model.

If you want to try this out online, `go here <https://mybinder.org/v2/gh/OMS-NetZero/FAIR/master?filepath=examples/basic_run_example.ipynb>`_.


Important: A note about calibrating and constraining
----------------------------------------------------

FaIR is naive. It will run whatever climate scenario and climate configuration you give it. If you violate the laws of physics, FaIR won't stop you. For simple climate models as for complex, garbage in leads to garbage out.  More subtle to spot are those analyses with simple climate models where the present day warming (or historical) is wrong or the climate is warming too slowly or too quickly. At least, plot a historical temperature reconstruction over your results and see if it looks right.

We provide some emissions scenarios (the CMIP6 SSPs, as produced by rcmip.org), and reasonably sensible default configurations (as of late 2022) for you to play with. These should not be taken to be best estimates or up-to-date projections of the historical or future climate. Do not interpret the CMIP6 calibrations from the `CMIP6 emissions driven <examples/cmip6_ssp_emissions_run.ipynb>`_ notebook as a true assessment of the uncertainty in future climate (or even of CMIP6 models, as only the climate response is calibrated and not e.g. the aerosol forcing).

We have produced IPCC AR6 Working Group 1 consistent probabilistic ensembles elsewhere in another repository. These parameter sets are calibrated to CMIP6 models, run in a large Monte Carlo ensemble, and constrained based on observed and assessed climate metrics. If you're writing a paper using FaIR, you should use these. They are not public yet, but give us time and we'll write this into a paper and release the code and calibrations. If you want them now, `get in touch <https://homepages.see.leeds.ac.uk/~mencsm/contact.htm>`_ and we'll give them to you.

Citation
--------

If you use FaIR in your work, please cite the following references depending on the version:

- **v2.0+:** Leach, N. J., Jenkins, S., Nicholls, Z., Smith, C. J., Lynch, J., Cain, M., Walsh, T., Wu, B., Tsutsui, J., and Allen, M. R.: FaIRv2.0.0: a generalized impulse response model for climate uncertainty and future scenario exploration, Geosci. Model Dev., 14, 3007–3036, https://doi.org/10.5194/gmd-14-3007-2021, 2021
- **v1.1-v1.6**: Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.3: A simple emissions-based impulse response and carbon cycle model, Geosci. Model Dev., https://doi.org/10.5194/gmd-11-2273-2018, 2018.
- **v1.0** (or the concept of the state-dependent impulse-response function for CO\ :sub:`2`): Millar, R. J., Nicholls, Z. R., Friedlingstein, P., and Allen, M. R.: A modified impulse-response representation of the global near-surface air temperature and atmospheric concentration response to carbon dioxide emissions, Atmos. Chem. Phys., 17, 7213-7228, https://doi.org/10.5194/acp-17-7213-2017, 2017.
