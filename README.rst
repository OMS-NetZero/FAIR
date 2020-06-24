| |Build Status|
| |Binder|
| |Docs Status|
| |Zenodo|
| |Codecov|

FaIR
====

Finite Amplitude Impulse-Response simple climate-carbon-cycle model

Installation
------------

#. Make sure you have Python 3.6+ and pip installed
#. From terminal/command prompt ``pip install fair``

Usage
-----

FaIR takes emissions of greenhouse gases, aerosol and ozone precursors,
and converts these into greenhouse gas concentrations, radiative forcing
and temperature change.

There are two ways to run FaIR:

#. Carbon dioxide emissions only with all other radiative forcings
   specified externally (specify ``useMultigas=False`` in the call to
   ``fair_scm``);
#. All species included in the RCP emissions datasets, with, optionally,
   solar and volcanic forcing still specified externally. For
   convenience, the RCP datasets are provided in the RCP subdirectory
   and can be imported:

::

    from fair.forward import fair_scm
    from fair.RCPs import rcp85
    emissions = rcp85.Emissions.emissions
    C,F,T = fair_scm(emissions=emissions)

The main engine of the model is the ``fair_scm`` function in
``forward.py``. This function can be imported into a Python script or
iPython session. The most important keyword to ``fair_scm`` is the
``emissions``. This should be either a (nt, 40) numpy array (in multigas
mode) or (nt,) numpy array (in CO2 only mode), where nt is the number of
model timesteps. The outputs are a tuple of ``(C, F, T)`` arrays which
are GHG concentrations ((nt, 31) in multigas mode, (nt,) in CO2-only
mode), forcing ((nt, 13) or (nt,)) and temperature change (nt,). The
index numbers corresponding to each species will be given in tables 1 to
3 of the revised version of the Smith et al. paper reference below (we
hope to make this object-oriented in the future). For now, note that the
input emissions follow the ordering of the RCP datasets, which are
included under ``fair/RCPs``, and the GHG concentrations output are in
the same order, except that we don't output the year, only use one
column for total CO2, and the short-lived species (input indices 5 to 11
inclusive) are not included, reducing the number of columns from 40 to
31. In multigas mode the forcing output indices are:

0. CO\ :sub:`2`\
1. CH\ :sub:`4`\
2. N\ :sub:`2`\ O
3. Minor GHGs (CFCs, HFCs etc)
4. Tropospheric ozone
5. Stratospheric ozone
6. Stratospheric water vapour from methane oxidation
7. Contrails
8. Aerosols
9. Black carbon on snow
10. Land use
11. Volcanic
12. Solar


For further information, see the example ipython notebook contained in
the GitHub repo at https://github.com/OMS-NetZero/FAIR.

References:
-----------

Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J.,
Passerello, G. A., and Regayre, L. A.: FAIR v1.3: A simple
emissions-based impulse response and carbon cycle model, Geosci. Model
Dev., https://doi.org/10.5194/gmd-11-2273-2018, 2018.

Millar, R. J., Nicholls, Z. R., Friedlingstein, P., and Allen, M. R.: A
modified impulse-response representation of the global near-surface air
temperature and atmospheric concentration response to carbon dioxide
emissions, Atmos. Chem. Phys., 17, 7213-7228,
https://doi.org/10.5194/acp-17-7213-2017, 2017.

.. |Build Status| image:: https://travis-ci.org/OMS-NetZero/FAIR.svg?branch=master
   :target: https://travis-ci.org/OMS-NetZero/FAIR
.. |Binder| image:: https://mybinder.org/badge.svg
   :target: https://mybinder.org/v2/gh/OMS-NetZero/FAIR/master?filepath=notebooks/Example-Usage.ipynb
.. |Docs Status| image:: https://readthedocs.org/projects/fair/badge/?version=latest
   :target: http://fair.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1247898.svg
   :target: https://doi.org/10.5281/zenodo.1247898
.. |Codecov| image:: https://codecov.io/gh/OMS-NetZero/FAIR/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/OMS-NetZero/FAIR
