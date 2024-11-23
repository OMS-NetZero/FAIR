Changelog
---------

v2.2.2
------

(`#170 <https://github.com/OMS-NetZero/FAIR/pull/170>`_) Fix to the radiative forcing scale factor in aerosol-radiation interactions (`#168 <https://github.com/OMS-NetZero/FAIR/issues/168>`_)

(`#169 <https://github.com/OMS-NetZero/FAIR/pull/169>`_) Removed logging messages which were too annoying and slowed down overrides in all but default cases (`#162 <https://github.com/OMS-NetZero/FAIR/issues/162>`_)

v2.2.1
------

(`#167 <https://github.com/OMS-NetZero/FAIR/pull/167>`_) Add support for python 3.13

(`#163 <https://github.com/OMS-NetZero/FAIR/pull/163>`_) Fixed bug with calibrated constained example (`#161 <https://github.com/OMS-NetZero/FAIR/issues/161>`_)

v2.2.0
------

(`#125 <https://github.com/OMS-NetZero/FAIR/pull/125>`_) Stratospheric water vapour forcing is now a function of concentration or emissions, rather than forcing

(`#158 <https://github.com/OMS-NetZero/FAIR/pull/158>`_) Add `fair` to `conda-forge`

(`#155 <https://github.com/OMS-NetZero/FAIR/pull/155>`_) Implement general CSV filling function

(`#153 <https://github.com/OMS-NetZero/FAIR/pull/153>`_) Add support for python 3.12

(`#134 <https://github.com/OMS-NetZero/FAIR/pull/134>`_) Increase code coverage to 100%

v2.1.4
------

(`#152 <https://github.com/OMS-NetZero/FAIR/pull/152>`_) Fix bug in reading in species_configs CSV file when either ch4 or aci are omitted

(`#150 <https://github.com/OMS-NetZero/FAIR/pull/150>`_) Reformat for black update

(`#149 <https://github.com/OMS-NetZero/FAIR/pull/149>`_) Add nightly tests


v2.1.3
------

(`#147 <https://github.com/OMS-NetZero/FAIR/pull/147>`_) Remove support for python 3.7

(`#146 <https://github.com/OMS-NetZero/FAIR/pull/146>`_) Restarts for ocean heat content working correctly

v2.1.2
------

(`#138 <https://github.com/OMS-NetZero/FAIR/pull/138>`_) Zenodo DOI error when using `pooch` fixed

v2.1.1
------

(`#118 <https://github.com/OMS-NetZero/FAIR/pull/118>`_) Development branch merged in

(`#130 <https://github.com/OMS-NetZero/FAIR/issues/130>`_) Fixed breaking change in xarray 2023.9.0

(`#127 <https://github.com/OMS-NetZero/FAIR/pull/127>`_) Add a calibrated, constrained example to the documentation

(`#121 <https://github.com/OMS-NetZero/FAIR/pull/121>`_) Make `FAIR.ghg_forcing_offset` an attribute

(`#120 <https://github.com/OMS-NetZero/FAIR/pull/120>`_) Fix a weird bug with scipy-1.10 and linear algebra. Remove support for python 3.6, add support for 3.11

v2.1.0
------

(`#112 <https://github.com/OMS-NetZero/FAIR/pull/112>`_) Write docs for v2.1

(`#111 <https://github.com/OMS-NetZero/FAIR/pull/111>`_) Large overhaul of FaIR including adding most features from v2.0.0-alpha, plus species-dependent methane lifetime and new interface

v1.6.4
------

(`#101 <https://github.com/OMS-NetZero/FAIR/pull/101>`_) Add the SSP emissions scenarios as a module interface, like RCPs

v1.6.3
------

(`#94 <https://github.com/OMS-NetZero/FAIR/pull/94>`_) Enable restarts with Geoffroy temperature function and GIR carbon cycle

v1.6.2
------

(`#99 <https://github.com/OMS-NetZero/FAIR/pull/99>`_) Add option to use ozone forcing dependent on N2O concentrations

v1.6.1
------

(`#87 <https://github.com/OMS-NetZero/FAIR/pull/87>`_) Add Geoffroy temperature and prescribed forcing to inverse

(`#84 <https://github.com/OMS-NetZero/FAIR/pull/84>`_) Fix concentration-driven runs for 45-species FaIR

(`#83 <https://github.com/OMS-NetZero/FAIR/pull/83>`_) Support scmdata >= 0.7.1

earlier
-------

(`#82 <https://github.com/OMS-NetZero/FAIR/pull/82>`_) Expand AR6 forcing diagnostics to get aerosol direct forcing by species

(`#78 <https://github.com/OMS-NetZero/FAIR/pull/78>`_) Optimise ``fair.tools.scmdf.scmdf_to_emissions`` a bit

(`#77 <https://github.com/OMS-NetZero/FAIR/pull/77>`_) Fixed ``fair.tools.scmdf.scmdf_to_emissions``

(`#76 <https://github.com/OMS-NetZero/FAIR/pull/76>`_) Added in the GIR carbon cycle as an option and ScmDataFrame reader, required for openscm-runner and iiasa-climate-assessment

(`#69 <https://github.com/OMS-NetZero/FAIR/pull/69>`_) Added in switch to directly specify tropospheric ozone forcing time series and added in AR6 radiative efficiencies, lifetimes and molecular weights

(`#68 <https://github.com/OMS-NetZero/FAIR/pull/68>`_) Attempting to infill CFCs before 1765 raises ValueError

(`#62 <https://github.com/OMS-NetZero/FAIR/pull/62>`_) Add inverse carbon cycle to diagnose emissions from concentrations

(`#60 <https://github.com/OMS-NetZero/FAIR/pull/60>`_) Refactor carbon cycle into separate function

(`#4 <https://github.com/OMS-NetZero/FAIR/issues/4>`_) Make ``iirf_h`` a keyword parameter

(`#57 <https://github.com/OMS-NetZero/FAIR/pull/57>`_) Update to the ensemble generator to use the variance of the input data rather than the standard deviation

(`#52 <https://github.com/OMS-NetZero/FAIR/pull/52>`_) Fixed bug that didn't pick up first year of CO2 concentrations and radiative forcing after restart_in. Changed internal book-keeping of CO2 concentrations to be absolute, rather than anomaly from pre-industrial, to be consistent with other GHGs.

(`#49 <https://github.com/OMS-NetZero/FAIR/issues/49>`_) Change name of model to FaIR in docs, README and example notebooks

(`#47 <https://github.com/OMS-NetZero/FAIR/pull/47>`_) Remove ``requirements.txt`` and add ``.readthedocs.yml`` to put install requirements all in one place
