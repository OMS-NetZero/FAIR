Changelog
---------

master
------

(`#84 <https://github.com/OMS-NetZero/FAIR/pull/84>`_) Fix concentration-driven runs for 45-species FaIR

(`#83 <https://github.com/OMS-NetZero/FAIR/pull/83>`_) Support scmdata >= 0.7.1

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
