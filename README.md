[![Build Status](https://travis-ci.org/OMS-NetZero/FAIR.svg?branch=master)](https://travis-ci.org/OMS-NetZero/FAIR)

# FAIR
Finite Amplitude Impulse-Response simple climate-carbon-cycle model 

## Installation
1. Make sure you have Python 2 and pip installed and are working in an environment that uses Python 2 by default
1. From terminal/command prompt `pip install fair`

## Usage
FAIR takes emissions of greenhouse gases, aerosol and ozone precursors, and converts these into greenhouse gas concentrations, radiative forcing and temperature change.

There are two ways to run FAIR:
1. Carbon dioxide emissions only with all other radiative forcings specified externally (same behaviour as FAIR 1.0 with the same keyword arguments);
1. All species included in the RCP emissions datasets, with, optionally, solar and volcanic forcing still specified externally. For convenience, the RCP datasets are provided in the RCP subdirectory and can be imported:

```
from fair.RCPs import rcp85
emissions = rcp85.Emissions.emissions
```

For further information, see the example ipython notebook contained in the GitHub repo at https://github.com/OMS-NetZero/FAIR.

## References:
Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.1: A simple emissions-based impulse response and carbon cycle model, Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2017-266, in review, 2017

Millar, R. J., Nicholls, Z. R., Friedlingstein, P., and Allen, M. R.: A modified impulse-response representation of the global near-surface air temperature and atmospheric concentration response to carbon dioxide emissions, Atmos. Chem. Phys., 17, 7213-7228, https://doi.org/10.5194/acp-17-7213-2017, 2017.

Bibtex references:  
@Article{gmd-2017-266,  
AUTHOR = {Smith, C. J. and Forster, P. M. and Allen, M. and Leach, N. and Millar, R. J. and Passerello, G. A. and Regayre, L. A.},  
TITLE = {FAIR v1.1: A simple emissions-based impulse response and carbon cycle model},  
JOURNAL = {Geoscientific Model Development Discussions},  
VOLUME = {2017},  
YEAR = {2017},  
PAGES = {1--45},  
URL = {https://www.geosci-model-dev-discuss.net/gmd-2017-266/},  
DOI = {10.5194/gmd-2017-266}  
}

@Article{acp-17-7213-2017,  
AUTHOR = {Millar, R. J. and Nicholls, Z. R. and Friedlingstein, P. and Allen, M. R.},  
TITLE = {A modified impulse-response representation of the global near-surface air temperature and atmospheric concentration response to carbon dioxide emissions},  
JOURNAL = {Atmospheric Chemistry and Physics},  
VOLUME = {17},  
YEAR = {2017},  
NUMBER = {11},  
PAGES = {7213--7228},  
URL = {https://www.atmos-chem-phys.net/17/7213/2017/},  
DOI = {10.5194/acp-17-7213-2017}  
}
