Introduction
============

The Finite Amplitude Impulse Response (FaIR) model is a simple
emissions-based climate model. It allows the user to input emissions of
greenhouse gases and short lived climate forcers in order to estimate
global mean atmospheric GHG concentrations, radiative forcing and
temperature anomalies.

The original FaIR model (v1.0) was developed to simulate the earth
system response to CO\ :sub:`2` \ emissions, with all non-CO\ :sub:`2` \
forcing implemented as an "external" source. It was developed by Richard
Millar, Zebedee Nicholls, Pierre Friedlingstein and Myles Allen. The
motivation for developing it and its formulation is documented in a 
`paper published in Atmospheric Chemistry and Physics in 2017
<https://www.atmos-chem-phys.net/17/7213/2017/acp-17-7213-2017.html>`_.

The emissions-based model (v1.3) extends FaIR by replacing all sources of
non-CO\ :sub:`2` \ forcing with relationships that are based on the
source emissions, with the exception of natural forcings (variations 
in solar irradiance and volcanic eruptions). It is useful for
assessing future policy commitments to anthropogenic emissions
(something which we can control) than to radiative forcing (something
which is less certain and which we can only partially control).

If you use this model, please cite the following references:

1. Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.3: A simple emissions-based impulse response and carbon cycle model, Geosci. Model Dev. https://doi.org/10.5194/gmd-2017-266, 2018.

2. Millar, R. J., Nicholls, Z. R., Friedlingstein, P., and Allen, M. R.: A modified impulse-response representation of the global near-surface air temperature and atmospheric concentration response to carbon dioxide emissions, Atmos. Chem. Phys., 17, 7213-7228, https://doi.org/10.5194/acp-17-7213-2017, 2017.

For assistance, `set up an issue on the GitHub site <https://github.com/OMS-NetZero/FAIR/issues>`_ or contact Chris Smith (c.j.smith1 <AT> leeds <DOT> ac <DOT> uk).

Note: The name was modified from FAIR to FaIR in version 1.3.5, to reduce
confusion with an existing evironmental policy decision-making model called
FAIR. If you're looking for that FAIR, `go here
<https://www.pbl.nl/en/publications/2005/The_FAIR_model-a-tool-to-analyse-environmental-and-costs-implications>`_.
