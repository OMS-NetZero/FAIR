"""General constants."""

import numpy as np

#: Time in years for a doubling to occur at a rate of 1% per year.
DOUBLING_TIME_1PCT = np.log(2) / np.log(1.01)

#: Axis of the data arrays referencing time.
TIME_AXIS = 0

#: Axis of the data arrays referencing scenario.
SCENARIO_AXIS = 1

#: Axis of the data arrays referencing climate or specie configuration.
CONFIG_AXIS = 2

#: Axis of emissions, concentration and forcing data arrays representing species.
SPECIES_AXIS = 3

#: Axis of atmopsheric gas box for the gas partition data array in full emissions to
#: concentration mode.
GASBOX_AXIS = 4
