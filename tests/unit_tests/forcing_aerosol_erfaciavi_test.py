
"""Module for aviation ERFaci tests."""

import numpy as np

from fair.forcing.aerosol.erfaciavi import calculate_erfaciavi_forcing


def test_erfaciavi_linear():

    emissions = (
        np.array([
            3.29469043281061,
            0.2987951381038286,
            0.009671277770680551,
        ])
        * np.ones((1, 1, 1, 3))
    )

    result = calculate_erfaciavi_forcing(
        emissions=emissions,
        forcing_scaling=1.0,
        linear_radiative_efficiency=(
            np.array([
                -0.0033613445378151,
                -0.0629032258064516,
                -0.5799075066509034,
            ])
            * np.ones((1, 1, 1, 3))
        ),
        linear_emissions_indices=np.array([True, True, True]),
    )

    print(result[0, 0, 0, 0])

    np.testing.assert_almost_equal(
        result[0, 0, 0, 0],
        -0.03547821431025829,
    )
    