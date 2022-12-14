"""Module for ERFaci tests."""

import numpy as np

from fair.forcing.aerosol.erfaci import logsum


def test_aerosol_erfaci_logsum():
    erfaci = logsum(
        np.array([100, 9, 36]) * np.ones((1, 1, 1, 3)),
        np.array([np.nan] * 3) * np.ones((1, 1, 1, 3)),
        np.array([2.5, 2, 15]) * np.ones((1, 1, 1, 3)),
        np.array([np.nan] * 3) * np.ones((1, 1, 1, 3)),
        np.ones((1, 1, 1, 3)),
        -2.010402,
        np.array([1.058066e-02, 4.426106e-03, 2.779413e-02]),
        np.array([True, True, True]),
        np.array([False, False, False]),
    )
    np.testing.assert_almost_equal(erfaci[0, 0, 0, 0], -1.52353137118739)
