"""Module for unit test of the gas cycle."""

import numpy as np

from fair.gas_cycle.forward import step_concentration
from fair.gas_cycle.inverse import unstep_concentration


def shaped(x):
    return np.asarray(x, dtype=float).reshape(1, 1, 1, 1, -1)


LIFETIME = shaped([1e9, 394.4, 36.54, 4.304])
PARTITION = shaped([0.2173, 0.224, 0.2824, 0.2763])
BOXES = shaped([1000.0, 500.0, 200.0, 50.0])


def step(alpha):
    return step_concentration(
        shaped([40.0]),
        BOXES,
        shaped([BOXES.sum()]),
        shaped([alpha]),
        shaped([278.3]),
        shaped([0.0]),
        shaped([0.1282]),
        LIFETIME,
        PARTITION,
        1,
    )


def test_step_concentration_near_permanent_box_precision():
    box0 = [step(alpha)[1][0, 0, 0, 0, 0] for alpha in (1.0, 1.0 + 1e-9)]
    assert abs(box0[1] - box0[0]) < 1e-12


def test_unstep_concentration_recovers_emissions():
    concentration, boxes, airborne = step(1.0)
    emissions, boxes_back, airborne_back = unstep_concentration(
        concentration,
        BOXES,
        shaped([BOXES.sum()]),
        shaped([1.0]),
        shaped([278.3]),
        shaped([0.0]),
        shaped([0.1282]),
        LIFETIME,
        PARTITION,
        1,
    )
    np.testing.assert_allclose(np.squeeze(emissions), 40.0, rtol=1e-12)
    np.testing.assert_allclose(np.squeeze(boxes_back), np.squeeze(boxes), rtol=1e-12)
    np.testing.assert_allclose(
        np.squeeze(airborne_back), np.squeeze(airborne), rtol=1e-12
    )
