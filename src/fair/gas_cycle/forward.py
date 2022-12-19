"""Module for the forward (emissions to concentration) model."""

import numpy as np

from ..constants import GASBOX_AXIS


def step_concentration(
    emissions,
    gasboxes_old,
    airborne_emissions_old,
    alpha_lifetime,
    baseline_concentration,
    baseline_emissions,
    concentration_per_emission,
    lifetime,
    partition_fraction,
    timestep,
):
    """Calculate concentrations from emissions of any greenhouse gas.

    Parameters
    ----------
    emissions : np.ndarray
        emissions rate (emissions unit per year) in timestep.
    gas_boxes_old : np.ndarray
        the greenhouse gas atmospheric burden in each lifetime box at the end of
        the previous timestep.
    airborne_emissions_old : np.ndarray
        The total airborne emissions at the beginning of the timestep. This is
        the concentrations above the pre-industrial control. It is also the sum
        of gas_boxes_old if this is an array.
    alpha_lifetime : np.ndarray
        scaling factor for `lifetime`. Necessary where there is a state-
        dependent feedback.
    baseline_concentration : np.ndarray
        original (possibly pre-industrial) concentration of gas(es) in question.
    baseline_emissions : np.ndarray or float
        original (possibly pre-industrial) emissions of gas(es) in question.
    concentration_per_emission : np.ndarray
        how much atmospheric concentrations grow (e.g. in ppm) per unit (e.g.
        GtCO2) emission.
    lifetime : np.ndarray
        atmospheric burden lifetime of greenhouse gas (yr). For multiple
        lifetimes gases, it is the lifetime of each box.
    partition_fraction : np.ndarray
        the partition fraction of emissions into each gas box. If array, the
        entries should be individually non-negative and sum to one.
    timestep : float
        emissions timestep in years.

    Notes
    -----
    Emissions are given in time intervals and concentrations are also reported
    on the same time intervals: the airborne_emissions values are on time
    boundaries and these are averaged before being returned.

    Where array input is taken, the arrays always have the dimensions of
    (scenario, species, time, gas_box). Dimensionality can be 1, but we
    retain the singleton dimension in order to preserve clarity of
    calculation and speed.

    Returns
    -------
    concentration_out : np.ndarray
        greenhouse gas concentrations at the centre of the timestep.
    gas_boxes_new : np.ndarray
        the greenhouse gas atmospheric burden in each lifetime box at the end of
        the timestep.
    airborne_emissions_new : np.ndarray
        airborne emissions (concentrations above pre-industrial control level)
        at the end of the timestep.
    """
    decay_rate = timestep / (alpha_lifetime * lifetime)
    decay_factor = np.exp(-decay_rate)

    # additions and removals
    gasboxes_new = (
        partition_fraction
        * (emissions - baseline_emissions)
        * 1
        / decay_rate
        * (1 - decay_factor)
        * timestep
        + gasboxes_old * decay_factor
    )

    airborne_emissions_new = np.sum(gasboxes_new, axis=GASBOX_AXIS)
    concentration_out = (
        baseline_concentration + concentration_per_emission * airborne_emissions_new
    )

    return concentration_out, gasboxes_new, airborne_emissions_new
