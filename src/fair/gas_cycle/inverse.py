"""Module for deriving emissions from concentration."""

import numpy as np

from ..constants import GASBOX_AXIS


def unstep_concentration(
    concentration,
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
    """
    Calculate emissions from concentrations of any greenhouse gas.

    Parameters
    ----------
    concentration : np.ndarray
        greenhouse gas concentration at the centre of the timestep.
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
        baseline (possibly pre-industrial) concentration of gas in question.
    baseline_emissions : np.ndarray
        baseline (possibly pre-industrial) emissions of gas in question.
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
    boundaries. Therefore it is not actually possible to provide the exact
    emissions that would reproduce the concentrations without using a slower
    root-finding mechanism (that was present in v1.4) and will always be half
    a time step out.

    Where array input is taken, the arrays always have the dimensions of
    (scenario, species, time, gas_box). Dimensionality can be 1, but we
    retain the singleton dimension in order to preserve clarity of
    calculation and speed.

    Returns
    -------
    emissions_out : np.ndarray
        emissions in timestep.
    gas_boxes_new : np.ndarray
        the greenhouse gas atmospheric burden in each lifetime box at the end of
        the timestep.
    airborne_emissions_new : np.ndarray
        airborne emissions (concentrations above pre-industrial control level)
        at the end of the timestep.
    """
    decay_rate = timestep / (alpha_lifetime * lifetime)  # [1]
    decay_factor = np.exp(-decay_rate)  # [1]

    airborne_emissions_new = (
        concentration - baseline_concentration
    ) / concentration_per_emission
    emissions = (
        airborne_emissions_new - np.sum(gasboxes_old * decay_factor, axis=GASBOX_AXIS)
    ) / (
        np.sum(
            partition_fraction / decay_rate * (1.0 - decay_factor) * timestep,
            axis=GASBOX_AXIS,
        )
    )

    gasboxes_new = (
        timestep
        * emissions[..., None]
        * partition_fraction
        * 1
        / decay_rate
        * (1.0 - decay_factor)
        + gasboxes_old * decay_factor
    )
    emissions_out = emissions + baseline_emissions

    return emissions_out, gasboxes_new, airborne_emissions_new
