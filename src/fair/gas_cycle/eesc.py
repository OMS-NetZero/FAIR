"""Module for calculaing equivalent effective stratospheric chlorine."""

import numpy as np

from ..constants import SPECIES_AXIS


def calculate_eesc(
    concentration,
    fractional_release,
    cl_atoms,
    br_atoms,
    cfc_11_index,
    halogen_indices,
    br_cl_ratio,
):
    """Calculate equivalent effective stratospheric chlorine.

    Parameters
    ----------
    concentration : np.ndarray
        concentrations in timestep
    fractional_release : np.ndarray
        fractional release describing the proportion of available ODS that
        actually contributes to ozone depletion.
    cl_atoms : np.ndarray
        Chlorine atoms in each species
    br_atoms : np.ndarray
        Bromine atoms in each species
    cfc_11_index : int or None
        array index along SPECIES_AXIS corresponding to CFC-11.
    halogen_indices : list of int
        provides a mapping of which halogen species corresponds to which
        index along the SPECIES_AXIS.
    br_cl_ratio : float
        how much more effective bromine is as an ozone depletor than chlorine.

    Returns
    -------
    eesc_out : np.ndarray
        equivalent effective stratospheric chlorine

    Notes
    -----
    At present, CFC-11 must be provided in the scenario, with species type cfc-11.
    """
    # EESC is in terms of CFC11-eq
    cfc11_fr = fractional_release[:, :, :, cfc_11_index]
    eesc_out = np.nansum(
        (
            cl_atoms * concentration * fractional_release / cfc11_fr
            + br_cl_ratio * br_atoms * concentration * fractional_release / cfc11_fr
        )
        * cfc11_fr,
        axis=SPECIES_AXIS,
        keepdims=True,
    )
    return eesc_out
