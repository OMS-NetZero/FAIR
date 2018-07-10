from __future__ import division

import numpy as np
from itertools import compress
from ..constants.general import M_ATMOS
from ..constants import molwt as molwt_builtin, preindconc,\
   lifetime as lt_builtin


def _lookup(species):
    # get list of GHGs from preindconc model
    named_gases = dir(preindconc)

    # strip out reserved names and the "aslist" list
    include = [i[:2]!='__' and i!='aslist' for i in named_gases]
    named_gases = list(compress(named_gases, include))

    # see if the species appears on the list of builtin gases and return
    # properties if so, otherwise raise error
    if species in named_gases:
        inout = {}
        inout['pi'] = preindconc
        inout['lt'] = lt_builtin
        inout['mw'] = molwt_builtin
        exec("C, L, M = pi."+species+", lt."+species+", mw."+species, inout)
        return inout['C'], inout['L'], inout['M']
    else:
        raise ValueError(species + ' is not in the list of recognised '+
          'greenhouse gases') 


def emissions(C=None, lifetime=None, molwt=None, species=None):
    """Calculate steady state background emissions from a given lifetime

    Keywords:
        C: concentrations, scalar or None
            if scalar, steady state concentrations to acheive.
            if None, use the pre-industrial concentrations for the specified
              gas.
        lifetime: scalar or None
            if scalar, use the given greenhouse gas atmospheric lifetime
            if None, use the default lifetime for the specified gas
        molwt: scalar or None
            if scalar, use molecular weight of given species
            if None, use the default molecular weight for the specified gas
        species: string or None
            Name of the greenhouse gas concentrations you want to calculate.
            See ..constants.lifetime module for used names.
            If None, use the values given in C, lifetime and wt.

    Any user-specified values for C, lifetime and molwt overrides the
    defaults. They must be specified if species is None.

    Returns:
        emissions (scalar)
    """

    # if not using a built-in gas, C, lifetime and molwt must be specified.
    if species is None:
        if any((C is None, lifetime is None, molwt is None)):
            raise ValueError('If species is not given then C, '+
              'lifetime and molwt must be specified.')
    else:
        # populate defaults but override if values given
        species = species.upper()
        C0, lifetime0, molwt0 = _lookup(species)
        if C is None: C=C0
        if lifetime is None: lifetime=lifetime0
        if molwt is None: molwt=molwt0

        # convert units for N2O
        if species == 'N2O':
            molwt = molwt * molwt_builtin.N2/molwt_builtin.N2O

    # now invert the emissions-concentrations relationship
    E = (C * (1.0 - np.exp(-1.0 / lifetime)) * M_ATMOS * molwt /
      molwt_builtin.AIR) * 1e-18

    return E
