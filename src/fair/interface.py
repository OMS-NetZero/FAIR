"""Convenience functions for filling FaIR `xarray` instances."""


def fill(var, data, **kwargs):
    """Fill a ``FAIR`` variable instance.

    One could work directly with the ``xarray`` DataArrays, but this function
    includes additional error checking and validation, and is slightly less
    cumbersome than the ``xarray`` method of allocation.

    Parameters
    ----------
    var : attr
        ``FAIR`` variable attribute to fill, for example
        fair.climate_configs["ocean_heat_capacity"]
    data : np.ndarray
        data to fill the variable with
    **kwargs :
        the dimensions represeted in ``data``

    Raises
    ------
    ValueError
        if a ``kwargs`` element provided doesn't correspond to a dimension name in
        ``var``
    """
    # does variable exist?
    for kwarg in kwargs:
        if kwarg not in var.coords:
            raise ValueError(
                f"{kwarg} is not a coordinate of {var.name}. Valid coordinates are "
                f"{var.coords._names}"
            )

    var.loc[kwargs] = data


def initialise(var, value, **kwargs):
    """Fill a ``fair`` variable instance with ``value`` in first ``timebound``.

    Otherwise identical to ``fill``.

    Parameters
    ----------
    var : attr
        ``FAIR`` variable attribute to fill for example
        fair.climate_configs["ocean_heat_capacity"]
    value : np.ndarray
        value to fill the first ``timebound`` with
    **kwargs :
        the dimensions represeted in ``data``
    """
    # check value is a scalar?
    fill(var[0, ...], value, **kwargs)
