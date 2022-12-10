"""Convenience functions for filling FaIR `xarray` instances."""


def fill(var, data, **kwargs):
    """Fill a FaIR variable instance.

    One could work directly with the `xarray` `DataArray`s, but this function
    includes additional error checking and validation, and is slightly less
    cumbersome than the `xarray` method of allocation.

    Parameters
    ----------
    var : `FAIR` variable attribute to fill
        for example fair.climate_configs["ocean_heat_capacity"]
    data : `np.ndarray` like
        data to fill the variable with
    **kwargs : `str`s
        the dimensions represeted in the `data`

    Raises
    ------
    ValueError :
        if a `kwarg` provided doesn't correspond to a dimension name in `var`
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
    """Fill a `fair` variable instance with `value` in first timebound.

    Otherwise identical to `fill`.

    Parameters
    ----------
    var : `FAIR` variable attribute to fill
        for example fair.climate_configs["ocean_heat_capacity"]
    value : `np.ndarray` like
        value to fill the first timebound with
    **kwargs : `str`s
        indices of the dimensions represeted in the `data`
    """
    # check value is a scalar?
    fill(var[0, ...], value, **kwargs)
