"""Custom exceptions used in FaIR."""

# TODO: some of what are currently ValueErrors should be defined as custom
# exceptions. I was being lazy.

# Then I got even lazier and purged the two remaining errors that weren't
# ValueErrors. We'll leave this module as a placeholder.


class DuplicateScenarioError(Exception):
    """Used where scenario label occurs more than once in an emissions/concentration/
    forcing file.
    """


class MetaAfterValueError(Exception):
    """Used where a meta column occurs to the right of time series in an emissions/
    concentration/forcing file.
    """


class MissingColumnError(Exception):
    """Used when a required column is missing in an emissions/concentration/forcing
    file.
    """


class MissingDataError(Exception):
    """Used when a defined scenario or species is not present in an emissions/
    concentration/forcing file.
    """


class MissingRegionError(Exception):
    """Used when a defined region is not present in an emissions/concentration/forcing
    file.
    """


class NonMonotonicError(Exception):
    """Used when the time axis is not monotonically increasing in an emissions/
    concentration/forcing file.
    """


class UnitParseError(Exception):
    """Used when the units defined in an emissions/concentration/forcing file for a
    particular species are unknown to fair.
    """
