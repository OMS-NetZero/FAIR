"""Custom exceptions used in FaIR."""

# TODO: some of what are currently ValueErrors should be defined as custom
# exceptions. I was being lazy.

# Then I got even lazier and purged the two remaining errors that weren't
# ValueErrors. We'll leave this module as a placeholder.

class MetaAfterValueError(Exception):
    pass

class MissingColumnError(Exception):
    pass

class MissingDataError(Exception):
    pass

class TimeNotMonotonicError(Exception):
    pass
