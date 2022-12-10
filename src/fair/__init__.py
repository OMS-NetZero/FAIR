"""Initialiser module for FaIR."""

from . import _version

__version__ = _version.get_versions()["version"]


from .energy_balance_model import EnergyBalanceModel, multi_ebm  # noqa: F401
from .fair import FAIR  # noqa: F401
