"""
Constants used in UConnRCMPy.
"""
import cantera

one_atm_in_torr = 760.0
"""Conversion between atmosphere and Torr."""
one_atm_in_bar = 1.01325
"""Conversion between atmosphere and bar."""
one_bar_in_pa = 1E5
"""Conversion between bar and Pascal."""
cantera_version = None
"""Installed version of Cantera software."""

try:
    cantera_version = map(int, cantera.__version__.split('.'))
except ValueError:
    cantera_version = map(int, cantera.__version__.split('.')[0:2])
    cantera_version.append(int(cantera.__version__.split('.')[2][0]))
