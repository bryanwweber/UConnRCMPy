"""
Module to test value of constants
"""

from ..constants import ONE_ATM_IN_TORR, ONE_ATM_IN_BAR, ONE_BAR_IN_PA


def test_one_atm_in_torr():
    assert ONE_ATM_IN_TORR == 760.0


def test_one_atm_in_bar():
    assert ONE_ATM_IN_BAR == 1.01325


def test_one_bar_in_pa():
    assert ONE_BAR_IN_PA == 1.0E5
