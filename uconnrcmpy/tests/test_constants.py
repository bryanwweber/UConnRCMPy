"""
Module to test value of constants
"""

from ..constants import one_atm_in_torr, one_atm_in_bar, one_bar_in_pa


def test_one_atm_in_torr():
    assert one_atm_in_torr == 760.0


def test_one_atm_in_bar():
    assert one_atm_in_bar == 1.01325


def test_one_bar_in_pa():
    assert one_bar_in_pa == 1.0E5
