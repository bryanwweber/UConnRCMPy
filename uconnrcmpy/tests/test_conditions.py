"""
Test module for the conditions module
"""
import pytest
from unittest import mock
import os
from ..conditions import Condition


class TestCondition(object):
    """
    """
    def test_create_condition_no_cti(self):
        cti_file = os.path.join(os.path.dirname(__file__), 'species.cti')
        with mock.patch('builtins.input', return_value=cti_file):
            Condition(plotting=False)

    def test_create_condition_with_cti(self):
        cti_file = os.path.join(os.path.dirname(__file__), 'species.cti')
        Condition(cti_file=cti_file, plotting=False)

    def test_add_experiment(self):
        datadir = os.path.dirname(__file__)
        cti_file = os.path.join(datadir, 'species.cti')
        reacfile = os.path.join(datadir, '00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt')
        nonrfile = os.path.join(datadir, 'NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt')
        c = Condition(cti_file=cti_file, plotting=False)
        c.add_experiment(reacfile, copy=False)
        c.add_experiment(nonrfile, copy=False)
