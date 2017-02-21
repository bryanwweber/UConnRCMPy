"""
Test module for the experiments module
"""
import numpy as np
import os
import pytest
from ..experiments import Experiment


@pytest.fixture(scope='module')
def files():
    datadir = os.path.dirname(__file__)
    cti_file = os.path.join(datadir, 'species.cti')
    reacfile = os.path.join(datadir, '00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt')
    nonrfile = os.path.join(datadir, 'NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt')
    return {'cti_file': cti_file, 'reacfile': reacfile, 'nonrfile': nonrfile}


def test_create_reactive_experiment(files):
    exp = Experiment(files['reacfile'],
                     cti_file=files['cti_file'], copy=False)
    assert np.isclose(exp.ignition_delay, 65.81)
    assert np.isclose(exp.T_EOC, 765.374)
    assert np.isclose(exp.first_stage, 63.80)
    assert np.isclose(exp.pressure_trace.p_EOC, 30.111)


def test_create_nonreactive_experiment(files):
    exp = Experiment(files['nonrfile'],
                     cti_file=files['cti_file'], copy=False)
    assert np.isclose(exp.ignition_delay, 0.0)
    assert np.isclose(exp.T_EOC, 0.0)
    assert np.isclose(exp.first_stage, 0.0)
    assert np.isclose(exp.pressure_trace.p_EOC, 30.184)
