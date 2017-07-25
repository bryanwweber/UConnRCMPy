"""
Test module for the traces module
"""
import pytest
import pkg_resources
import os
import platform
from pathlib import Path
from numpy import isclose
from ..traces import VoltageTrace, ExperimentalPressureTrace, AltExperimentalPressureTrace, PressureFromVolume, TemperatureFromPressure, VolumeFromPressure
from ..traces import ONE_ATM_IN_TORR, ONE_ATM_IN_BAR, ONE_BAR_IN_PA


class TestConstants(object):
    def test_one_atm_in_torr(self):
        assert ONE_ATM_IN_TORR == 760.0

    def test_one_atm_in_bar(self):
        assert ONE_ATM_IN_BAR == 1.01325

    def test_one_bar_in_pa(self):
        assert ONE_BAR_IN_PA == 1.0E5


class TestVoltageTrace(object):
    """
    """
    @pytest.fixture(scope='class')
    def voltage_trace(self, request):
        file_path = os.path.join(request.param)
        filename = pkg_resources.resource_filename(__name__, file_path)

        return VoltageTrace(Path(filename))

    @pytest.mark.parametrize('voltage_trace', [
        '00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt',
        'NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt'
    ], indirect=['voltage_trace'])
    def test_repr(self, voltage_trace):
        if platform.system() in ['Darwin', 'Linux']:
            path_str = 'PosixPath'
        else:
            path_str = 'WindowsPath'

        repr_str = "VoltageTrace(file_path={}('{}'))".format(path_str, voltage_trace.file_path)
        assert repr(voltage_trace) == repr_str

    @pytest.mark.parametrize('voltage_trace, frequency', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 2802.6561771647121),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', 1008.4673825016182),
    ], indirect=['voltage_trace'])
    def test_filter_frequency(self, voltage_trace, frequency):
        assert isclose(voltage_trace.filter_frequency, frequency)

    @pytest.mark.parametrize('voltage_trace, frequency', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 1000.0),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', 500.0),
    ], indirect=['voltage_trace'])
    def test_change_filter_frequency(self, voltage_trace, frequency):
        voltage_trace.change_filter_freq(frequency)
        assert isclose(voltage_trace.filter_frequency, frequency)


class TestExperimentalPressureTrace(object):
    """
    """
    @pytest.fixture(scope='class')
    def pressure_trace(self, request):
        file_path = os.path.join(request.param)
        filename = pkg_resources.resource_filename(__name__, file_path)

        v = VoltageTrace(Path(filename))
        return ExperimentalPressureTrace(v, 1146, 100)

    @pytest.mark.parametrize('pressure_trace, frequency, is_reactive, EOC_idx, p_EOC', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 100E3, True, 16702, 30.110942557214699),
    ], indirect=['pressure_trace'])
    def test_create(self, pressure_trace, frequency, is_reactive, EOC_idx, p_EOC):
        assert pressure_trace.frequency == frequency
        assert pressure_trace.is_reactive == is_reactive
        assert pressure_trace.EOC_idx == EOC_idx
        assert isclose(pressure_trace.p_EOC, p_EOC)

    @pytest.mark.parametrize('pressure_trace, is_reactive, p_EOC', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 'True', '30.110942557214699'),
    ], indirect=['pressure_trace'])
    def test_repr(self, pressure_trace, is_reactive, p_EOC):
        repr_str = "ExperimentalPressureTrace(p_EOC={}, is_reactive={})".format(p_EOC, is_reactive)
        assert repr(pressure_trace) == repr_str

    @pytest.mark.parametrize('pressure_trace, linear_fit', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', [-0.00744538, 1.52787853]),
    ], indirect=['pressure_trace'])
    def test_pressure_fit(self, pressure_trace, linear_fit):
        fit = pressure_trace.pressure_fit()
        assert isclose(fit[0], linear_fit[0])
        assert isclose(fit[1], linear_fit[1])

    @pytest.mark.parametrize('pressure_trace, time, p_EOC, EOC_idx', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', -4.0, 14.729642457576805, 16302),
    ], indirect=['pressure_trace'])
    def test_change_EOC_time(self, pressure_trace, time, p_EOC, EOC_idx):
        pressure_trace.change_EOC_time(time)
        assert isclose(pressure_trace.p_EOC, p_EOC)
        assert pressure_trace.EOC_idx == EOC_idx

    @pytest.mark.parametrize('pressure_trace, time', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', -10000.0),
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 10000.0),
    ], indirect=['pressure_trace'])
    def test_change_EOC_time_fail(self, pressure_trace, time):
        with pytest.raises(ValueError):
            pressure_trace.change_EOC_time(time)
