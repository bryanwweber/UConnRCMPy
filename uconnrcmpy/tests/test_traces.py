"""
Test module for the traces module
"""
import pytest
import pkg_resources
import os
import platform
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
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

    def test_filtering(self):
        nyquist_freq = 1000
        time = np.linspace(0, 1, 2*nyquist_freq+1)
        noise_freq = 5000
        signal = np.sin(2*np.pi*noise_freq*time)
        with patch.object(VoltageTrace, '__init__', lambda x, y: None, None):
            vt = VoltageTrace(None)
            vt.frequency = 2*nyquist_freq
            vt.filter_frequency = 500
            filt_sig = vt.filtering(signal)
        assert all(np.isclose(filt_sig, 0.0))

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
        assert np.isclose(voltage_trace.filter_frequency, frequency)

    @pytest.mark.parametrize('voltage_trace, frequency', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 1000.0),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', 500.0),
    ], indirect=['voltage_trace'])
    def test_change_filter_frequency(self, voltage_trace, frequency):
        voltage_trace.change_filter_freq(frequency)
        assert np.isclose(voltage_trace.filter_frequency, frequency)


class TestExperimentalPressureTrace(object):
    """
    """
    @pytest.fixture(scope='function')
    def voltage_trace(self):
        vt = SimpleNamespace()
        vt.filtered_voltage = np.linspace(0.0, 1.0, 501)
        vt.time = np.arange(0, 501)
        vt.signal = np.vstack((vt.time, vt.filtered_voltage)).T
        vt.frequency = 1
        return vt

    def test_pressure_trace(self, voltage_trace):
        P0 = 760.0  # Torr
        factor = 100
        ep = ExperimentalPressureTrace(voltage_trace, P0, factor)
        pressure = (voltage_trace.filtered_voltage - voltage_trace.filtered_voltage[0])
        pressure *= factor
        pressure += P0/760.0*101325/1.0E5
        assert all(np.isclose(ep.pressure, pressure))
        assert all(np.isclose(ep.raw_pressure, pressure))
        assert all(np.isclose(ep.time, np.arange(0, 501)))

    def test_calculate_derivative(self):
        dep_var = np.linspace(0.0, 1.0, 501)
        indep_var = np.arange(0, 501)
        ddt = ExperimentalPressureTrace.calculate_derivative(None, dep_var, indep_var)
        assert all(np.isclose(ddt[:-2], 0.002))
        assert all(np.isclose(ddt[-2:], 0.0))

    @pytest.fixture(scope='class')
    def pressure_trace(self, request):
        file_path = os.path.join(request.param)
        filename = pkg_resources.resource_filename(__name__, file_path)
        filename = Path(filename)

        v = VoltageTrace(filename)
        spl = filename.name.split('-')
        P0 = float(spl[1].strip('t'))
        factor = int(spl[2].strip('x'))
        return ExperimentalPressureTrace(v, P0, factor)

    @pytest.mark.parametrize('pressure_trace, frequency, is_reactive, EOC_idx, p_EOC', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 100E3, True, 16702, 30.110942557214699),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', 100E3, False, 16546, 30.184215533917),
    ], indirect=['pressure_trace'])
    def test_create(self, pressure_trace, frequency, is_reactive, EOC_idx, p_EOC):
        assert pressure_trace.frequency == frequency
        assert pressure_trace.is_reactive == is_reactive
        assert pressure_trace.EOC_idx == EOC_idx
        assert np.isclose(pressure_trace.p_EOC, p_EOC)

    @pytest.mark.parametrize('pressure_trace, is_reactive, p_EOC', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 'True', '30.110942557214699'),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', 'False', '30.184215533916884'),
    ], indirect=['pressure_trace'])
    def test_repr(self, pressure_trace, is_reactive, p_EOC):
        repr_str = "ExperimentalPressureTrace(p_EOC={}, is_reactive={})".format(p_EOC, is_reactive)
        assert repr(pressure_trace) == repr_str

    @pytest.mark.parametrize('pressure_trace, linear_fit', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', [-0.00744538, 1.52787853]),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', [0.060223508, 1.51117102]),
    ], indirect=['pressure_trace'])
    def test_pressure_fit(self, pressure_trace, linear_fit):
        fit = pressure_trace.pressure_fit()
        assert np.isclose(fit[0], linear_fit[0])
        assert np.isclose(fit[1], linear_fit[1])

    @pytest.mark.parametrize('pressure_trace, time, p_EOC, EOC_idx', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', -4.0, 14.729642457576805, 16302),
        ('NR_00_in_00_mm_333K-1137t-100x-21-Jul-15-1251.txt', -4.0, 15.049653946103636, 16146),
    ], indirect=['pressure_trace'])
    def test_change_EOC_time(self, pressure_trace, time, p_EOC, EOC_idx):
        pressure_trace.change_EOC_time(time)
        assert np.isclose(pressure_trace.p_EOC, p_EOC)
        assert pressure_trace.EOC_idx == EOC_idx

    @pytest.mark.parametrize('pressure_trace, time', [
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', -10000.0),
        ('00_in_00_mm_333K-1146t-100x-21-Jul-15-1226.txt', 10000.0),
    ], indirect=['pressure_trace'])
    def test_change_EOC_time_fail(self, pressure_trace, time):
        with pytest.raises(ValueError):
            pressure_trace.change_EOC_time(time)
