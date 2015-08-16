"""
Classes related to pressure traces.
"""

# System imports
from pathlib import Path

# Third party imports
import numpy as np
import cantera as ct
from scipy import signal as sig

# Local imports
from .constants import (cantera_version,
                        one_atm_in_bar,
                        one_atm_in_torr,
                        one_bar_in_pa,
                        )
from .utilities import ParsedFilename

__pdoc__ = {
    'ParsedFilename': None,
    'cantera_version': None,
    'one_atm_in_bar': None,
    'one_atm_in_torr': None,
    'one_bar_in_pa': None,
    'k': None,
}


class ExperimentalPressureTrace(object):
    """Generic class for pressure traces"""

    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, value):
        self._frequency = value

    def file_loader(self):
        """
        Load a signal trace from a text file.
        """
        self.signal = np.genfromtxt(str(self.file_path))

    def smoothing(self, data, span=21):
        """
        Smooth the input `data` using a moving average of width `span`.
        """
        window = np.ones(span)/span
        output = sig.fftconvolve(data, window, mode='same')
        midpoint = (span - 1)/2
        output[:midpoint] = output[midpoint]
        return output

    def filtering(self, data, cutoff_hz=10000):
        """
        Filter the input `data` using a low-pass filter with cutoff at 10 kHz
        """
        nyquist_freq = self.frequency/2.0
        n_taps = 2**14
        low_pass_filter = sig.firwin(
            n_taps,
            cutoff_hz/nyquist_freq,
            window='blackman',
        )
        return sig.fftconvolve(data, low_pass_filter, mode='same')

    def process_pressure_trace(self):
        self.file_loader()

        self.time = self.signal[:, 0]
        """The time loaded from the signal trace."""
        self.frequency = np.rint(1/self.time[1])
        """The sampling frequency of the pressure trace."""

        self.filtered_voltage = self.filtering(self.signal[:, 1])
        self.smoothed_voltage = self.smoothing(self.filtered_voltage)

        initial_pressure_in_bar = self.pin*one_atm_in_bar/one_atm_in_torr
        self.pressure = (self.smoothed_voltage - self.smoothed_voltage[0])
        """The smoothed pressure trace."""
        self.pressure *= self.factor
        self.pressure += initial_pressure_in_bar

        self.find_EOC()
        self.dpdt = self.derivative(self.pressure, self.time)
        """The raw derivative calculated from the smoothed pressure."""
        self.smdp = self.smoothing(self.dpdt, span=151)
        """The smoothed derivative."""
        self.ztim = self.time - self.time[self.p_EOC_idx]
        """A time array where the zero point is at the EOC."""

    def pressure_fit(self):
        """
        Fit a line to the part of the pressure trace before compression
        starts.
        """
        beg_compress = np.floor(self.p_EOC_idx - 0.08*self.frequency)
        time = np.linspace(0, (beg_compress - 1)/self.frequency, beg_compress)
        fit_pres = self.pressure[:beg_compress]
        fit_pres[0:9] = fit_pres[10]
        self.linear_fit = np.polyfit(time, fit_pres, 1)

    def find_EOC(self):
        """
        Find the end of compression point and pressure of the pressure
        trace. If the pressure is close to the initial pressure, assume
        the case is non-reactive and set the pressure at the end of
        compression and the index to the max pressure point.
        """
        self.max_p = np.amax(self.pressure)
        self.max_p_idx = np.argmax(self.pressure)
        min_p_idx = self.max_p_idx - 100
        while self.pressure[min_p_idx] >= self.pressure[min_p_idx - 100]:
            min_p_idx -= 1

        self.p_EOC = np.amax(self.pressure[0:min_p_idx])
        self.p_EOC_idx = np.argmax(self.pressure[0:min_p_idx])
        diff = abs(self.pressure[self.p_EOC_idx] - self.pressure[15])
        if diff < 5:
            self.p_EOC, self.p_EOC_idx = self.max_p, self.max_p_idx

    def derivative(self, dep_var, indep_var):
        """
        Calculate the derivative of the `dep_var` with respect to the
        `indep_var` using a second order forward difference. Set any
        points where the derivative is infinite to zero.
        """
        m = len(dep_var)
        ddt = np.zeros(m)
        for i in range(m-2):
            ddt[i] = (-dep_var[i+2] + 4*dep_var[i+1] -
                      3*dep_var[i])/(2*(indep_var[i+1] -
                                        indep_var[i]))
        ddt[np.isinf(ddt)] = 0
        return ddt


class PressureFromTemperature(ExperimentalPressureTrace):
    """Class for pressure trace computed from a temperature trace."""

    def __init__(self, temperature, P_in):
        """Create a pressure trace given a temperature trace.

        The required method to set the temperature and entropy of the
        Cantera Solution to set the state are not implemented, so this
        class is a stub for now.
        """
        pass
        # gas = ct.Solution('species.cti')
        # gas.TP = temperature[0], P_in
        # initial_entropy = gas.entropy_mass
        # self.pressure = np.zeros((len(temperature)))
        # for i, v in enumerate(temperature):
        #     gas.ST = initial_entropy, temperature[i]


class PressureFromVolume(ExperimentalPressureTrace):
    """ Class for pressure trace computed from a volume trace."""

    def __init__(self, volume, p_initial, T_initial=None):
        """Create a pressure trace given a volume trace.

        Compute a pressure trace given a `volume` trace. Also requires
        inputs of initial pressure `p_initial`, and if Cantera is less
        than version 2.2.1, `T_initial`. If Cantera is greater than or
        equal to version 2.2.1, it is possible to set the state by
        pressure and density, so compute the density as the inverse of
        the initial volume.
        """
        gas = ct.Solution('species.cti')
        if cantera_version[2] >= 1 and cantera_version[1] >= 2:
            gas.DP = 1.0/volume[0], p_initial*one_bar_in_pa
        elif T_initial is None:
            raise OSError
        else:
            gas.TP = T_initial, p_initial
        initial_volume = gas.volume_mass
        initial_entropy = gas.entropy_mass
        self.pressure = np.zeros((len(volume)))
        for i, v in enumerate(volume):
            gas.SV = initial_entropy, v*initial_volume
            self.pressure[i] = gas.P/one_bar_in_pa


class ReactivePressureTrace(ExperimentalPressureTrace, ParsedFilename):
    """Class for reactive pressure traces."""

    def __init__(self, file_path=None):
        """
        Load and process a reactive pressure trace from the data file.
        """
        if file_path is None:
            fname = input('Filename: ')
            try:
                self.file_path = Path(fname).resolve()
            except FileNotFoundError:
                fname += '.txt'
                self.file_path = Path(fname).resolve()
        else:
            self.file_path = file_path.resolve()

        super().__init__(self.file_path)
        self.process_pressure_trace()


class NonReactivePressureTrace(ExperimentalPressureTrace, ParsedFilename):
    """Class for non-reactive pressure traces."""

    def __init__(self, file_path=None):
        """
        Load and process a non-reactive pressure trace from the data
        file.
        """
        if file_path is None:
            fname = input('Non-reactive filename: ')
            try:
                self.file_path = Path(fname).resolve()
            except FileNotFoundError:
                fname += '.txt'
                self.file_path = Path(fname).resolve()
        else:
            self.file_path = file_path.resolve()

        super().__init__(self.file_path)
        self.process_pressure_trace()


class SimulatedPressureTrace(object):
    """Class for pressure traces derived from simulations."""

    def __init__(self, filename='export.csv', data=None):
        """
        Load the pressure trace from the simulation file. The default
        filename is `export.csv`, which can be overridden by passing
        the new filename to the constructor. The data is expected to be
        in csv format with a header row of names. The header for the
        pressure is expected to be `'Pressure_(bar)'` and the header
        for the time is expected to be `'Time_(sec)'`.
        """
        if data is None:
            self.data = np.genfromtxt(filename, delimiter=',', names=True)
        else:
            self.data = data
        """The data from the simulation file."""

        self.pres = self.data['Pressure_(bar)']
        """The simulated pressure trace."""
        self.time = self.data['Time_(sec)']
        """The simulated time trace."""

        self.dpdt = self.derivative(self.pres, self.time)
        """The derivative calculated from the simulated pressure trace."""

    def derivative(self, dep_var, indep_var):
        """
        Calculate the derivative of the `dep_var` with respect to the
        `indep_var`. The derivative is calculated by computing the
        first order Lagrange polynomial fit to the point under
        consideration and its nearest neighbors. The Lagrange
        polynomial is used because of the unequal spacing of the
        simulated data.
        """
        m = len(dep_var)
        ddt = np.zeros(m)
        for i in range(1, m-2):
            x = indep_var[i]
            x_min = indep_var[i-1]
            x_plu = indep_var[i+1]
            y = dep_var[i]
            y_min = dep_var[i-1]
            y_plu = dep_var[i+1]
            ddt[i] = (y_min*(x - x_plu)/((x_min - x)*(x_min - x_plu)) +
                      y*(2*x - x_min - x_plu)/((x - x_min)*(x - x_plu)) +
                      y_plu*(x - x_min)/((x_plu - x_min)*(x_plu - x)))

        return ddt


for k in ExperimentalPressureTrace.__dict__.keys():
    if k != '__init__':
        __pdoc__['PressureFromTemperature.{}'.format(k)] = None
        __pdoc__['PressureFromVolume.{}'.format(k)] = None
        __pdoc__['ReactivePressureTrace.{}'.format(k)] = None
        __pdoc__['NonReactivePressureTrace.{}'.format(k)] = None
        __pdoc__['SimulatedPressureTrace.{}'.format(k)] = None

for k in SimulatedPressureTrace.__dict__.keys():
    if 'SimulatedPressureTrace.{}'.format(k) in __pdoc__:
        del __pdoc__['SimulatedPressureTrace.{}'.format(k)]
