from __future__ import print_function
import numpy as np
import cantera as ct
from .constants import (cantera_version,
                        one_atm_in_bar,
                        one_atm_in_torr,
                        one_bar_in_pa,
                        )


class PressureTrace(object):
    """Generic class for pressure traces"""

    def file_loader(self, filename):
        """
        Load a voltage trace from a text file.

        Load a voltage trace from a text file. Check if the file exists
        and if not, try again after adding the proper file extension.
        """
        self.voltage = None
        try:
            self.voltage = np.genfromtxt(filename)
        except OSError:
            filename += '.txt'
            self.voltage = np.genfromtxt(filename)
        if self.voltage is None:
            raise OSError('Data file not found')

    def smoothing(self, pres, span=21):
        window = np.ones(span)/span
        return np.convolve(pres, window, 'same')

    def pressure_fit(self):
        beg_compress = np.floor(self.p_EOC_idx - 0.08*self.sampfreq)
        time = np.linspace(0, (beg_compress - 1)/self.sampfreq, beg_compress)
        fit_pres = self.pressure[:beg_compress]
        fit_pres[0:9] = fit_pres[10]
        self.linear_fit = np.polyfit(time, fit_pres, 1)

    def compress(self):
        self.maxp = np.amax(self.pressure)
        self.maxpi = np.argmax(self.pressure)
        minpi = self.maxpi - 100
        while self.pressure[minpi] >= self.pressure[minpi - 100]:
            minpi -= 1

        self.p_EOC = np.amax(self.pressure[0:minpi])
        self.p_EOC_idx = np.argmax(self.pressure[0:minpi])
        diff = abs(self.pressure[self.p_EOC_idx] - self.pressure[15])
        if diff < 5:
            self.p_EOC, self.p_EOC_idx = self.maxp, self.maxpi

    def derivative(self):
        """

        """
        m = len(self.pressure)
        self.dpdt = np.zeros(m)
        for i in range(m-2):
            self.dpdt[i] = (-self.pressure[i+2] + 4*self.pressure[i+1] -
                            3*self.pressure[i])/(2*(self.time[i+1] -
                                                    self.time[i]))
        self.dpdt[np.isinf(self.dpdt)] = 0


class PressureFromTemperature(PressureTrace):
    """Class for pressure trace computed from a temperature trace."""

    def __init__(self, temperature, P_in):
        """Create a pressure trace given a temperature trace.

        The required method to set the temperature and entropy of the
        Solution to set the state are not implemented, so this method
        is a stub for now.
        """
        pass
        # gas = ct.Solution('species.cti')
        # gas.TP = temperature[0], P_in
        # initial_entropy = gas.entropy_mass
        # self.pressure = np.zeros((len(temperature)))
        # for i, v in enumerate(temperature):
        #     gas.ST = initial_entropy, temperature[i]


class PressureFromVolume(PressureTrace):
    """ Class for pressure trace computed from a volume trace."""

    def __init__(self, volume, P_in, T_in=None):
        gas = ct.Solution('species.cti')
        if cantera_version[2] >= 1 and cantera_version[1] >= 2:
            gas.DP = 1.0/volume[0], P_in*one_bar_in_pa
        elif T_in is None:
            raise OSError
        else:
            gas.TP = T_in, P_in
        initial_volume = gas.volume_mass
        initial_entropy = gas.entropy_mass
        self.pressure = np.zeros((len(volume)))
        for i, v in enumerate(volume):
            gas.SV = initial_entropy, v*initial_volume
            self.pressure[i] = gas.P/one_bar_in_pa


class ReactivePressureTrace(PressureTrace):
    """Class for reactive pressure traces."""

    def __init__(self, filename):
        self.file_loader(filename)

        file_info = ParsedFilename(filename)
        initial_pressure_in_bar = file_info.pin*one_atm_in_bar/one_atm_in_torr
        self.pres = (self.voltage[:, 1] - self.voltage[0, 1])*file_info.factor
        self.pres += initial_pressure_in_bar
        self.time = self.voltage[:, 0]

        self.smoothing()
        self.compress()
        self.derivative()


class NonReactivePressureTrace(PressureTrace):
    """Class for non-reactive pressure traces."""

    def __init__(self, filename):
        self.file_loader(filename)

        file_info = ParsedFilename(filename)
        initial_pressure_in_bar = file_info.pin*one_atm_in_bar/one_atm_in_torr
        self.pres = (self.voltage[:, 1] - self.voltage[0, 1])*file_info.factor
        self.pres += initial_pressure_in_bar
        self.time = self.voltage[:, 0]

        self.smoothing()
        self.compress()


class SimulatedPressureTrace(PressureTrace):
    """Class for pressure traces derived from simulations."""

    def __init__(self, filename='export.csv'):
        self.data = np.genfromtxt(filename, delimiter=',', names=True)

        self.pres = self.data['Pressure_(bar)']
        self.time = self.data['Time_(sec)']

    def derivative(self):
        m = len(self.pres)
        self.dpdt = np.zeros(m)
        for i in range(1, m-2):
            x = self.time[i]
            x_min = self.time[i-1]
            x_plu = self.time[i+1]
            y = self.pressure[i]
            y_min = self.pressure[i-1]
            y_plu = self.pressure[i+1]
            self.dpdt[i] = (y_min*(x - x_plu)/((x_min - x)*(x_min - x_plu)) +
                            y*(2*x - x_min - x_plu)/((x - x_min)*(x - x_plu)) +
                            y_plu*(x - x_min)/((x_plu - x_min)*(x_plu - x)))

