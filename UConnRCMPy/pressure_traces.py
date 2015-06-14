from __future__ import print_function
import numpy as np
from datetime import datetime
import win32clipboard
import cantera as ct

one_atm_in_torr = 760.0
one_atm_in_bar = 1.01325
try:
    cantera_version = map(int, ct.__version__.split('.'))
except ValueError:
    cantera_version = map(int, ct.__version__.split('.')[0:2])
    cantera_version.append(int(ct.__version__.split('.')[2][0]))


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

    def smoothing(self, span=21):
        window = np.ones(span)/span
        self.smooth = np.convolve(self.pressure, window, 'same')

    def pressure_to_temperature(self):
        gas = ct.Solution('species.cti')
        gas.TP = self.T_in, self.pressure[0]*1E5
        initial_entropy = gas.entropy_mass
        self.temperature = np.zeros((len(self.pressure)))
        for i, p in enumerate(self.pressure):
            gas.SP = initial_entropy, p*1E5
            self.temperature[i] = gas.T

    def pressure_to_volume(self, V_initial=1):
        gas = ct.Solution('species.cti')
        gas.TP = self.T_in, self.pressure[0]*1E5
        initial_entropy = gas.entropy_mass
        initial_density = gas.density
        self.volume = np.zeros((len(self.pressure)))
        for i, p in enumerate(self.pressure):
            gas.SP = initial_entropy, p*1E5
            self.volume[i] = V_initial*initial_density/gas.density

    def volume_to_pressure(self, P_in):
        gas = ct.Solution('species.cti')
        gas.TP = self.T_in, P_in
        initial_volume = gas.volume_mass
        initial_entropy = gas.entropy_mass
        pressure = np.zeros((len(self.volume)))
        for i, v in enumerate(self.volume):
            gas.SV = initial_entropy, v*initial_volume
            pressure[i] = gas.P/1E5
        return pressure

    def pressure_fit(pressure, pci, sampfreq):
        beg_compress = np.floor(pci - 0.08*sampfreq)
        pressure[:9] = pressure[10]
        time = np.linspace(0, (beg_compress - 1)/sampfreq, beg_compress)
        fit_pres = pressure[:beg_compress]
        polyn = np.polyfit(time, fit_pres, 1)
        return polyn


class ReactivePressureTrace(PressureTrace):
    pass


class NonReactivePressureTrace(PressureTrace):
    pass


def compress(pressure):
    maxp = np.amax(pressure)
    maxpi = np.argmax(pressure)
    minpi = maxpi - 100
    while pressure[minpi] >= pressure[minpi - 100]:
        minpi -= 1

    pc = np.amax(pressure[0:minpi])
    pci = np.argmax(pressure[0:minpi])
    diff = abs(pressure[pci] - pressure[15])
    if diff < 5:
        pc, pci = maxp, maxpi
    return pc, pci


class ParsedFilename(object):
    """
    Parse a filename for experimental conditions.

    Given a filename of the form
    [NR_]XX_in_YY_mm_ZZZK-AAAAt-BBBx-Mon-DD-YY-Time.txt,
    parse the filename for the experimental conditions.
    XX = inches of spacers
    YY = millimeters of shims
    ZZZ = Initial temperature in Kelvin
    AAAA = Initial pressure in Torr
    BBB = Multiplication factor set on the charge amplifier
    Mon-DD-YY-Time = Month, day, year, and time of experiment
    """

    def __init__(self, filename):
        self.fname = filename.lstrip('NR_')
        self.fname = self.fname.rstrip('.txt')
        self.name_split = self.fname.split('_')
        self.spacers = int(self.name_split[0])/10
        self.shims = int(self.name_split[2])
        self.self.name_split = self.name_split[4].split('-', maxsplit=3)
        self.Tin = int(self.name_split[0][:-1])
        self.pin = int(self.name_split[1][:-1])
        self.factor = int(self.name_split[2][:-1])
        self.data_date = datetime.strptime(self.name_split[3], '%d-%b-%y-%H%M')
        self.timeofday = self.data_date.strftime('%H%M')
        self.date = self.data_date.strftime('%d-%b-%H%M')


def file_loader(filename):
    data = None
    try:
        data = np.genfromtxt(filename)
    except OSError:
        filename += '.txt'
        data = np.genfromtxt(filename)
    if data is not None:
        return data
    else:
        raise OSError('Data file not found')


def smoothing(pressure, span=21):
    window = np.ones(span)/span
    smooth = np.convolve(pressure, window, 'same')
    return smooth


def derivative(time, pressure):
    """

    """
    m = len(pressure)
    dpdt = np.zeros(m)
    for i in range(m-2):
        dpdt[i] = (-pressure[i+2] + 4*pressure[i+1] -
                   3*pressure[i])/(2*(time[i+1]-time[i]))
    dpdt[np.isinf(dpdt)] = 0
    dpdt[-1:] = dpdt[-2]
    return dpdt


def copy(text):
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardText(text)
    win32clipboard.CloseClipboard()


def pressure_to_temperature(pressure, T_in):
    gas = ct.Solution('species.cti')
    gas.TP = T_in, pressure[0]*1E5
    initial_entropy = gas.entropy_mass
    temperature = np.zeros((len(pressure)))
    for i, p in enumerate(pressure):
        gas.SP = initial_entropy, p*1E5
        temperature[i] = gas.T
    return temperature


def pressure_to_volume(pressure, T_in, V_initial=1):
    gas = ct.Solution('species.cti')
    gas.TP = T_in, pressure[0]*1E5
    initial_entropy = gas.entropy_mass
    initial_density = gas.density
    volume = np.zeros((len(pressure)))
    for i, p in enumerate(pressure):
        gas.SP = initial_entropy, p*1E5
        volume[i] = V_initial*initial_density/gas.density
    return volume


def volume_to_pressure(volume, P_in, T_in):
    gas = ct.Solution('species.cti')
    gas.TP = T_in, P_in
    initial_volume = gas.volume_mass
    initial_entropy = gas.entropy_mass
    pressure = np.zeros((len(volume)))
    for i, v in enumerate(volume):
        gas.SV = initial_entropy, v*initial_volume
        pressure[i] = gas.P/1E5
    return pressure


def pressure_fit(pressure, pci, sampfreq):
    beg_compress = np.floor(pci - 0.08*sampfreq)
    pressure[:9] = pressure[10]
    time = np.linspace(0, (beg_compress - 1)/sampfreq, beg_compress)
    fit_pres = pressure[:beg_compress]
    polyn = np.polyfit(time, fit_pres, 1)
    return polyn
