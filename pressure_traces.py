from __future__ import print_function
import numpy as np
from datetime import datetime
import win32clipboard
import cantera as ct


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


def filename_parse(filename):
    filename = filename.lstrip('NR_')
    filename = filename.rstrip('.txt')
    name_split = filename.split('_')
    spacers = int(name_split[0])/10
    shims = int(name_split[2])
    name_split = name_split[4].split('-', maxsplit=3)
    Tin = int(name_split[0][:-1])
    pin = int(name_split[1][:-1])
    factor = int(name_split[2][:-1])
    data_date = datetime.strptime(name_split[3], '%d-%b-%y-%H%M')
    return spacers, shims, Tin, pin, factor, data_date


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
