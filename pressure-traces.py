from __future__ import print_function
import numpy as np
import re
import calendar
from datetime import datetime
import matplotlib.pyplot as plt
import sys
import pickle


def compress(pressure):
    maxp = np.max(pressure)
    maxpi = np.argmax(pressure)
    minpi = maxpi - 100
    while pressure[minpi] >= pressure[minpi - 100]:
        minpi -= 1

    pc = np.max(pressure[0:minpi])
    pci = np.argmax(pressure[0:minpi])
    diff = abs(pressure[pci] - pressure[15])
    if diff < 5:
        pc, pci = maxp, maxpi
    return pc, pci, maxp, maxpi


def filename_parse(filename):
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


def smoothing(pressure, span):
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

if __name__ == '__main__':
    filename = '00_in_00_mm_301K-1055t-40x-25-Feb-15-1132'
    data = file_loader(filename)
    spacers, shims, Tin, pin, factor, data_date = filename_parse(filename)
    pres = data[:, 1]*factor + pin*1.01325/760
    time = data[:, 0]
    smpr = smoothing(pres, 20)
    pc, pci, maxp, maxpi = compress(smpr)
    ztim = time - time[pci]
    dpdt = derivative(time, pres)
    smdp = smoothing(dpdt, 5)
    freq = np.rint(1/time[1])
    offt = 0.001
    iofig = np.argmax(smdp[pci+offt*freq:pci+offt*freq+200000])
    timeofig = time[iofig+freq*offt]*1000
    print(timeofig)
    # with open('myplot.pickle', 'rb') as f:
    #     pickle.load(f)
    # print(sys.executable)
    # spacers, shims, Tin, pin, data_date = filename_parse(filename)
    ax = plt.subplot(111)
    plt.plot(ztim, smpr)
    # plt.plot(ztim,smdp/1000)
    with open('myplot.pickle', 'wb') as f:
        pickle.dump(ax, f)
    plt.show()
