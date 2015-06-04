# -*- coding: utf-8 -*-
"""
Experimental ignition delay processor
"""
# System imports

# Third party imports
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from .pressure_traces import (copy, compress, derivative, file_loader,
                              ParsedFilename, pressure_to_temperature,
                              smoothing)


def ignitiondelayexp():
    """
    Calculate the ignition delay from an experimental pressure trace

    This function calculates the ignition delay, given an experimental
    pressure trace from an RCM experiment.
    """
    filename = input('Filename: ')
    data = file_loader(filename)
    file_info = ParsedFilename(filename)
    pres = data[:, 1]*file_info.factor + file_info.pin*1.01325/760
    time = data[:, 0]
    smpr = smoothing(pres)
    pc, pci = compress(smpr)
    ztim = time - time[pci]
    freq = np.rint(1/time[1])
    dpdt = derivative(time, smpr)
    smdp = smoothing(dpdt, 5)

    # offt is an offset from the EOC to ensure that if ignition is
    # weak, the peak in dP/dt from the compression stroke is not
    # treated as the ignition event.
    offt = 0.001

    # The index of the maximum of the smoothed derivative trace is taken
    # as the point of ignition
    iofig = np.argmax(smdp[(pci + offt*freq):(pci + offt*freq + 200000)])

    # Computing first stage ignition is disabled at the moment
    firststage = 0.0

    # Take the offset into account when calculating the time of ignition
    timeofig = time[iofig + freq*offt]*1000
    tempp = pres[(pci - 0.03*freq):(pci)]
    temperature = pressure_to_temperature(tempp, file_info.Tin)
    T_EOC = max(temperature)

    # Copy the relevant information to the clipboard for pasting into a
    # spreadsheet
    copy('\t'.join(map(str, [
        file_info.timeofday, file_info.pin, file_info.Tin, pc, timeofig,
        firststage, T_EOC, file_info.spacers, file_info.shims
        ])))

    # Plot the smoothed pressure and overlay future runs
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.plot(ztim, smpr, label=file_info.date)
    m = plt.get_current_fig_manager()
    m.window.showMaximized()

    # Plot the raw and smoothed pressure and the smoothed derivative
    # on a new figure every time
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(ztim, smpr)
    ax2.plot(ztim, smdp/1000)
    ax2.plot(ztim, pres)
    m = plt.get_current_fig_manager()
    m.window.showMaximized()
