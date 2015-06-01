# -*- coding: utf-8 -*-
"""
Nonreactive pressure trace processor
"""
# System imports

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from pressure_traces import (compress, copy, file_loader,
                             ParsedFilename, smoothing)


def nonreactive():
    """
    This file compares the nonreactive trace with the corresponding reactive
    trace for the condition. It provides a method of user-set offset of the
    nonreactive trace to allow aligning of the compression stroke of the
    nonreactive and reactive cases.
    """
    nonrfile = input('Non-reactive filename: ')

    # Here plot_reactive is set based on the defined variables at
    # runtime. If plot_reactive exists in the list, we don't need
    # to plot the reactive trace again, so we don't even process it.
    # If the variable doesn't exist, we want to plot the reactive
    # trace.
    try:
        plot_reactive
    except NameError:
        plot_reactive = True
    else:
        plot_reactive = False

    # Process the reactive pressure trace
    if plot_reactive:
        reacfile = input('Reactive filename: ')
        redata = file_loader(reacfile)
        rfile_info = ParsedFilename(reacfile)
        repres = redata[:, 1]*rfile_info.factor + rfile_info.rpin*1.01325/760
        smrepr = smoothing(repres)
        _, rpci = compress(smrepr)
        ztimere = redata[:, 0] - redata[rpci, 0]
        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(ztimere, smrepr)

    # Process the non-reactive pressure trace
    nfile_info = ParsedFilename(nonrfile)
    nrdata = file_loader(nonrfile)
    nrpres = nrdata[:, 1]*nfile_info.nfactor + nfile_info.npin*1.01325/760
    nrsmpr = smoothing(nrpres)
    maxnr = np.amax(nrsmpr)
    maxnri = np.argmax(nrsmpr)
    ztimenr = nrdata[:, 0] - nrdata[maxnri, 0]
    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    ax.plot(ztimenr, nrsmpr)
    m = plt.get_current_fig_manager()
    m.window.showMaximized()
    timeofday = nfile_info.data_date.strftime('%H%M')
    copy('\t'.join(map(str, [
        timeofday, nfile_info.npin, nfile_info.Tin, maxnr, 'NR', 'NR', '',
        nfile_info.spacers, nfile_info.shims
        ])))
