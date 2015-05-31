# -*- coding: utf-8 -*-
"""
Created on Thu May 21 17:55:00 2015

@author: weber
"""
# System imports
from glob import glob

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml

# Local imports
from .pressure_traces import copy, derivative


def compare_to_sim():
    simdata = np.genfromtxt('export.csv', delimiter=',', skip_header=1)
    simtime = simdata[:, 0]
    # simvolume = simdata[:, 1]
    simtemperature = simdata[:, 2]
    simpressure = simdata[:, 3]

    flist = glob('*pressure.txt')
    if not len(flist) == 1:
        flist = [input('Input the experimental pressure trace file name: ')]
    expdata = np.genfromtxt(flist[0])
    exptime = expdata[:, 0]
    exppressure = expdata[:, 1]

    fig = plt.figure(2)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(exptime, exppressure)
    ax.plot(simtime, simpressure)
    m = plt.get_current_fig_manager()
    m.window.showMaximized()

    maxT = np.amax(simtemperature)
    if maxT > 1200:
        with open('volume-trace.yaml') as yaml_file:
            y = yaml.load(yaml_file)

        comptime = y['comptime']
        dpdt = derivative(simtime, simpressure)
        ign_delay = simtime[np.argmax(dpdt)]*1000 - comptime
        print('{:.6f}'.format(ign_delay))
        copy(str(ign_delay))
    else:
        print('{:.0f}'.format(maxT))
        copy(str(maxT))
