"""
Compare simulated pressure trace to the corresponding reactive trace
"""

# System imports
from glob import glob

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml

# Local imports
from .pressure_traces import copy, derivative


def compare_to_sim(simname='export.csv'):
    """
    Compare a reactive pressure trace to the corresponding simulation.

    This function reads the experimental pressure trace from output in
    the volumetrace function and the output from a CHEMKIN-Pro simulation
    and plots them to verify that the pressure traces match. It copies either
    the simulated ignition delay or the EOC temperature to the clipboard,
    depending on the maximum temperature.
    """

    # Load the simulation data file
    simdata = np.genfromtxt(simname, delimiter=',', skip_header=1)
    simtime = simdata[:, 0]
    # simvolume = simdata[:, 1]
    simtemperature = simdata[:, 2]
    simpressure = simdata[:, 3]

    # Load the experimental pressure trace. Try the glob function first
    # and if it fails, ask the user for help.
    flist = glob('*pressure.txt')
    if not len(flist) == 1:
        flist = [input('Input the experimental pressure trace file name: ')]
    expdata = np.genfromtxt(flist[0])
    exptime = expdata[:, 0]
    exppressure = expdata[:, 1]

    # Plot the pressure traces together
    fig = plt.figure(2)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(exptime, exppressure)
    ax.plot(simtime, simpressure)
    m = plt.get_current_fig_manager()
    m.window.showMaximized()

    # Compute the maximum temperature; if ignition has occurred, judged
    # by maxT > 1200 K, compute the ignition delay. Otherwise, the maxT
    # is Tc. In either case, copy the value to the clipboard.
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
