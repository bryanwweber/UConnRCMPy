import numpy as np
import cantera as ct

from .constants import cantera_version, one_bar_in_pa


class VolumeTrace(object):
    """Base class for volume traces."""
    pass


class VolumeFromPressure(VolumeTrace):

    def __init__(self, pressure, v_in, T_in):
        gas = ct.Solution('species.cti')
        if cantera_version[2] >= 1 and cantera_version[1] >= 2:
            gas.DP = 1.0/v_in, pressure[0]*one_bar_in_pa
        elif T_in is None:
            raise OSError
        else:
            gas.TP = T_in, pressure[0]*one_bar_in_pa
        initial_entropy = gas.entropy_mass
        initial_density = gas.density
        self.volume = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p*one_bar_in_pa
            self.volume[i] = v_in*initial_density/gas.density
