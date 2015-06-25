import numpy as np
import cantera as ct
from .constants import one_bar_in_pa


class TemperatureTrace(object):
    """Class for temperature traces."""
    pass


class TemperatureFromPressure(TemperatureTrace):

    def __init__(self, pressure, T_in):
        gas = ct.Solution('species.cti')
        gas.TP = T_in, pressure[0]*one_bar_in_pa
        initial_entropy = gas.entropy_mass
        self.temperature = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p*one_bar_in_pa
            self.temperature[i] = gas.T
