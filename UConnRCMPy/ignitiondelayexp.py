# -*- coding: utf-8 -*-
"""
Experimental ignition delay processor
"""
# System imports

# Third party imports
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from .pressure_traces import ReactivePressureTrace
from .temperature_traces import TemperatureFromPressure
from .utilities import copy


class ExperimentalIgnitionDelay(ReactivePressureTrace):
    """Class storing one experimental ignition delay case."""

    def __init__(self):
        # First, initialize the pressure trace by calling the __init__
        # of the ReactivePressureTrace class we're inheriting from
        super().__init__()

        # offset_points is an offset from the EOC to ensure that if
        # ignition is weak, the peak in dP/dt from the compression
        # stroke is not treated as the ignition event. Define points
        # to start looking for and stop looking for ignition.
        offset_points = 0.005*self.frequency
        start_point = self.p_EOC_idx + offset_points
        end_point = self.p_EOC_idx + offset_points + 100000

        # The index of the maximum of the smoothed derivative trace
        # is taken as the point of ignition
        self.idx_of_ig = np.argmax(self.smdp[start_point:end_point])

        # Take the offset into account when calculating the ignition
        # delay. The index of ignition is already relative to zero,
        # so we don't need to subtract any times
        self.ignition_delay = self.time[self.idx_of_ig + offset_points]*1000
        self.first_stage = 0.0

        pres_to_temp_start_idx = self.p_EOC_idx - 0.03*self.frequency
        tempp = self.pressure[(pres_to_temp_start_idx):(self.p_EOC_idx)]
        temperature_trace = TemperatureFromPressure(tempp, self.Tin)
        self.T_EOC = np.amax(temperature_trace.temperature)

        # Copy the relevant information to the clipboard for pasting
        # into a spreadsheet
        copy('\t'.join(map(str, [
            self.time_of_day, self.pin, self.Tin, self.p_EOC,
            self.ignition_delay, self.first_stage, self.T_EOC, self.spacers,
            self.shims])))

        # Plot the smoothed pressure and overlay future runs
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(1, 1, 1)
        ax1.plot(self.ztim, self.pressure, label=self.date)
        m = plt.get_current_fig_manager()
        m.window.showMaximized()

        # Plot the raw and smoothed pressure and the smoothed derivative
        # on a new figure every time
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1, 1, 1)
        ax2.plot(self.ztim, self.pressure)
        ax2.plot(self.ztim, self.smdp/1000)
        ax2.plot(self.ztim, self.pres)
        m = plt.get_current_fig_manager()
        m.window.showMaximized()
