"""Experiment Module"""

# System imports
from pathlib import Path

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from .utilities import parse_file_name, copy
from .traces import (VoltageTrace,
                     ExperimentalPressureTrace,
                     TemperatureFromPressure,
                     )


class Condition(object):
    def __init__(self):
        self.reactive_experiments = []
        self.nonreactive_experiments = []

    def add_experiment(self):
        pass


class Experiment(object):
    """Class containing a single RCM experiment"""

    def __init__(self, file_path=None):
        if file_path is None:
            filename = input('Filename: ')
            try:
                self.file_path = Path(filename).resolve()
            except FileNotFoundError:
                filename += '.txt'
                self.file_path = Path(filename).resolve()
        else:
            self.file_path = file_path.resolve()

        self.experiment_parameters = parse_file_name(self.file_path)
        self.voltage_trace = VoltageTrace(self.file_path)
        self.pressure_trace = ExperimentalPressureTrace(self.voltage_trace,
                                                        self.experiment_parameters['pin'],
                                                        self.experiment_parameters['factor'],
                                                        )

        self.ignition_delay, self.first_stage = self.calculate_ignition_delay()

        try:
            self.T_EOC = self.calculate_EOC_temperature()
        except RuntimeError:
            self.T_EOC = 0
            print('Exception in computing the temperature at EOC')
        self.plot_figures()

        # Copy the relevant information to the clipboard for pasting
        # into a spreadsheet
        copy('\t'.join(map(str, [
            self.experiment_parameters['time_of_day'], self.experiment_parameters['pin'],
            self.experiment_parameters['Tin'], self.pressure_trace.p_EOC, self.ignition_delay,
            self.first_stage, self.T_EOC, self.experiment_parameters['spacers'],
            self.experiment_parameters['shims']])))

    def calculate_ignition_delay(self):

        # offset_points is an offset from the EOC to ensure that if
        # ignition is weak, the peak in dP/dt from the compression
        # stroke is not treated as the ignition event. Define points
        # to start looking for and stop looking for ignition.
        offset_points = 0.002*self.pressure_trace.frequency
        start_point = self.pressure_trace.EOC_idx + offset_points
        end_point = self.pressure_trace.EOC_idx + offset_points + 100000

        # The index of the maximum of the smoothed derivative trace
        # is taken as the point of ignition
        idx_of_ig = np.argmax(self.pressure_trace.smoothed_derivative[start_point:end_point])

        # Take the offset into account when calculating the ignition
        # delay. The index of ignition is already relative to zero,
        # so we don't need to subtract any times. Stored in milliseconds
        ignition_delay = self.pressure_trace.time[idx_of_ig + offset_points]*1000

        try:
            end_first_stage = start_point + idx_of_ig - offset_points
            idx_of_first_stage = np.argmax(
                self.pressure_trace.smoothed_derivative[start_point:end_first_stage]
            )
            first_stage = self.pressure_trace.time[idx_of_first_stage + offset_points]*1000
        except ValueError:
            first_stage = 0.0

        return ignition_delay, first_stage

    def calculate_EOC_temperature(self):

        pres_to_temp_start_idx = self.pressure_trace.EOC_idx - 0.03*self.pressure_trace.frequency
        tempp = self.pressure_trace.pressure[(pres_to_temp_start_idx):(self.p_EOC_idx)]
        temperature_trace = TemperatureFromPressure(tempp, self.Tin)
        return np.amax(temperature_trace.temperature)

    def plot_figures(self):
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
        m = plt.get_current_fig_manager()
        m.window.showMaximized()
