"""Experiment Module"""

# System imports
from pathlib import Path

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml

# Local imports
from .utilities import parse_file_name, copy
from .traces import (VoltageTrace,
                     ExperimentalPressureTrace,
                     TemperatureFromPressure,
                     )


class Condition(object):
    def __init__(self, plotting=True):
        self.reactive_experiments = {}
        self.nonreactive_experiments = {}
        if plotting:
            self.plotting = plotting
            self.all_runs_figure = plt.figure('Reactive Pressure Trace Comparison')
            self.all_runs_axis = self.all_runs_figure.add_subplot(1, 1, 1)
            m = plt.get_current_fig_manager()
            m.window.showMaximized()
            self.nonreactive_figure = None

    def add_experiment(self, file_name=None):
        exp = Experiment(file_name)
        if exp.pressure_trace.is_reactive:
            self.reactive_experiments[exp.experiment_parameters['date']] = exp
            if self.plotting:
                self.plot_reactive_figures(exp)
        else:
            self.nonreactive_experiments[exp.experiment_parameters['date']] = exp
            if self.plotting:
                self.plot_nonreactive_figure(exp)

    def plot_reactive_figures(self, exp):
        # Plot the smoothed pressure and overlay future runs
        self.all_runs_axis.plot(
            exp.pressure_trace.zeroed_time,
            exp.pressure_trace.pressure,
            label=exp.experiment_parameters['date'],
        )

        exp.plot_pressure_trace()

    def load_yaml(self):
        """
        Load the yaml file called `volume-trace.yml` containing the
        information needed to construct the volume trace. All the
        following data are required unless otherwise noted. The format
        of the yaml file is

            variable: value

        * `nonrfile`: File name of the non-reactive pressure trace.
                      Type: String
        * `reacfile`: File name of the reactive pressure trace.
                      Type: String
        * `comptime`: Length of time of the compression stroke.
                      Type: Integer
        * `nonrend`: End time used for the produced volume trace.
                      Type: Integer
        * `reacend`: End time of the output reactive pressure trace.
                      Type: Integer
        * `reacoffs`: Offset in number of points from EOC for the
                      reactive case. Optional, defaults to zero.
                      Type: Integer
        * `nonroffs`: Offset in number of points from EOC for the
                      non-reactive case. Optional, defaults to zero.
                      Type: Integer
        """
        with open('volume-trace.yaml') as yaml_file:
            return yaml.load(yaml_file)

        # self.nonrfile = Path(self.yaml_data['nonrfile'])
        # self.reacfile = Path(self.yaml_data['reacfile'])
        # self.comptime = self.yaml_data['comptime']
        # self.nonrend = self.yaml_data['nonrend']
        # self.reacend = self.yaml_data['reacend']
        # self.reacoffs = self.yaml_data.get('reacoffs', 0)
        # self.nonroffs = self.yaml_data.get('nonroffs', 0)

    def plot_nonreactive_figure(self, exp):
        if self.nonreactive_figure is None:
            self.nonreactive_figure = plt.figure('Non-Reactive Pressure Trace Comparison')
            self.nonreactive_axis = self.nonreactive_figure.add_subplot(1, 1, 1)
            m = plt.get_current_fig_manager()
            m.window.showMaximized()

            yaml_data = self.load_yaml()
            reactive_parameters = parse_file_name(Path(yaml_data['reacfile']))
            reactive_case = self.reactive_experiments[reactive_parameters['date']]
            self.nonreactive_axis.plot(
                reactive_case.pressure_trace.zeroed_time,
                reactive_case.pressure_trace.pressure,
                label=reactive_parameters['date'],
            )

        self.nonreactive_axis.plot(
            exp.pressure_trace.zeroed_time,
            exp.pressure_trace.pressure,
            label=exp.experiment_parameters['date'],
        )


class Experiment(object):
    """Contains all the information of a single RCM experiment.

    Parameters
    ----------
    file_path : :class:`pathlib.Path`, optional
        If an argument is supplied, it should be an instance of
        :class:`~pathlib.Path`. If no argument is supplied, the
        filename is read from the standard input as a string and
        resolved into a :class:`~pathlib.Path` object.

    Attributes
    ----------
    file_path : :class:`pathlib.Path`
        Object storing the file path
    experiment_parameters : :class:`dict`
        Stores the parameters of the experiment parsed from the
        filename by :func:`parse_file_name`.
    voltage_trace : :class:`VoltageTrace`
        Stores the experimental voltage signal and related traces
    pressure_trace : :class:`ExperimentalPressureTrace`
        Stores the experimental pressure trace and its derivative
    ignition_delay : :class:`float`
        The overall ignition delay of the experiment. Will be zero for
        a non-reactive experiment.
    first_stage : :class:`float`
        The first stage ignition delay of the experiment. Will be zero
        for a non-reactive experiment of if there is no first-stage
        ignition.
    T_EOC : :class:`float`
        The temperature estimated at the end of compression
    """

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
        if self.pressure_trace.is_reactive:
            self.ignition_delay, self.first_stage = self.calculate_ignition_delay()

            try:
                self.T_EOC = self.calculate_EOC_temperature()
            except RuntimeError:
                self.T_EOC = 0
                print('Exception in computing the temperature at EOC')
        else:
            self.ignition_delay, self.first_stage, self.T_EOC = 0, 0, 0

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
        tempp = self.pressure_trace.pressure[(pres_to_temp_start_idx):(self.pressure_trace.EOC_idx)]
        temperature_trace = TemperatureFromPressure(tempp, self.experiment_parameters['Tin'])
        return np.amax(temperature_trace.temperature)

    def plot_pressure_trace(self):
        # Plot the smoothed pressure and the smoothed derivative
        # on a new figure every time
        fig2 = plt.figure(self.experiment_parameters['date'])
        ax2 = fig2.add_subplot(1, 1, 1)
        ax2.plot(self.pressure_trace.zeroed_time, self.pressure_trace.pressure)
        ax2.plot(self.pressure_trace.zeroed_time, self.pressure_trace.smoothed_derivative/1000)
        m = plt.get_current_fig_manager()
        m.window.showMaximized()
