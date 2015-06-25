"""
Create a volume trace for a given condition
"""
import numpy as np
import matplotlib.pyplot as plt
import yaml
from .pressure_traces import (ReactivePressureTrace,
                              NonReactivePressureTrace,
                              PressureFromVolume,
                              )
from .volume_traces import VolumeFromPressure
from .temperature_traces import TemperatureFromPressure
from .utilities import copy


class VolumeTraceBuilder(object):
    """Class to build a volume trace from an experiment."""

    def __init__(self):
        self.__call__()

    def load_yaml(self):
        with open('volume-trace.yaml') as yaml_file:
            self.yaml_data = yaml.load(yaml_file)
            self.nonrfile = self.yaml_data['nonrfile']
            self.reacfile = self.yaml_data['reacfile']
            self.comptime = self.yaml_data['comptime']
            self.nonrend = self.yaml_data['nonrend']
            self.reacend = self.yaml_data['reacend']
            self.reacoffs = self.yaml_data.get('reacoffs', 0)
            self.nonroffs = self.yaml_data.get('nonroffs', 0)

    def create_volume_trace(self):
        self.load_yaml()
        self.reactive_trace = ReactivePressureTrace(self.reacfile)
        self.nonreactive_trace = NonReactivePressureTrace(self.nonrfile)
        self.reactive_trace.pressure_fit()
        self.reactive_line = np.polyval(self.reactive_trace.linear_fit,
                                        self.reactive_trace.time)

        self.nonreactive_end_idx = (
            self.nonreactive_trace.p_EOC_idx +
            self.nonrend/1000.0*self.nonreactive_trace.frequency
        )
        self.reactive_end_idx = (
            self.reactive_trace.p_EOC_idx +
            self.reacend/1000.0*self.reactive_trace.frequency
        )
        self.reactive_start_point = (
            self.reactive_trace.p_EOC_idx -
            self.comptime/1000.0*self.reactive_trace.frequency
        )
        self.stroke_pressure = self.reactive_trace.pressure[
            (self.reactive_start_point):(self.reactive_trace.p_EOC_idx +
                                         1 + self.reacoffs)
        ]
        self.post_pressure = self.nonreactive_trace.pressure[
            (self.nonreactive_trace.p_EOC_idx + self.nonroffs):(
                self.nonreactive_end_idx + self.nonroffs - self.reacoffs
            )
        ]
        self.print_pressure = self.reactive_trace.pressure[
            (self.reactive_start_point):(self.reactive_end_idx)
        ]
        self.n_print_pts = len(self.print_pressure)
        self.time = np.arange(-self.comptime/1000.0, self.nonrend/1000.0,
                              1/self.reactive_trace.frequency)
        self.stroke_volume = VolumeFromPressure(self.stroke_pressure, 1.0,
                                                self.reactive_trace.Tin).volume
        self.stroke_temperature = TemperatureFromPressure(
            self.stroke_pressure,
            self.reactive_trace.Tin
        ).temperature
        self.post_volume = VolumeFromPressure(
            self.post_pressure,
            self.stroke_volume[-1],
            self.stroke_temperature[-1],
        ).volume
        self.volume = np.concatenate((self.stroke_volume,
                                      self.post_volume[1:]))

        self.computed_pressure = PressureFromVolume(
            self.volume[::5],
            self.stroke_pressure[0]*1E5,
            self.reactive_trace.Tin,
        ).pressure

    def write_output(self):
        volout = np.vstack(
            (self.time[::5] + self.comptime/1000, self.volume[::5])
        ).transpose()
        presout = np.vstack(
            (self.time[:self.n_print_pts:5] + self.comptime/1000,
             self.print_pressure[::5])
        ).transpose()
        np.savetxt('volume.csv', volout, delimiter=',')
        np.savetxt('Tc__P0__T0_{}K_pressure.txt'.format(
            self.reactive_trace.Tin
            ), presout, delimiter='\t')

    def plot_figure(self):
        self.fig = plt.figure('Pressure Trace Comparison')
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax.cla()
        self.ax.plot(self.reactive_trace.ztim, self.reactive_trace.pressure)
        self.ax.plot(self.time[:self.n_print_pts:5], self.print_pressure[::5])
        self.ax.plot(self.time[::5], self.computed_pressure)
        self.ax.plot(self.reactive_trace.ztim, self.reactive_line)
        m = plt.get_current_fig_manager()
        m.window.showMaximized()

    def __call__(self):
        self.create_volume_trace()
        self.plot_figure()

        print('{:.4f}'.format(self.stroke_pressure[0]))
        copy('{:.4f}'.format(self.stroke_pressure[0]))

        self.write_output()
