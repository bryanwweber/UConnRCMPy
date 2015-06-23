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
        pass

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
            self.reactive_trace._EOC_idx -
            self.comptime/1000.0*self.reactive_trace.frequency
        )
        self.stroke_pressure = self.reactive_trace.pressure[
            (self.reactive_start_point):(self.reactive_trace.p_EOC_idx + 1 + self.reacoffs)
        ]
        self.post_pressure = self.nonreactive_trace.pressure[
            (self.nonreactive_trace.p_EOC_idx + self.nonroffs):(self.nonreactive_end_idx + self.nonroffs - self.reacoffs)
        ]
        self.print_pressure = self.reactive_trace.pressure[
            (self.reactive_start_point):(self.reactive_end_idx)
        ]
        self.time = np.linspace(-self.comptime/1000, self.nonrend/1000,
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
        )
        self.volume = np.concatenate((self.stroke_volume,
                                      self.post_volume[1:]))

        self.computed_pressure = PressureFromVolume(
            self.volume,
            self.stroke_pressure[0]*1E5,
            self.reactive_trace.Tin,
        )
        self.plot_figure()

        print('{:.4f}'.format(stroke_pressure[0]))
        copy('{:.4f}'.format(stroke_pressure[0]))

        volout = np.vstack((self.time[::5] + self.comptime/1000, self.volume[::5])).transpose()
        presout = np.vstack((self.time[:len(self.print_pressure):5] + self.comptime/1000, self.print_pressure[::5])).transpose()
        np.savetxt('volume.csv', volout, delimiter=',')
        np.savetxt('Tc__P0__T0_{}K_pressure.txt'.format(
            self.reactive_trace.Tin
            ), presout, delimiter='\t')

    def plot_figure(self):
        self.fig = plt.figure(1)
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax.cla()
        self.ax.plot(self.reactive_trace.ztim, self.reactive_trace.pressure)
        self.ax.plot(self.time[:len(self.print_pressure)], self.print_pressure)
        self.ax.plot(self.time, self.computed_pressure)
        self.ax.plot(self.reactive_trace.ztim, self.reactive_line)
        m = plt.get_current_fig_manager()
        m.window.showMaximized()


with open('volume-trace.yaml') as yaml_file:
    y = yaml.load(yaml_file)

nonrfile = y['nonrfile']
reacfile = y['reacfile']
comptime = y['comptime']
nonrend = y['nonrend']
reacend = y['reacend']
try:
    reacoffs = y['reacoffs']
except KeyError:
    reacoffs = 0

try:
    nonroffs = y['nonroffs']
except KeyError:
    nonroffs = 0

nonrdata = np.fromfile(nonrfile, sep=' ')
nonrdata = np.reshape(nonrdata, (len(nonrdata)/2, 2))
_, _, _, npin, nfactor, _ = filename_parse(nonrfile)
reacdata = np.fromfile(reacfile, sep=' ')
reacdata = np.reshape(reacdata, (len(reacdata)/2, 2))
_, _, reacTin, rpin, rfactor, _ = filename_parse(reacfile)
reactime = reacdata[:, 0]
reacpres = reacdata[:, 1]*rfactor + rpin*1.01325/760
nonrtime = nonrdata[:, 0]
nonrpres = nonrdata[:, 1]*nfactor + npin*1.01325/760

reacsmpr = smoothing(reacpres)
nonrsmpr = smoothing(nonrpres)
reacpc, reacpci = compress(reacsmpr)
nonrpci = np.argmax(nonrsmpr)
reacztim = reactime - reactime[reacpci]
reacfreq = np.rint(1/reactime[1])
nonrfreq = np.rint(1/nonrtime[1])
nonrendidx = nonrend/1000*nonrfreq
reacendidx = reacend/1000*reacfreq
startpoint = comptime/1000*reacfreq

reacp = pressure_fit(reacsmpr, reacpci, reacfreq)
line = np.polyval(reacp, reactime)
stroke_pressure = reacsmpr[(reacpci - startpoint):(reacpci + 1 + reacoffs)]
post_pressure = nonrsmpr[(nonrpci + nonroffs):(nonrpci +
                         nonrendidx + nonroffs - reacoffs)]
print_pressure = reacsmpr[(reacpci - startpoint):(reacpci + reacendidx)]
time = np.arange(-comptime/1000, nonrend/1000, 1/reacfreq)
stroke_volume = pressure_to_volume(stroke_pressure, reacTin)
stroke_temperature = pressure_to_temperature(stroke_pressure, reacTin)
compressed_volume = stroke_volume[-1]
compressed_temperature = stroke_temperature[-1]
post_volume = pressure_to_volume(post_pressure,
                                 compressed_temperature,
                                 compressed_volume)
volume = np.concatenate((stroke_volume, post_volume[1:]))
pressure = volume_to_pressure(volume, stroke_pressure[0]*1E5, reacTin)

fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
ax.cla()
ax.plot(reacztim, reacsmpr)
ax.plot(time[:len(print_pressure)], print_pressure)
ax.plot(time, pressure)
ax.plot(reacztim, line)
m = plt.get_current_fig_manager()
m.window.showMaximized()

print('{:.4f}'.format(stroke_pressure[0]))
copy('{:.4f}'.format(stroke_pressure[0]))

volout = np.vstack((time[::5] + comptime/1000, volume[::5])).transpose()
presout = np.vstack((time[:len(print_pressure):5] + comptime/1000, print_pressure[::5])).transpose()
np.savetxt('volume.csv', volout, delimiter=',')
np.savetxt('Tc__P0__T0_{}K_pressure.txt'.format(reacTin), presout, delimiter='\t')
