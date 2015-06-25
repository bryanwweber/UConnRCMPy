"""
Create a volume trace for a given condition
"""
# System imports

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml

# Local imports
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
        """
        Create a volume trace corresponding to a given experiment.
        The volume trace will be output in csv format into a file
        called `volume.csv`. In addition, the experimental (reactive)
        pressure trace will be output into a file with blanks in the
        name for the user to fill in the data.
        """
        self.__call__()

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
            self.yaml_data = yaml.load(yaml_file)
            self.nonrfile = self.yaml_data['nonrfile']
            self.reacfile = self.yaml_data['reacfile']
            self.comptime = self.yaml_data['comptime']
            self.nonrend = self.yaml_data['nonrend']
            self.reacend = self.yaml_data['reacend']
            self.reacoffs = self.yaml_data.get('reacoffs', 0)
            self.nonroffs = self.yaml_data.get('nonroffs', 0)

    def create_volume_trace(self):
        """
        Create the volume trace based on the information in the loaded
        yaml file.
        """
        self.load_yaml()
        self.reactive_trace = ReactivePressureTrace(self.reacfile)
        """The reactive pressure trace from the experiment."""
        self.nonreactive_trace = NonReactivePressureTrace(self.nonrfile)
        """The non-reactive pressure trace from the experiment."""
        self.reactive_trace.pressure_fit()
        self.reactive_line = np.polyval(self.reactive_trace.linear_fit,
                                        self.reactive_trace.time)
        """
        The line fit through the initial section of the reactive
        pressure trace.
        """

        self.nonreactive_end_idx = (
            self.nonreactive_trace.p_EOC_idx +
            self.nonrend/1000.0*self.nonreactive_trace.frequency
        )
        """
        The index of the desired end point of the non-reactive pressure
        trace.
        """
        self.reactive_end_idx = (
            self.reactive_trace.p_EOC_idx +
            self.reacend/1000.0*self.reactive_trace.frequency
        )
        """
        The index of the desired end point of the reactive pressure
        trace.
        """
        self.reactive_start_point = (
            self.reactive_trace.p_EOC_idx -
            self.comptime/1000.0*self.reactive_trace.frequency
        )
        """
        The index of the desired start point of the reactive pressure
        trace.
        """
        self.stroke_pressure = self.reactive_trace.pressure[
            (self.reactive_start_point):(self.reactive_trace.p_EOC_idx +
                                         1 + self.reacoffs)
        ]
        """
        The pressure during the compression stroke from the reactive
        pressure trace."""
        self.post_pressure = self.nonreactive_trace.pressure[
            (self.nonreactive_trace.p_EOC_idx + self.nonroffs):(
                self.nonreactive_end_idx + self.nonroffs - self.reacoffs
            )
        ]
        """
        The pressure after the compression stroke from the non-reactive
        pressure trace. The end time is set by the yaml variable
        `nonrend`.
        """
        self.print_pressure = self.reactive_trace.pressure[
            (self.reactive_start_point):(self.reactive_end_idx)
        ]
        """
        The pressure trace that will be printed to the output file,
        starting at the beginning of compression and ending at the
        `reacend` time.
        """
        self.n_print_pts = len(self.print_pressure)
        """Number of points in the printed pressure trace."""
        self.time = np.arange(-self.comptime/1000.0, self.nonrend/1000.0,
                              1/self.reactive_trace.frequency)
        """Time array for the output traces."""
        self.stroke_volume = VolumeFromPressure(self.stroke_pressure, 1.0,
                                                self.reactive_trace.Tin).volume
        """
        Volume trace for the compression stroke computed from the
        `stroke_pressure` pressure trace.
        """
        self.stroke_temperature = TemperatureFromPressure(
            self.stroke_pressure,
            self.reactive_trace.Tin
        ).temperature
        """
        Temperature trace for the compression stroke computed from the
        `stroke_pressure` pressure trace.
        """
        self.post_volume = VolumeFromPressure(
            self.post_pressure,
            self.stroke_volume[-1],
            self.stroke_temperature[-1],
        ).volume
        """
        Volume trace for the post compression time computed from
        the `post_pressure` array and the temperature at the end of
        compression.
        """
        # The post_volume array is indexed from the second element to
        # eliminate the duplicated element from the end of the stroke
        # volume array.
        self.volume = np.concatenate((self.stroke_volume,
                                      self.post_volume[1:]))
        """The final volume trace created by concatenating the
        `stroke_volume` and `post_volume` arrays.
        """

        self.computed_pressure = PressureFromVolume(
            self.volume[::5],
            self.stroke_pressure[0]*1E5,
            self.reactive_trace.Tin,
        ).pressure
        """
        A pressure trace computed from the computed volume trace
        to ensure that the volume trace was produced properly.
        """

    def write_output(self):
        """
        Write the output files from the volume and pressure traces.
        The traces are sampled every 5 points to reduce the amount
        of required computational time. This has negligible effect on
        the computed pressure trace. The volume trace is written in
        `csv` format to the file `volume.csv` and the experimental
        pressure trace is written to the file
        `Tc__P0__T0_XXXK_pressure.txt`, where `XXX` represent the
        initial temperature of the experiment. The user should fill in
        the missing values into the file name after simulations are
        completed and the values are known.
        """
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
        """
        Plot the original reactive pressure trace, the pressure trace
        that will be output to the file, and the pressure trace
        computed from the volume trace, to ensure that everything
        matches and the compression time is set properly. Also plot the
        line fitted to the initial part of the full experimental trace
        to help judge the compression time.
        """
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
        """
        When the instance of the class is called, run the analysis
        again.
        """
        self.create_volume_trace()
        self.plot_figure()

        print('{:.4f}'.format(self.stroke_pressure[0]))
        copy('{:.4f}'.format(self.stroke_pressure[0]))

        self.write_output()
