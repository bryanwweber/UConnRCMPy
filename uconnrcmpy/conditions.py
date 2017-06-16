"""Data Processing Module"""

# System imports
from pathlib import Path
from glob import glob
import platform
import sys

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml
import cantera as ct

# Local imports
from pyperclip import copy
from .traces import (TemperatureFromPressure,
                     VolumeFromPressure,
                     PressureFromVolume,
                     )
from .experiments import Experiment, AltExperiment
from .simulations import Simulation


class Condition(object):
    """Class containing all the experiments at a condition.

    Parameters
    ----------
    cti_file : `str` or `pathlib.Path`, optional
        The location of the CTI file to use for Cantera. If it is
        not specified, the default is ``./species.cti``, but if that
        file cannot be found, the user is prompted for a filename.
    plotting : `bool`, optional
        Set to True to enable plotting when experiments are added

    Attributes
    ----------
    reactive_experiments : `dict`
        Dictionary of instances of class `Experiment` containing
        reactive experiments. Indexed by the file name of the
        underlying voltage trace from the experiment.
    nonreactive_experiments : `dict`
        Dictionary of instances of class `Experiment` containing
        nonreactive experiments. Indexed by the file name of the
        underlying voltage trace from the experiment.
    reactive_case : `Experiment`
        Single instance of class `Experiment` that is the closest to
        the mean of all of the experiments.
    nonreactive_case : `Experiment`
        Single instance of class `Experiment` that best matches the
        `~Condition.reactive_case`.
    presout : `numpy.ndarray`
        Pressure array that will be saved to the output file
    volout : `numpy.ndarray`
        Volume array that will be saved to the output file
    nonreactive_sim : `Simulation`
        Instance containing the nonreactive simulation
    reactive_sim : `Simulation`
        Instance containing the reactive simulation
    plotting : `bool`
        Set to True when plotting of experiments is enabled
    all_runs_figure : `matplotlib.figure.Figure`
        Figure showing all the runs at a condition
    nonreactive_figure : `matplotlib.figure.Figure`
        Figure showing a comparison of the nonreactive pressure trace
        with the reactive_case
    pressure_comparison_figure : `matplotlib.figure.Figure`
        Comparison figure of the reactive_case pressure and the
        pressure calculated by the volume trace routine
    simulation_figure : `matplotlib.figure.Figure`
        Comparison of the simulation with the reactive_case pressure
    """

    def __init__(self, cti_file=None, plotting=True):
        self.reactive_experiments = {}
        self.nonreactive_experiments = {}
        self.reactive_case = None
        self.reactive_file = None
        self.nonreactive_case = None
        self.nonreactive_file = None
        self.presout = None
        self.volout = None
        self.nonreactive_sim = None
        self.reactive_sim = None
        self.plotting = plotting
        if self.plotting:
            self.all_runs_figure = None
            self.all_runs_lines = {}
            self.nonreactive_figure = None
            self.pressure_comparison_figure = None
            self.simulation_figure = None

        self.output_attributes = [
            'reactive_file', 'nonreactive_file', 'nonreactive_end_time', 'reactive_end_time',
            'reactive_compression_time', 'nonreactive_offset_points', 'reactive_offset_points'
        ]
        if cti_file is None:
            path_args = {'strict': True} if sys.version_info >= (3, 6) else {}
            try:
                cti_file = str(Path('./species.cti').resolve(**path_args))
            except FileNotFoundError:
                cti_file = str(Path(input('Input the name of the CTI file: ')).resolve())

        self.cti_file = cti_file

        with open(str(cti_file), 'r') as in_file:
            self.cti_source = in_file.read()

        ct.suppress_thermo_warnings(False)
        ct.Solution(self.cti_file)
        ct.suppress_thermo_warnings()

    def __repr__(self):
        return 'Condition(cti_file, plotting={self.plotting!r})'.format(self=self)

    def summary(self):
        summary_str = [('Date-Time     P0 (Torr)  T0 (Kelvin)  Pc (bar)  τ (ms)  τ1 (ms)\n'
                        '---------     ---------  -----------  --------  ------  -------')]
        spacer = len(summary_str[0].split('\n')[0])*'-'
        summary_str.append('Reactive Experiments')
        summary_str.append(spacer)
        fmt = '{date}  {P0:^9d}  {T0:^11d}  {Pc:^8.2f}  {tau_ig:^6.2f}  {tau_first:^7.2f}'
        if not self.reactive_experiments:
            summary_str.append('No reactive experiments in this Condition.\n'
                               'Add some with add_experiment().')
        else:
            for key, expt in self.reactive_experiments.items():
                summary_str.append(fmt.format(
                    date=expt.experiment_parameters['date'],
                    P0=expt.experiment_parameters['pin'],
                    T0=expt.experiment_parameters['Tin'],
                    Pc=expt.pressure_trace.p_EOC,
                    tau_ig=expt.ignition_delay,
                    tau_first=expt.first_stage,
                ))
        summary_str.append(spacer)
        summary_str.append('Non-Reactive Experiments')
        summary_str.append(spacer)
        if not self.nonreactive_experiments:
            summary_str.append('No non-reactive experiments in this Condition.\n'
                               'Add some with add_experiment().')
        else:
            for key, expt in self.nonreactive_experiments.items():
                summary_str.append(fmt.format(
                    date=expt.experiment_parameters['date'],
                    P0=expt.experiment_parameters['pin'],
                    T0=expt.experiment_parameters['Tin'],
                    Pc=expt.pressure_trace.p_EOC,
                    tau_ig=expt.ignition_delay,
                    tau_first=expt.first_stage,
                ))
        return '\n'.join(summary_str)

    def __call__(self):
        print(self.summary())

    @property
    def reactive_file(self):
        """`pathlib.Path`: File name of the case that is the closest to the mean of all
        the experiments. Can be set with a `str` or `~pathlib.Path`.
        """
        return self._reactive_file

    @reactive_file.setter
    def reactive_file(self, value):
        if value is not None:
            self._reactive_file = Path(value).resolve()
        else:
            self._reactive_file = value

    @property
    def nonreactive_file(self):
        """`pathlib.Path`: File name of the non-reactive case that best matches
        the pressure trace of the `~Condition.reactive_case`. Can be set with a
        `str` or `~pathlib.Path`.
        """
        return self._nonreactive_file

    @nonreactive_file.setter
    def nonreactive_file(self, value):
        if value is not None:
            self._nonreactive_file = Path(value).resolve()
        else:
            self._nonreactive_file = value

    @property
    def reactive_end_time(self):
        """`float`: The end time for the output of the reactive case to a file, relative to
        the end of compression in milliseconds
        """
        if self.reactive_case is not None:
            return self.reactive_case.output_end_time
        else:
            raise AttributeError('The reactive_case has not been set for this Condition yet.')

    @reactive_end_time.setter
    def reactive_end_time(self, value):
        if self.reactive_case is not None:
            max_time = self.reactive_case.pressure_trace.zeroed_time[-1]*1000
            if value > max_time:
                raise ValueError('The reactive end time cannot be after the end of the pressure '
                                 'trace. The maximum time is: {} ms.'.format(max_time))

            self.reactive_case.output_end_time = value
        else:
            raise AttributeError('The reactive_case has not been set for this Condition yet.')

    @property
    def nonreactive_end_time(self):
        """`float`: The end time for the output of the nonreactive case to a file, relative to
        the end of compression in milliseconds
        """
        if self.nonreactive_case is not None:
            return self.nonreactive_case.output_end_time
        else:
            raise AttributeError('The nonreactive_case has not been set for this Condition yet.')

    @nonreactive_end_time.setter
    def nonreactive_end_time(self, value):
        if self.nonreactive_case is not None:
            max_time = self.nonreactive_case.pressure_trace.zeroed_time[-1]*1000
            if value > max_time:
                raise ValueError('The nonreactive end time cannot be after the end of the pressure'
                                 ' trace. The maximum time is: {} ms.'.format(max_time))

            self.nonreactive_case.output_end_time = value
        else:
            raise AttributeError('The nonreactive_case has not been set for this Condition yet.')

    @property
    def reactive_compression_time(self):
        """`float`: The compression time for the reactive case, in milliseconds
        """
        if self.reactive_case is not None:
            return self.reactive_case.compression_time
        else:
            raise AttributeError('The reactive_case has not been set for this Condition yet.')

    @reactive_compression_time.setter
    def reactive_compression_time(self, value):
        if self.reactive_case is not None:
            max_time = np.abs(self.reactive_case.pressure_trace.zeroed_time[0]*1000)
            if value > max_time:
                raise ValueError('The compression time cannot be longer than the time before'
                                 'compression. The maximum time is: {} ms'.format(max_time))

            self.reactive_case.compression_time = value
        else:
            raise AttributeError('The reactive_case has not been set for this Condition yet.')

    @property
    def reactive_offset_points(self):
        """`int`: The number of points to offset the end of compression, to better
        match the non-reactive and reactive pressure traces together"""
        if self.reactive_case is not None:
            return self.reactive_case.offset_points
        else:
            raise AttributeError('The reactive_case has not been set for this Condition yet.')

    @reactive_offset_points.setter
    def reactive_offset_points(self, value):
        if self.reactive_case is not None:
            self.reactive_case.offset_points = value
        else:
            raise AttributeError('The reactive_case has not been set for this Condition yet.')

    @property
    def nonreactive_offset_points(self):
        """`int`: The number of points to offset the end of compression, to better
        match the non-reactive and reactive pressure traces together"""
        if self.nonreactive_case is not None:
            return self.nonreactive_case.offset_points
        else:
            raise AttributeError('The nonreactive_case has not been set for this Condition yet.')

    @nonreactive_offset_points.setter
    def nonreactive_offset_points(self, value):
        if self.nonreactive_case is not None:
            self.nonreactive_case.offset_points = value
        else:
            raise AttributeError('The nonreactive_case has not been set for this Condition yet.')

    def add_reactive_case(self):
        """Add the `~Condition.reactive_case` to the `Condition` if
        it hasn't been set yet. If the `~Condition.reactive_case`
        exists but the name of the `~Condition.reactive_file` is
        different than the name of the existing case, the case is
        replaced. If neither of these is true, nothing is done.
        """
        if (self.reactive_case is None or
                self.reactive_case.file_path.name != self.reactive_file.name):
            if self.reactive_file in self.reactive_experiments:
                self.reactive_case = self.reactive_experiments[self.reactive_file]
            else:
                self.reactive_case = Experiment(self.reactive_file, cti_file=self.cti_file)
                self.reactive_experiments[self.reactive_file] = self.reactive_case

    def add_nonreactive_case(self):
        """Add the `~Condition.nonreactive_case` to the `Condition` if
        it hasn't been set yet. If the `~Condition.nonreactive_case`
        exists but the name of the `~Condition.nonreactive_file` is
        different than the name of the existing case, the case is
        replaced. If neither of these is true, nothing is done.
        """
        if (self.nonreactive_case is None or
                self.nonreactive_case.file_path.name != self.nonreactive_file.name):
            if self.nonreactive_file in self.nonreactive_experiments:
                self.nonreactive_case = self.nonreactive_experiments[self.nonreactive_file]
            else:
                self.nonreactive_case = Experiment(self.nonreactive_file, cti_file=self.cti_file)
                self.nonreactive_experiments[self.nonreactive_file] = self.nonreactive_case

    def add_experiment(self, file_name=None, **kwargs):
        """Add an experiment to the Condition.

        Parameters
        ----------
        file_name : `pathlib.Path` or `str` or `None`
            Filename of the file with the voltage trace of the
            experiment to be added.
        """
        exp = Experiment(file_name, cti_file=self.cti_file, **kwargs)
        if exp.pressure_trace.is_reactive:
            self.reactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_reactive_figures(exp)
        else:
            self.nonreactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_nonreactive_figure(exp)

    def change_EOC_time(self, experiment, time, is_reactive=True):
        """Change the EOC time manually for an experiment

        Parameters
        ----------
        experiment : `Experiment` or `str` or `pathlib.Path`
            The experiment to be modified, either as a filename
            or as the instance of an `Experiment` class.
        time : `float`
            The value for the EOC time
        is_reactive : `boolean`
            The experiment is reactive or not
        """
        if isinstance(experiment, (str, Path)):
            if str(experiment) in self.reactive_experiments:
                experiment = self.reactive_experiments[str(experiment)]
            elif str(experiment) in self.nonreactive_experiments:
                experiment = self.nonreactive_experiments[str(experiment)]
            else:
                raise ValueError('{} could not be found in the Condition. '
                                 'Did you add it with add_experiment?'.format(str(experiment)))

        experiment.change_EOC_time(time, is_reactive)
        if self.plotting:
            self.all_runs_lines[experiment.file_path.name].remove()
            self.all_runs_lines[experiment.file_path.name], = self.all_runs_axis.plot(
                experiment.pressure_trace.zeroed_time*1000.0,
                experiment.pressure_trace.pressure,
                label=experiment.experiment_parameters['date'],
            )
            self.all_runs_axis.legend(loc='best')

    def change_filter_freq(self, experiment, value):
        """Change the cutoff frequency of the filter for an experiment

        Parameters
        ----------
        experiment : `Experiment` or `str` or `pathlib.Path`
            The experiment to be modified, either as a filename
            or as the instance of an `Experiment` class.
        value : `float`
            The value for the cutoff frequency
        """
        if isinstance(experiment, (str, Path)):
            if str(experiment) in self.reactive_experiments:
                experiment = self.reactive_experiments[str(experiment)]
            elif str(experiment) in self.nonreactive_experiments:
                experiment = self.nonreactive_experiments[str(experiment)]
            else:
                raise ValueError('{} could not be found in the Condition. '
                                 'Did you add it with add_experiment?'.format(str(experiment)))

        experiment.change_filter_freq(value)
        if self.plotting:
            self.all_runs_lines[experiment.file_path.name].remove()
            self.all_runs_lines[experiment.file_path.name], = self.all_runs_axis.plot(
                experiment.pressure_trace.zeroed_time*1000.0,
                experiment.pressure_trace.pressure,
                label=experiment.experiment_parameters['date'],
            )
            self.all_runs_axis.legend(loc='best')

    def plot_reactive_figures(self, exp):
        """Plot the reactive pressure trace on figures.

        Plots the reactive pressure trace on two figures. The first is
        a comparison of all of the processed reactive pressure traces.
        The second plot is a comparison of the pressure trace and the
        derivative for a single experiment.

        Parameters
        ----------
        exp : `Experiment`
            Experiment to plot
        """
        # Plot the filtered pressure and overlay future runs
        if self.all_runs_figure is None:
            self.all_runs_figure = plt.figure('Reactive Pressure Trace Comparison')
            self.all_runs_axis = self.all_runs_figure.add_subplot(1, 1, 1)
            self.all_runs_axis.set_ylabel('Pressure [bar]')
            self.all_runs_axis.set_xlabel('Time [ms]')
            if platform.system() == 'Windows':
                m = plt.get_current_fig_manager()
                m.window.showMaximized()

        self.all_runs_lines[exp.file_path.name], = self.all_runs_axis.plot(
            exp.pressure_trace.zeroed_time*1000.0,
            exp.pressure_trace.pressure,
            label=exp.experiment_parameters['date'],
        )
        self.all_runs_axis.legend(loc='best')

        exp.plot_pressure_trace()

    def load_yaml(self):
        """
        Load the yaml file called ``volume-trace.yaml`` containing the
        information needed to construct the volume trace. This method
        reads the YAML file and sets the appropriate attributes of the
        `Condition` instance. All the following data are required unless
        otherwise noted. The format of the yaml file is::

            variable: value

        * ``nonreactive_file``: File name of the non-reactive pressure trace.
            Type: String
        * ``reactive_file``: File name of the reactive pressure trace.
            Type: String
        * ``reactive_compression_time``: Length of time of the compression stroke, in ms.
            Type: Float
        * ``nonreactive_end_time``: End time used for the produced volume trace, in ms.
            Type: Float
        * ``reactive_end_time``: End time of the output reactive pressure trace, in ms.
            Type: Float
        * ``reactive_offset_points``: Offset in number of points from EOC for the
            reactive case. Optional, defaults to zero.
            Type: Integer
        * ``nonreactive_offset_points``: Offset in number of points from EOC for the
            non-reactive case. Optional, defaults to zero.
            Type: Integer
        """
        with open('volume-trace.yaml') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        for attribute in ['nonreactive_file', 'reactive_file']:
            if attribute in yaml_data:
                setattr(self, attribute, yaml_data[attribute])

        self.add_reactive_case()
        self.add_nonreactive_case()

        for attribute in self.output_attributes:
            if attribute in yaml_data:
                setattr(self, attribute, yaml_data[attribute])

    def load_old_yaml(self):
        """
        Load the yaml file called ``volume-trace.yaml`` containing the
        information needed to construct the volume trace. All the
        following data are required unless otherwise noted. The format
        of the yaml file is::

            variable: value

        * ``nonrfile``: File name of the non-reactive pressure trace.
                        Type: String
        * ``reacfile``: File name of the reactive pressure trace.
                        Type: String
        * ``comptime``: Length of time of the compression stroke, in ms.
                        Type: Float
        * ``nonrend``: End time used for the produced volume trace, in ms.
                        Type: Float
        * ``reacend``: End time of the output reactive pressure trace, in ms.
                        Type: Float
        * ``reacoffs``: Offset in number of points from EOC for the
                        reactive case. Optional, defaults to zero.
                        Type: Integer
        * ``nonroffs``: Offset in number of points from EOC for the
                        non-reactive case. Optional, defaults to zero.
                        Type: Integer
        """
        with open('volume-trace.yaml') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        attributes = [('reacfile', 'reactive_file'),
                      ('nonrfile', 'nonreactive_file'),
                      ('nonrend', 'nonreactive_end_time'),
                      ('reacend', 'reactive_end_time'),
                      ('comptime', 'reactive_compression_time'),
                      ('nonroffs', 'nonreactive_offset_points'),
                      ('reacoffs', 'reactive_offset_points')]

        for yaml_key, attribute in attributes:
            if yaml_key in yaml_data:
                setattr(self, attribute, yaml_data[yaml_key])

    def plot_nonreactive_figure(self, exp):
        """Plot the nonreactive pressure traces on a figure.

        Plots the nonreactive pressure traces in comparison with the
        reactive pressure trace. Adds additional traces each time.

        Parameters
        ----------
        exp : `Experiment`
            Experiment to plot
        """
        if self.nonreactive_figure is None:
            self.nonreactive_figure = plt.figure('Non-Reactive Pressure Trace Comparison')
            self.nonreactive_axis = self.nonreactive_figure.add_subplot(1, 1, 1)
            self.nonreactive_axis.set_ylabel('Pressure [bar]')
            self.nonreactive_axis.set_xlabel('Time [ms]')
            if platform.system() == 'Windows':
                m = plt.get_current_fig_manager()
                m.window.showMaximized()

            if self.reactive_file is None:
                self.reactive_file = input('Reactive filename: ')

            self.add_reactive_case()

            self.nonreactive_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time*1000.0,
                self.reactive_case.pressure_trace.pressure,
                label=self.reactive_case.experiment_parameters['date'],
            )

            linear_fit = self.reactive_case.pressure_trace.pressure_fit()
            reactive_line = np.polyval(linear_fit, self.reactive_case.pressure_trace.time)
            self.nonreactive_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time*1000.0,
                reactive_line,
                label='Linear Fit to Initial Pressure',
            )

        self.nonreactive_axis.plot(
            exp.pressure_trace.zeroed_time*1000.0,
            exp.pressure_trace.pressure,
            label=exp.experiment_parameters['date'],
        )
        self.nonreactive_axis.legend(loc='best')

    def create_volume_trace(self):
        """
        Create the volume trace based on the information in the loaded
        yaml file.
        """
        if self.reactive_file is None:
            self.reactive_file = input('Reactive filename: ')

        self.add_reactive_case()

        if self.nonreactive_file is None:
            self.nonreactive_file = input('Non-Reactive filename: ')

        self.add_nonreactive_case()

        # Set the offsets to zero by default for the initial run through
        # The user should modify these after the first run if necessary
        for attribute in ['reactive_offset_points', 'nonreactive_offset_points']:
            if not hasattr(self, attribute) or getattr(self, attribute) is None:
                setattr(self, attribute, 0)

        for attribute in self.output_attributes:
            if not hasattr(self, attribute) or getattr(self, attribute) is None:
                temp_val = input('Specify a value for the {}: '.format(attribute))
                setattr(self, attribute, float(temp_val))

        nonreactive_end_idx = int(
            self.nonreactive_case.pressure_trace.EOC_idx + self.nonreactive_offset_points +
            self.nonreactive_end_time/1000.0*self.nonreactive_case.pressure_trace.frequency
        )

        reactive_end_idx = int(
            self.reactive_case.pressure_trace.EOC_idx + self.reactive_offset_points +
            self.reactive_end_time/1000.0*self.reactive_case.pressure_trace.frequency
        )

        reactive_start_point = int(
            self.reactive_case.pressure_trace.EOC_idx + self.reactive_offset_points -
            self.reactive_compression_time/1000.0*self.reactive_case.pressure_trace.frequency
        )

        stroke_pressure = self.reactive_case.pressure_trace.pressure[
            (reactive_start_point):(self.reactive_case.pressure_trace.EOC_idx +
                                    1 + self.reactive_offset_points)
        ]

        post_pressure = self.nonreactive_case.pressure_trace.pressure[
            (self.nonreactive_case.pressure_trace.EOC_idx + self.nonreactive_offset_points):(
                nonreactive_end_idx
            )
        ]

        print_pressure = self.reactive_case.pressure_trace.pressure[
            reactive_start_point:reactive_end_idx
        ]

        n_print_pts = len(print_pressure)
        time = np.arange(-self.reactive_compression_time/1000.0, self.nonreactive_end_time/1000.0,
                         1/self.reactive_case.pressure_trace.frequency)
        stroke_volume = VolumeFromPressure(
            stroke_pressure,
            1.0,
            self.reactive_case.experiment_parameters['Tin'],
            chem_file=str(self.cti_file),
        ).volume
        stroke_temperature = TemperatureFromPressure(
            stroke_pressure,
            self.reactive_case.experiment_parameters['Tin'],
            chem_file=str(self.cti_file),
        ).temperature

        post_volume = VolumeFromPressure(
            post_pressure,
            stroke_volume[-1],
            stroke_temperature[-1],
            chem_file=str(self.cti_file),
        ).volume

        # The post_volume array is indexed from the second element to
        # eliminate the duplicated element from the end of the stroke
        # volume array.
        volume = np.concatenate((stroke_volume, post_volume[1:]))

        computed_pressure = PressureFromVolume(
            volume[::5],
            stroke_pressure[0]*1E5,
            self.reactive_case.experiment_parameters['Tin'],
            chem_file=str(self.cti_file),
        ).pressure

        copy('{:.4f}'.format(stroke_pressure[0]))

        self.volout = np.vstack(
            (time[::5] + self.reactive_compression_time/1000, volume[::5])
        ).transpose()

        self.presout = np.vstack(
            (time[:n_print_pts:5] + self.reactive_compression_time/1000,
             print_pressure[::5])
        ).transpose()

        self.write_output(
            self.volout,
            self.presout,
            self.reactive_case.experiment_parameters['Tin'],
        )

        self.write_yaml()

        if self.plotting:
            if self.pressure_comparison_figure is None:
                self.pressure_comparison_figure = plt.figure('Pressure Trace Comparison')
                self.pressure_comparison_axis = self.pressure_comparison_figure.add_subplot(1, 1, 1)
                if platform.system() == 'Windows':
                    m = plt.get_current_fig_manager()
                    m.window.showMaximized()

            self.pressure_comparison_axis.cla()
            self.pressure_comparison_axis.set_ylabel("Pressure [bar]")
            self.pressure_comparison_axis.set_xlabel("Time [ms]")
            plot_time = (self.reactive_case.pressure_trace.zeroed_time -
                         self.reactive_offset_points/self.reactive_case.pressure_trace.frequency)
            self.pressure_comparison_axis.plot(
                plot_time*1000.0,
                self.reactive_case.pressure_trace.pressure,
                label="Reactive Pressure",
            )
            self.pressure_comparison_axis.plot(time[:n_print_pts:5]*1000.0, print_pressure[::5],
                                               label="Output Pressure")
            self.pressure_comparison_axis.plot(time[::5]*1000.0, computed_pressure,
                                               label="Computed Pressure")
            linear_fit = self.reactive_case.pressure_trace.pressure_fit(
                comptime=self.reactive_compression_time/1000,
            )
            reactive_line = np.polyval(linear_fit, self.reactive_case.pressure_trace.time)
            self.pressure_comparison_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time*1000.0,
                reactive_line,
                label="Linear Fit to Initial Pressure",
            )
            self.pressure_comparison_axis.legend(loc='best')

    def write_yaml(self):
        """Write the volume-trace.yaml output file for storage of parameters.
        The YAML file format is detailed in `~Condition.load_yaml`.
        """

        yaml_data = {}
        for attribute in self.output_attributes:
            val = getattr(self, attribute)
            if isinstance(val, Path):
                yaml_data[attribute] = val.name
            else:
                yaml_data[attribute] = val

        with open('volume-trace.yaml', 'w') as yaml_file:
            yaml.dump(yaml_data, yaml_file, default_flow_style=False)

    def write_output(self, volout, presout, Tin):
        """
        Write the output files from the volume and pressure traces.
        The traces are sampled every 5 points to reduce the amount
        of required computational time. This has negligible effect on
        the computed pressure trace. The volume trace is written in
        ``csv`` format to the file ``volume.csv`` and the experimental
        pressure trace is written to the file
        ``Tc__P0__T0_XXXK_pressure.txt``, where ``XXX`` represents the
        initial temperature of the experiment. The user should fill the
        missing values into the file name after simulations are
        completed and the values are known.
        """
        np.savetxt('volume.csv', volout, delimiter=',')
        np.savetxt('Tc__P0__T0_{}K_pressure.txt'.format(Tin), presout, delimiter='\t')

    def run_simulation(self, run_reactive=False, run_nonreactive=True,
                       end_temp=2500.0, end_time=0.2):
        """Run the simulations for this condition.

        Parameters
        ----------
        run_reactive : `bool`, optional
            True to run the reactive case. False by default.
        run_nonreactive : `bool`, optional
            True to run the nonreactive case. True by default.
        end_temp : `float`, optional
            Temperature at which the simulation is ended
        end_time : `float`, optional
            Time at which the simulation is ended.
        """
        def process_choice(sim_type):
            choice = input('Are you sure you want to overwrite the {sim_type} simulation? '
                           'Input y or n: '.format(sim_type=sim_type))
            if choice.startswith('n'):
                return False
            elif choice.startswith('y'):
                return True
            else:
                raise IOError('Invalid input')

        if run_nonreactive:
            if self.nonreactive_sim is None:
                self.nonreactive_sim = Simulation(
                    initial_temperature=self.reactive_case.experiment_parameters['Tin'],
                    initial_pressure=self.presout[0, 1]*1E5,
                    volume=self.volout,
                    is_reactive=False,
                    end_temp=end_temp,
                    end_time=end_time,
                    chem_file=str(self.cti_file),
                )
            else:
                if process_choice('nonreactive'):
                    self.nonreactive_sim = Simulation(
                        initial_temperature=self.reactive_case.experiment_parameters['Tin'],
                        initial_pressure=self.presout[0, 1]*1E5,
                        volume=self.volout,
                        is_reactive=False,
                        end_temp=end_temp,
                        end_time=end_time,
                        chem_file=str(self.cti_file),
                    )
                else:
                    print('Nothing was done')

        if run_reactive:
            if self.reactive_sim is None:
                self.reactive_sim = Simulation(
                    initial_temperature=self.reactive_case.experiment_parameters['Tin'],
                    initial_pressure=self.presout[0, 1]*1E5,
                    volume=self.volout,
                    is_reactive=True,
                    end_temp=end_temp,
                    end_time=end_time,
                    chem_file=str(self.cti_file),
                )
            else:
                if process_choice('reactive'):
                    self.reactive_sim = Simulation(
                        initial_temperature=self.reactive_case.experiment_parameters['Tin'],
                        initial_pressure=self.presout[0, 1]*1E5,
                        volume=self.volout,
                        is_reactive=True,
                        end_temp=end_temp,
                        end_time=end_time,
                        chem_file=str(self.cti_file),
                    )
                else:
                    print('Nothing was done')

    def compare_to_sim(self, run_reactive=False, run_nonreactive=True):
        """Compare the experiments to the simulations.

        Run the simulations for this condition, and if plotting is on,
        generate the comparison plot.

        Parameters
        ----------
        run_reactive : `bool`, optional
            True to run the reactive comparison. False by default.
        run_nonreactive : `bool`, optional
            True to run the nonreactive comparison. True by default.
        """
        if self.presout is None:
            # Load the experimental pressure trace. Try the glob function first
            # and if it fails, ask the user for help.
            flist = glob('*pressure.txt')
            if not len(flist) == 1:
                flist = [input('Input the experimental pressure trace file name: ')]
            self.presout = np.genfromtxt(flist[0])

        self.run_simulation(run_reactive, run_nonreactive)

        # Plot the pressure traces together
        compression_time = self.reactive_compression_time

        if self.plotting:
            if self.simulation_figure is None:
                self.simulation_figure = plt.figure('Simulation Comparison')
                self.simulation_axis = self.simulation_figure.add_subplot(1, 1, 1)
                if platform.system() == 'Windows':
                    m = plt.get_current_fig_manager()
                    m.window.showMaximized()

            self.simulation_axis.cla()
            self.simulation_axis.set_xlabel('Time [ms]')
            self.simulation_axis.set_ylabel('Pressure [bar]')

            self.simulation_axis.plot(
                self.presout[:, 0]*1000 - compression_time,
                self.presout[:, 1],
                label="Experimental Pressure"
            )

            if self.nonreactive_sim is not None:
                self.simulation_axis.plot(
                    self.nonreactive_sim.time*1000 - compression_time,
                    self.nonreactive_sim.pressure,
                    label="Simulated Non-Reactive Pressure"
                )

            if self.reactive_sim is not None:
                self.simulation_axis.plot(
                    self.reactive_sim.time*1000 - compression_time,
                    self.reactive_sim.pressure,
                    label="Simulated Reactive Pressure"
                )
                der = self.reactive_sim.derivative
                self.simulation_axis.plot(
                    self.reactive_sim.time*1000 - compression_time,
                    der/np.amax(der)*np.amax(self.reactive_sim.pressure),
                    label="Scaled Derivative of Simulated Pressure"
                )

            self.simulation_axis.legend(loc='best')

        print_str = ''
        copy_str = ''

        if self.nonreactive_sim is not None:
            T_EOC = np.amax(self.nonreactive_sim.temperature)
            print_str += '{:.0f}'.format(T_EOC)
            copy_str += '{}'.format(T_EOC)

        if self.reactive_sim is not None:
            if self.nonreactive_sim is not None:
                print_str += '\t'
                copy_str += '\t\t\t\t'
            ignition_idx = np.argmax(self.reactive_sim.derivative)
            ignition_delay = self.reactive_sim.time[ignition_idx]*1000 - compression_time
            print_str += '{:.6f}'.format(ignition_delay)
            copy_str += '{}'.format(ignition_delay)

        if self.nonreactive_sim is not None and self.reactive_sim is not None:
            # Use np.where to find the index of EOC for reactive simulation
            # Use the given index to find T_EOC_reactive
            # Compare to T_EOC (nonreactive case) using np.isclose
            EOC_reactive_idx = np.where(self.reactive_sim.time >= compression_time/1000)
            T_EOC_reactive = self.reactive_sim.temperature[EOC_reactive_idx[0][0]]
            if not np.isclose(T_EOC, T_EOC_reactive, rtol=0.001):
                print('WARNING: Temperatures at Reactive EOC and non-reactive EOC did not match')
                print('T_EOC_nonreactive = {} K'.format(T_EOC))
                print('T_EOC_reactive = {} K'.format(T_EOC_reactive))
        print(print_str)
        copy(copy_str)


class AltCondition(Condition):
    """Class containing all of the alternate experiments at a condition
    """

    def add_experiment(self, file_name=None):
        """Add an experiment to the Condition.

        Parameters
        ----------
        file_name : `str` or `None`
            Filename of the file with the voltage trace of the
            experiment to be added.
        """
        exp = AltExperiment(file_name, cti_file=self.cti_file)
        if exp.pressure_trace.is_reactive:
            self.reactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_reactive_figures(exp)
        else:
            self.nonreactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_nonreactive_figure(exp)

    def add_reactive_case(self):
        """Add the `~AltCondition.reactive_case` to the `AltCondition` if
        it hasn't been set yet. If the `~AltCondition.reactive_case`
        exists but the name of the `~AltCondition.reactive_file` is
        different than the name of the existing case, the case is
        replaced. If neither of these is true, nothing is done.
        """
        if (self.reactive_case is None or
                self.reactive_case.file_path.name != self.reactive_file.name):
            if self.reactive_file in self.reactive_experiments:
                self.reactive_case = self.reactive_experiments[self.reactive_file]
            else:
                self.reactive_case = AltExperiment(self.reactive_file, cti_file=self.cti_file)
                self.reactive_experiments[self.reactive_file] = self.reactive_case

    def add_nonreactive_case(self):
        """Add the `~AltCondition.nonreactive_case` to the `AltCondition` if
        it hasn't been set yet. If the `~AltCondition.nonreactive_case`
        exists but the name of the `~AltCondition.nonreactive_file` is
        different than the name of the existing case, the case is
        replaced. If neither of these is true, nothing is done.
        """
        if (self.nonreactive_case is None or
                self.nonreactive_case.file_path.name != self.nonreactive_file.name):
            if self.nonreactive_file in self.nonreactive_experiments:
                self.nonreactive_case = self.nonreactive_experiments[self.nonreactive_file]
            else:
                self.nonreactive_case = AltExperiment(self.nonreactive_file, cti_file=self.cti_file)
                self.nonreactive_experiments[self.nonreactive_file] = self.nonreactive_case

    def __repr__(self):
        return 'AltCondition(plotting={!r})'.format(self.plotting)


def process_folder(cti_file, path='.', plot=False):
    """Process a folder of experimental files.

    Process a folder containing files with reactive experiments to
    calculate the ignition delays and copy a table with the results
    to the clipboard.

    Parameters
    ----------
    cti_file : `str` or `pathlib.Path`
        File containing the CTI for Cantera
    path : `str`, optional
        Path to folder to be analyzed. Defaults to the current folder.
    plot : `bool`, optional
        True to enable plotting. False by default.
    """
    p = Path(path)
    result = []

    cti_file = Path(cti_file).resolve()

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    for f in p.glob('[0-3]*.txt'):
        print(f)
        case = Experiment(f.resolve(), cti_file=cti_file, copy=False)
        result.append('\t'.join(map(str, [
            case.experiment_parameters['time_of_day'], case.experiment_parameters['pin'],
            case.experiment_parameters['Tin'], case.pressure_trace.p_EOC, case.ignition_delay,
            case.first_stage, case.T_EOC, case.experiment_parameters['spacers'],
            case.experiment_parameters['shims'], case.voltage_trace.filter_frequency])))
        if plot:
            ax.plot(case.pressure_trace.zeroed_time, case.pressure_trace.pressure,
                    label=case.experiment_parameters['date'])

    copy('\n'.join(sorted(result)))
    print('Finished')


def process_alt_folder(cti_file, path='.', plot=False):
    """Process a folder of alternative experimental files.

    Process a folder containing files with reactive experiments to
    calculate the ignition delays and copy a table with the results
    to the clipboard.

    Parameters
    ----------
    cti_file : `str` or `pathlib.Path`
        File containing the CTI for Cantera
    path : `str`, optional
        Path to folder to be analyzed. Defaults to the current folder.
    plot : `bool`, optional
        True to enable plotting. False by default.
    """
    p = Path(path)
    result = []

    cti_file = Path(cti_file).resolve()

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    for f in p.glob('Fuel_MF*.txt'):
        if 'NR' in f.name:
            continue
        print(f)
        case = AltExperiment(f.resolve(), cti_file=cti_file, copy=False)
        result.append('\t'.join(map(str, [
            case.experiment_parameters['date_year'], case.experiment_parameters['time_of_day'],
            case.experiment_parameters['pin'], case.experiment_parameters['Tin'],
            case.pressure_trace.p_EOC, case.ignition_delay, case.first_stage, case.T_EOC,
            case.experiment_parameters['spacers'], case.experiment_parameters['shims'],
            case.pressure_trace.filter_frequency, f.name])))
        if plot:
            ax.plot(case.pressure_trace.zeroed_time, case.pressure_trace.pressure,
                    label=case.experiment_parameters['date'])

    copy('\n'.join(sorted(result)))
    print('Finished')


if __name__ == '__main__':
    process_folder()
