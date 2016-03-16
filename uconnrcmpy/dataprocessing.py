"""Data Processing Module"""

# System imports
from pathlib import Path
from glob import glob
import platform
from datetime import datetime

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml
import cantera as ct
from cansen.profiles import VolumeProfile

# Local imports
from .utilities import copy
from .traces import (VoltageTrace,
                     ExperimentalPressureTrace,
                     AltExperimentalPressureTrace,
                     TemperatureFromPressure,
                     VolumeFromPressure,
                     PressureFromVolume,
                     )
from .constants import cantera_version


class Condition(object):
    """Class containing all the experiments at a condition.

    Parameters
    ----------
    plotting : `bool`, optional
        Set to True to enable plotting when experiments are added

    Attributes
    ----------
    reactive_experiments : `dict`
        Dictionary of instances of class `Experiment` containing
        reactive experiments. Indexed by the DD-Mon-YYYY-Time of the
        experiment.
    nonreactive_experiments : `dict`
        Dictionary of instances of class `Experiment` containing
        nonreactive experiments. Indexed by the DD-Mon-YYYY-Time of the
        experiment.
    reactive_case : `Experiment`
        Single instance of class `Experiment` that is the closest to
        the mean of all of the experiments.
    nonreactive_case : `Experiment`
        Single instance of class `Experiment` that best matches the
        ``reactive_case``.
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

    def __init__(self, plotting=True):
        self.reactive_experiments = {}
        self.nonreactive_experiments = {}
        self.reactive_case = None
        self.nonreactive_case = None
        self.presout = None
        self.volout = None
        self.nonreactive_sim = None
        self.reactive_sim = None
        self.plotting = plotting
        if self.plotting:
            self.all_runs_figure = None
            self.nonreactive_figure = None
            self.pressure_comparison_figure = None
            self.simulation_figure = None

    def add_experiment(self, file_name=None):
        """Add an experiment to the Condition.

        Parameters
        ----------
        file_name : `str` or `None`
            Filename of the file with the voltage trace of the
            experiment to be added.
        """
        exp = Experiment(file_name)
        if exp.pressure_trace.is_reactive:
            self.reactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_reactive_figures(exp)
        else:
            self.nonreactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_nonreactive_figure(exp)

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
        # Plot the smoothed pressure and overlay future runs
        if self.all_runs_figure is None:
            self.all_runs_figure = plt.figure('Reactive Pressure Trace Comparison')
            self.all_runs_axis = self.all_runs_figure.add_subplot(1, 1, 1)
            if platform.system() == 'Windows':
                m = plt.get_current_fig_manager()
                m.window.showMaximized()

        self.all_runs_axis.plot(
            exp.pressure_trace.zeroed_time,
            exp.pressure_trace.pressure,
            label=exp.experiment_parameters['date'],
        )

        exp.plot_pressure_trace()

    def load_yaml(self):
        """
        Load the yaml file called ``volume-trace.yml`` containing the
        information needed to construct the volume trace. All the
        following data are required unless otherwise noted. The format
        of the yaml file is::

            variable: value

        * ``nonrfile``: File name of the non-reactive pressure trace.
                        Type: String
        * ``reacfile``: File name of the reactive pressure trace.
                        Type: String
        * ``comptime``: Length of time of the compression stroke.
                        Type: Integer
        * ``nonrend``: End time used for the produced volume trace.
                        Type: Integer
        * ``reacend``: End time of the output reactive pressure trace.
                        Type: Integer
        * ``reacoffs``: Offset in number of points from EOC for the
                        reactive case. Optional, defaults to zero.
                        Type: Integer
        * ``nonroffs``: Offset in number of points from EOC for the
                        non-reactive case. Optional, defaults to zero.
                        Type: Integer
        """
        with open('volume-trace.yaml') as yaml_file:
            return yaml.load(yaml_file)

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
            if platform.system() == 'Windows':
                m = plt.get_current_fig_manager()
                m.window.showMaximized()

            reactive_file = self.load_yaml()['reacfile']
            self.reactive_case = self.reactive_experiments[reactive_file]
            self.nonreactive_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time,
                self.reactive_case.pressure_trace.pressure,
                label=self.reactive_case.experiment_parameters['date'],
            )

        self.nonreactive_axis.plot(
            exp.pressure_trace.zeroed_time,
            exp.pressure_trace.pressure,
            label=exp.experiment_parameters['date'],
        )

    def create_volume_trace(self):
        """
        Create the volume trace based on the information in the loaded
        yaml file.
        """
        yaml_data = self.load_yaml()
        if self.reactive_case is None:
            try:
                self.reactive_case = self.reactive_experiments[yaml_data['reacfile']]
            except KeyError:
                plotting_old = self.plotting
                self.plotting = False
                self.add_experiment(Path(yaml_data['reacfile']))
                self.reactive_case = self.reactive_experiments[yaml_data['reacfile']]
                self.plotting = plotting_old

        if self.nonreactive_case is None:
            try:
                self.nonreactive_case = self.nonreactive_experiments[yaml_data['nonrfile']]
            except KeyError:
                plotting_old = self.plotting
                self.plotting = False
                self.add_experiment(Path(yaml_data['nonrfile']))
                self.nonreactive_case = self.nonreactive_experiments[yaml_data['nonrfile']]
                self.plotting = plotting_old

        linear_fit = self.reactive_case.pressure_trace.pressure_fit(
            comptime=yaml_data['comptime']/1000,
        )
        reactive_line = np.polyval(linear_fit, self.reactive_case.pressure_trace.time)

        nonreactive_end_idx = (
            self.nonreactive_case.pressure_trace.EOC_idx +
            yaml_data['nonrend']/1000.0*self.nonreactive_case.pressure_trace.frequency
        )

        reactive_end_idx = (
            self.reactive_case.pressure_trace.EOC_idx +
            yaml_data['reacend']/1000.0*self.reactive_case.pressure_trace.frequency
        )

        reactive_start_point = (
            self.reactive_case.pressure_trace.EOC_idx -
            yaml_data['comptime']/1000.0*self.reactive_case.pressure_trace.frequency
        )

        stroke_pressure = self.reactive_case.pressure_trace.pressure[
            (reactive_start_point):(self.reactive_case.pressure_trace.EOC_idx +
                                    1 + yaml_data.get('reacoffs', 0))
        ]

        post_pressure = self.nonreactive_case.pressure_trace.pressure[
            (self.nonreactive_case.pressure_trace.EOC_idx + yaml_data.get('nonroffs', 0)):(
                nonreactive_end_idx + yaml_data.get('nonroffs', 0) - yaml_data.get('reacoffs', 0)
            )
        ]

        print_pressure = self.reactive_case.pressure_trace.pressure[
            (reactive_start_point):(reactive_end_idx)
        ]

        n_print_pts = len(print_pressure)
        time = np.arange(-yaml_data['comptime']/1000.0, yaml_data['nonrend']/1000.0,
                         1/self.reactive_case.pressure_trace.frequency)
        stroke_volume = VolumeFromPressure(stroke_pressure, 1.0,
                                           self.reactive_case.experiment_parameters['Tin']).volume
        stroke_temperature = TemperatureFromPressure(
            stroke_pressure,
            self.reactive_case.experiment_parameters['Tin'],
        ).temperature

        post_volume = VolumeFromPressure(
            post_pressure,
            stroke_volume[-1],
            stroke_temperature[-1],
        ).volume

        # The post_volume array is indexed from the second element to
        # eliminate the duplicated element from the end of the stroke
        # volume array.
        volume = np.concatenate((stroke_volume, post_volume[1:]))

        computed_pressure = PressureFromVolume(
            volume[::5],
            stroke_pressure[0]*1E5,
            self.reactive_case.experiment_parameters['Tin'],
        ).pressure

        copy('{:.4f}'.format(stroke_pressure[0]))

        self.volout = np.vstack(
            (time[::5] + yaml_data['comptime']/1000, volume[::5])
        ).transpose()

        self.presout = np.vstack(
            (time[:n_print_pts:5] + yaml_data['comptime']/1000,
             print_pressure[::5])
        ).transpose()

        self.write_output(
            self.volout,
            self.presout,
            self.reactive_case.experiment_parameters['Tin'],
        )

        if self.plotting:
            if self.pressure_comparison_figure is None:
                self.pressure_comparison_figure = plt.figure('Pressure Trace Comparison')
                self.pressure_comparison_axis = self.pressure_comparison_figure.add_subplot(1, 1, 1)
                if platform.system() == 'Windows':
                    m = plt.get_current_fig_manager()
                    m.window.showMaximized()

            self.pressure_comparison_axis.cla()
            self.pressure_comparison_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time,
                self.reactive_case.pressure_trace.pressure,
            )
            self.pressure_comparison_axis.plot(time[:n_print_pts:5], print_pressure[::5])
            self.pressure_comparison_axis.plot(time[::5], computed_pressure)
            self.pressure_comparison_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time,
                reactive_line,
            )

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
            choice = input('Are you sure you want to overwrite the {sim_type} simulation?'
                           'Input y or n:'.format(sim_type=sim_type))
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
        compression_time = self.load_yaml()['comptime']

        if self.plotting:
            if self.simulation_figure is None:
                self.simulation_figure = plt.figure('Simulation Comparison')
                self.simulation_axis = self.simulation_figure.add_subplot(1, 1, 1)
                if platform.system() == 'Windows':
                    m = plt.get_current_fig_manager()
                    m.window.showMaximized()

            self.simulation_axis.plot(
                self.presout[:, 0]*1000 - compression_time,
                self.presout[:, 1],
            )

            if self.nonreactive_sim is not None:
                self.simulation_axis.plot(
                    self.nonreactive_sim.time*1000 - compression_time,
                    self.nonreactive_sim.pressure,
                )

            if self.reactive_sim is not None:
                self.simulation_axis.plot(
                    self.reactive_sim.time*1000 - compression_time,
                    self.reactive_sim.pressure,
                )
                der = self.reactive_sim.derivative
                self.simulation_axis.plot(
                    self.reactive_sim.time*1000 - compression_time,
                    der/np.amax(der)*np.amax(self.reactive_sim.pressure),
                )

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
        exp = AltExperiment(file_name)
        if exp.pressure_trace.is_reactive:
            self.reactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_reactive_figures(exp)
        else:
            self.nonreactive_experiments[exp.file_path.name] = exp
            if self.plotting:
                self.plot_nonreactive_figure(exp)


class Simulation(object):
    """Contains a single simulation of the experiment.

    Parameters
    ----------
    initial_temperature : `float`
        The initial temperature of the simulation
    initial_pressure : `float`
        The initial pressure of the simulation
    volume : `numpy.ndarray` or `None`
        The volume trace to be used for the simulation. Must be
        supplied, but if the input value is `None`, the volume trace
        will be read from the file ``volume.csv``. The first column
        should be the time, the second column should be the volume.
    is_reactive : `bool`
        If the simulation should be reactive or non-reactive. If `False`
        sets the Cantera reaction rate multiplier to 0.0 via the
        `~cantera.Kinetics.set_multiplier` function.
    end_temp : `float`, optional
        Reactor temperature at which the simulation will be ended
    end_time : `float`, optional
        Time at which the simulation will be ended
    chem_file : `str`, optional
        String filename of the chemistry file to use

    Attributes
    ----------
    time : `numpy.ndarray`
        Array of simulated time values
    temperature : `numpy.ndarray`
        Array of simulated temperature values
    pressure : `numpy.ndarray`
        Array of simulated pressure values
    input_volume : `numpy.ndarray`
        Array of input volume values
    simulated_volume : `numpy.ndarray`
        Array of simulated volume values
    """

    def __init__(self, initial_temperature, initial_pressure, volume, is_reactive,
                 end_temp=2500., end_time=0.2, chem_file='species.cti'):

        if volume is None:
            volume = np.genfromtxt('volume.csv', delimiter=',')
            keywords = {'vproTime': volume[:, 0], 'vproVol': volume[:, 1]}
        else:
            keywords = {'vproTime': volume[:, 0], 'vproVol': volume[:, 1]}

        self.time = []
        self.temperature = []
        self.pressure = []
        self.input_volume = volume
        self.simulated_volume = []

        gas = ct.Solution(chem_file)
        gas.TP = initial_temperature, initial_pressure
        if not is_reactive:
            gas.set_multiplier(0)
        reac = ct.IdealGasReactor(gas)
        env = ct.Reservoir(ct.Solution('air.xml'))
        ct.Wall(reac, env, A=1.0, velocity=VolumeProfile(keywords))
        netw = ct.ReactorNet([reac])
        netw.set_max_time_step(keywords['vproTime'][1])
        self.time.append(netw.time)
        self.temperature.append(reac.T)
        self.pressure.append(gas.P/1E5)
        self.simulated_volume.append(reac.volume)

        while reac.T < end_temp and netw.time < end_time:
            if cantera_version[1] > 2:
                netw.step()
            else:
                netw.step(1)
            self.time.append(netw.time)
            self.temperature.append(reac.T)
            self.pressure.append(gas.P/1E5)
            self.simulated_volume.append(reac.volume)

        self.time = np.array(self.time)
        self.pressure = np.array(self.pressure)
        self.temperature = np.array(self.temperature)
        self.simulated_volume = np.array(self.simulated_volume)
        self.derivative = self.calculate_derivative(self.pressure, self.time)

    def calculate_derivative(self, dep_var, indep_var):
        """Calculate the derivative.

        Parameters
        ----------
        dep_var : `numpy.ndarray`
            Dependent variable (e.g., the pressure)
        indep_var : `numpy.ndarray`
            Independent variable (e.g., the time)

        Returns
        -------
        `numpy.ndarray`
            1-D array containing the derivative

        Notes
        -----
        The derivative is calculated by computing the first-order
        Lagrange polynomial fit to the point under consideration and
        its nearest neighbors. The Lagrange polynomial is used because
        of the unequal spacing of the simulated data.
        """
        m = len(dep_var)
        ddt = np.zeros(m)
        for i in range(1, m-2):
            x = indep_var[i]
            x_min = indep_var[i-1]
            x_plu = indep_var[i+1]
            y = dep_var[i]
            y_min = dep_var[i-1]
            y_plu = dep_var[i+1]
            ddt[i] = (y_min*(x - x_plu)/((x_min - x)*(x_min - x_plu)) +
                      y*(2*x - x_min - x_plu)/((x - x_min)*(x - x_plu)) +
                      y_plu*(x - x_min)/((x_plu - x_min)*(x_plu - x)))

        return ddt


class Experiment(object):
    """Contains all the information of a single RCM experiment.

    Parameters
    ----------
    file_path : `pathlib.Path`, optional
        If an argument is supplied, it should be an instance of
        `~pathlib.Path`. If no argument is supplied, the
        filename is read from the standard input as a string and
        resolved into a `~pathlib.Path` object.

    Attributes
    ----------
    file_path : `pathlib.Path`
        Object storing the file path
    experiment_parameters : `dict`
        Stores the parameters of the experiment parsed from the
        filename by `~uconnrcmpy.utilities.parse_file_name`.
    voltage_trace : `~uconnrcmpy.traces.VoltageTrace`
        Stores the experimental voltage signal and related traces
    pressure_trace : `~uconnrcmpy.traces.ExperimentalPressureTrace`
        Stores the experimental pressure trace and its derivative
    ignition_delay : `float`
        The overall ignition delay of the experiment. Will be zero for
        a non-reactive experiment.
    first_stage : `float`
        The first stage ignition delay of the experiment. Will be zero
        for a non-reactive experiment of if there is no first-stage
        ignition.
    T_EOC : `float`
        The temperature estimated at the end of compression
    """

    def __init__(self, file_path=None):
        self.resolve_file_path(file_path)
        self.experiment_parameters = self.parse_file_name(self.file_path)
        self.voltage_trace = VoltageTrace(self.file_path)
        self.pressure_trace = ExperimentalPressureTrace(self.voltage_trace,
                                                        self.experiment_parameters['pin'],
                                                        self.experiment_parameters['factor'],
                                                        )
        self.process_pressure_trace()
        self.copy_to_clipboard()

    def parse_file_name(self, file_path):
        """Parse the file name of an experimental trace.

        Parameters
        ----------
        file_path : :class:`pathlib.Path`
            The Path object that contains the path of the current experimental data file.
            The filename associated with the path should be in the following format::

                [NR_]XX_in_YY_mm_ZZZK-AAAAt-BBBx-DD-Mon-YY-Time.txt

            where:

            - ``[NR_]``: Optional non-reactive indicator
            - ``XX``: inches of spacers
            - ``YY``: millimeters of shims
            - ``ZZZ``: Initial temperature in Kelvin
            - ``AAAA``: Initial pressure in Torr
            - ``BBB``: Multiplication factor set on the charge amplifier
            - ``DD-Mon-YY-Time``: Day, Month, Year, and Time of experiment

        Returns
        -------
        :class:`dict`
            Dictionary containing the parameters of the experiment with the
            following names:

            - ``spacers``: Inches of spacers
            - ``shims``: Millimeters of shims
            - ``Tin``: Initial temperature in Kelvin
            - ``pin``: Initial pressure in Torr
            - ``factor``: Multiplication factor set on the charge amplifier
            - ``time_of_day``: Time of day of the experiment
            - ``date``: Date of the experiment with the time
            - ``date_year``: Date of the experiment with the year
        """
        name_parts = {}
        fname = file_path.name.lstrip('NR_')
        fname = fname.rstrip('.txt')
        name_split = fname.split('_')
        name_split_end = name_split[4].split('-', maxsplit=3)
        data_date = datetime.strptime(name_split_end[3], '%d-%b-%y-%H%M')

        name_parts['spacers'] = int(name_split[0])/10
        name_parts['shims'] = int(name_split[2])
        name_parts['Tin'] = int(name_split_end[0][:-1])
        name_parts['pin'] = int(name_split_end[1][:-1])
        name_parts['factor'] = int(name_split_end[2][:-1])
        name_parts['time_of_day'] = data_date.strftime('%H%M')
        name_parts['date'] = data_date.strftime('%d-%b-%H%M')
        name_parts['date_year'] = data_date.strftime('%d-%b-%y')
        return name_parts

    def resolve_file_path(self, file_path=None):
        if file_path is None:
            filename = input('Filename: ')
            try:
                self.file_path = Path(filename).resolve()
            except FileNotFoundError:
                filename += '.txt'
                self.file_path = Path(filename).resolve()
        else:
            self.file_path = file_path.resolve()

    def process_pressure_trace(self):
        if self.pressure_trace.is_reactive:
            self.ignition_delay, self.first_stage = self.calculate_ignition_delay()

            try:
                self.T_EOC = self.calculate_EOC_temperature()
            except RuntimeError as e:
                self.T_EOC = 0
                print('Exception in computing the temperature at EOC', e)
        else:
            self.ignition_delay, self.first_stage, self.T_EOC = 0, 0, 0

    def copy_to_clipboard(self):
        # Copy the relevant information to the clipboard for pasting
        # into a spreadsheet
        copy('\t'.join(map(str, [
            self.experiment_parameters['time_of_day'], self.experiment_parameters['pin'],
            self.experiment_parameters['Tin'], self.pressure_trace.p_EOC, self.ignition_delay,
            self.first_stage, self.T_EOC, self.experiment_parameters['spacers'],
            self.experiment_parameters['shims']])))

    def calculate_ignition_delay(self):
        """Calculate the ignition delay from the pressure trace.
        """
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
        if platform.system() == 'Windows':
            m = plt.get_current_fig_manager()
            m.window.showMaximized()


class AltExperiment(Experiment):
    """Contains all the information of a single alternate RCM experiment
    """
    def __init__(self, file_path=None):
        self.resolve_file_path(file_path)
        self.experiment_parameters = self.parse_file_name(self.file_path)
        self.pressure_trace = AltExperimentalPressureTrace(self.file_path)
        self.process_pressure_trace()
        self.copy_to_clipboard()

    def parse_file_name(self, file_path):
        """Parse the file name of an experimental trace, where the name is in an alternate format.

        Parameters
        ----------
        file_path : :class:`pathlib.Path`
            The Path object that contains the path of the current experimental data file.
            The filename associated with the path should be in the following format::

                'Fuel_FF_EqRatio_P.PP_PercentAr_AR_PercentN2_N2_[NR_]XX_in_YY_mm_ZZZK_AAAAtorr-DD-Mon-YY-Time_Endplug_CC.txt'

            where:

            - ``[NR_]``: Optional non-reactive indicator
            - ``FF``: Fuel symbol
            - ``P.PP``: Equivalence ratio to two decimal places
            - ``AR``: Percent of argon in the oxidizer
            - ``N2``: Percent of nitrogen in the oxidizer
            - ``XX``: inches of spacers
            - ``YY``: millimeters of shims
            - ``ZZZ``: Initial temperature in Kelvin
            - ``AAAA``: Initial pressure in Torr
            - ``BBB``: Multiplication factor set on the charge amplifier
            - ``DD-Mon-YY-Time``: Day, Month, Year, and Time of experiment
            - ``CC``: ``HC`` or ``LC`` to indicate the high clearance or
                      low clearance endplug, respectively

        Returns
        -------
        :class:`dict`
            Dictionary containing the parameters of the experiment with the
            following names:

            - ``spacers``: Inches of spacers
            - ``shims``: Millimeters of shims
            - ``Tin``: Initial temperature in Kelvin
            - ``pin``: Initial pressure in Torr
            - ``time_of_day``: Time of day of the experiment
            - ``date``: Date of the experiment with the time
            - ``date_year``: Date of the experiment with the year
        """
        name_parts = {}
        fname = file_path.name.rstrip('.txt')
        name_split = fname.split('_')
        try:
            name_split.remove('NR')
        except ValueError:
            pass
        name_split_end = name_split[13].split('-', maxsplit=1)
        data_date = datetime.strptime(name_split_end[1], '%d-%b-%y-%H%M')

        name_parts['spacers'] = int(name_split[8])/10
        name_parts['shims'] = int(name_split[10])
        name_parts['Tin'] = int(name_split[12][:-1])
        name_parts['pin'] = int(name_split_end[0].rstrip('torr'))
        name_parts['time_of_day'] = data_date.strftime('%H%M')
        name_parts['date'] = data_date.strftime('%d-%b-%H%M')
        name_parts['date_year'] = data_date.strftime('%d-%b-%y')
        return name_parts


def process_folder(path='.', plot=False):
    """Process a folder of experimental files.

    Process a folder containing files with reactive experiments to
    calculate the ignition delays and copy a table with the results
    to the clipboard.

    Parameters
    ----------
    path : `str`, optional
        Path to folder to be analyzed. Defaults to the current folder.
    plot : `bool`, optional
        True to enable plotting. False by default.
    """
    p = Path(path)
    result = []

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    for f in p.glob('[0-3]*.txt'):
        print(f)
        case = Experiment(f.resolve())
        result.append('\t'.join(map(str, [
            case.experiment_parameters['time_of_day'], case.experiment_parameters['pin'],
            case.experiment_parameters['Tin'], case.pressure_trace.p_EOC, case.ignition_delay,
            case.first_stage, case.T_EOC, case.experiment_parameters['spacers'],
            case.experiment_parameters['shims']])))
        if plot:
            ax.plot(case.pressure_trace.zeroed_time, case.pressure_trace.pressure,
                    label=case.experiment_parameters['date'])

    copy('\n'.join(sorted(result)))
    print('Finished')


def process_alt_folder(path='.', plot=False):
    """Process a folder of alternative experimental files.

    Process a folder containing files with reactive experiments to
    calculate the ignition delays and copy a table with the results
    to the clipboard.

    Parameters
    ----------
    path : `str`, optional
        Path to folder to be analyzed. Defaults to the current folder.
    plot : `bool`, optional
        True to enable plotting. False by default.
    """
    p = Path(path)
    result = []

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    for f in p.glob('Fuel_MF*.txt'):
        if 'NR' in f.name:
            continue
        print(f)
        case = AltExperiment(f.resolve())
        result.append('\t'.join(map(str, [
            case.experiment_parameters['date_year'], case.experiment_parameters['time_of_day'],
            case.experiment_parameters['pin'], case.experiment_parameters['Tin'],
            case.pressure_trace.p_EOC, case.ignition_delay, case.first_stage, case.T_EOC,
            case.experiment_parameters['spacers'], case.experiment_parameters['shims'], f.name])))
        if plot:
            ax.plot(case.pressure_trace.zeroed_time, case.pressure_trace.pressure,
                    label=case.experiment_parameters['date'])

    copy('\n'.join(sorted(result)))
    print('Finished')

if __name__ == '__main__':
    process_folder()
