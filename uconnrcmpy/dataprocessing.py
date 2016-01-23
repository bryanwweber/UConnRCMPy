"""Experiment Module"""

# System imports
from pathlib import Path
from glob import glob

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml
import cantera as ct
from cansen.profiles import VolumeProfile

# Local imports
from .utilities import parse_file_name, copy
from .traces import (VoltageTrace,
                     ExperimentalPressureTrace,
                     TemperatureFromPressure,
                     VolumeFromPressure,
                     PressureFromVolume,
                     )
from .constants import cantera_version


class Condition(object):
    def __init__(self, plotting=True):
        self.reactive_experiments = {}
        self.nonreactive_experiments = {}
        self.reactive_case = None
        self.nonreactive_case = None
        self.presout = None
        self.volout = None
        self.nonreactive_sim = None
        self.reactive_sim = None
        if plotting:
            self.plotting = plotting
            self.all_runs_figure = None
            self.nonreactive_figure = None
            self.pressure_comparison_figure = None
            self.simulation_figure = None

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
        if self.all_runs_figure is None:
            self.all_runs_figure = plt.figure('Reactive Pressure Trace Comparison')
            self.all_runs_axis = self.all_runs_figure.add_subplot(1, 1, 1)
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

    def plot_nonreactive_figure(self, exp):
        if self.nonreactive_figure is None:
            self.nonreactive_figure = plt.figure('Non-Reactive Pressure Trace Comparison')
            self.nonreactive_axis = self.nonreactive_figure.add_subplot(1, 1, 1)
            m = plt.get_current_fig_manager()
            m.window.showMaximized()

            reactive_file = self.load_yaml()['reacfile']
            reactive_parameters = parse_file_name(Path(reactive_file))
            self.reactive_case = self.reactive_experiments[reactive_parameters['date']]
            self.nonreactive_axis.plot(
                self.reactive_case.pressure_trace.zeroed_time,
                self.reactive_case.pressure_trace.pressure,
                label=reactive_parameters['date'],
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
            reactive_parameters = parse_file_name(Path(yaml_data['reacfile']))
            try:
                self.reactive_case = self.reactive_experiments[reactive_parameters['date']]
            except KeyError:
                self.reactive_case = Experiment(Path(yaml_data['reacfile']))
                reactive_case_date = self.reactive_case.experiment_parameters['date']
                self.reactive_experiments[reactive_case_date] = self.reactive_case

        if self.nonreactive_case is None:
            nonreactive_parameters = parse_file_name(Path(yaml_data['nonrfile']))
            try:
                self.nonreactive_case = self.nonreactive_experiments[nonreactive_parameters['date']]
            except KeyError:
                self.nonreactive_case = Experiment(Path(yaml_data['nonrfile']))
                nonreactive_case_date = self.nonreactive_case.experiment_parameters['date']
                self.nonreactive_experiments[nonreactive_case_date] = self.nonreactive_case

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
        `csv` format to the file `volume.csv` and the experimental
        pressure trace is written to the file
        `Tc__P0__T0_XXXK_pressure.txt`, where `XXX` represent the
        initial temperature of the experiment. The user should fill in
        the missing values into the file name after simulations are
        completed and the values are known.
        """
        np.savetxt('volume.csv', volout, delimiter=',')
        np.savetxt('Tc__P0__T0_{}K_pressure.txt'.format(Tin), presout, delimiter='\t')

    def run_simulation(self, run_reactive=False, run_nonreactive=True, end_temp=2500, end_time=0.2):
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


class Simulation(object):
    """Class for simulations of experiments."""

    def __init__(self, initial_temperature, initial_pressure, volume, is_reactive,
                 end_temp=2500., end_time=0.2):

        if volume is None:
            data = np.genfromtxt('volume.csv', delimiter=',')
            keywords = {'vproTime': data[:, 0], 'vproVol': data[:, 1]}
        else:
            keywords = {'vproTime': volume[:, 0], 'vproVol': volume[:, 1]}

        self.time = []
        self.temperature = []
        self.pressure = []
        self.input_volume = volume
        self.simulated_volume = []

        gas = ct.Solution('species.cti')
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
        dep_var : :class:`numpy.ndarray`
            Dependent variable (e.g., the pressure)
        indep_var : :class:`numpy.ndarray`
            Independent variable (e.g., the time)

        Returns
        -------
        :class:`numpy.ndarray`
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


def process_folder(path='.', plot=False):
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
            ax.plot(case.zeroed_time, case.pressure, label=case.experiment_parameters['date'])

    copy('\n'.join(sorted(result)))
    print('Finished')

if __name__ == '__main__':
    process_folder()
