"""Experiments module"""

# System imports
from datetime import datetime
from pathlib import Path
import platform

# Third-part imports
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from pyperclip import copy
from .traces import (VoltageTrace,
                     ExperimentalPressureTrace,
                     AltExperimentalPressureTrace,
                     TemperatureFromPressure,
                     )


class Experiment(object):
    """Contains all the information of a single RCM experiment.

    Parameters
    ----------
    file_path : `str` or `pathlib.Path`, optional
        If an argument is supplied, it should be an instance of
        `~pathlib.Path` or `str`. If no argument is supplied, the
        filename is read from the standard input as a string and
        resolved into a `~pathlib.Path` object.
    cti_file : `str` or `pathlib.Path`, optional
        Location of the CTI file for Cantera
    cti_source : `str`, optional
        String containing the source of a CTI file for Cantera
    copy : `bool`, optional
        Boolean indicating whether values should be copied to
        the clipboard.

    Attributes
    ----------
    compression_time : `float`
        The compression time in milliseconds, from the end of compression to the
        approximate start of piston motion
    file_path : `pathlib.Path`
        Object storing the file path
    experiment_parameters : `dict`
        Stores the parameters of the experiment parsed from the
        filename by `~uconnrcmpy.Experiment.parse_file_name`.
    voltage_trace : `~uconnrcmpy.traces.VoltageTrace`
        Stores the experimental voltage signal and related traces
    pressure_trace : `~uconnrcmpy.traces.ExperimentalPressureTrace`
        Stores the experimental pressure trace and its derivative
    output_end_time : `float`
        The end time for the output to a file, relative to the end of compression in milliseconds
    offset_points : `float`
        The number of points to offset the end of compression, to better
        match the non-reactive and reactive pressure traces together
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

    def __init__(self, file_path=None, cti_file=None, cti_source=None, copy=True):
        self.resolve_file_path(file_path)
        self.experiment_parameters = self.parse_file_name(self.file_path)
        self.voltage_trace = VoltageTrace(self.file_path)
        self.pressure_trace = ExperimentalPressureTrace(self.voltage_trace,
                                                        self.experiment_parameters['pin'],
                                                        self.experiment_parameters['factor'],
                                                        )
        self.compression_time = None
        self.output_end_time = None
        self.offset_points = None
        if cti_source is None and cti_file is None:
            raise ValueError('One of cti_file or cti_source must be specified')
        elif cti_source is not None and cti_file is not None:
            raise ValueError('Only one of cti_file or cti_source can be specified')
        elif cti_source is None:
            self.cti_file = Path(cti_file).resolve()
            with open(str(cti_file), 'r') as in_file:
                self.cti_source = in_file.read()
        else:
            self.cti_source = cti_source
        self.process_pressure_trace()
        if copy:
            self.copy_to_clipboard()

    def __repr__(self):
        return 'Experiment(file_path={self.file_path!r})'.format(self=self)

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
        name_parts['date'] = data_date.strftime('%d-%b-%H:%M')
        name_parts['date_year'] = data_date.strftime('%d-%b-%y')
        return name_parts

    def resolve_file_path(self, file_path=None):
        if file_path is None:
            self.file_path = Path(input('Filename: '))
        else:
            self.file_path = Path(file_path)
        try:
            self.file_path = self.file_path.resolve()
        except FileNotFoundError:
            self.file_path = self.file_path.with_suffix('.txt')
            self.file_path = self.file_path.resolve()

    def process_pressure_trace(self):
        if self.pressure_trace.is_reactive:
            self.ignition_delay, self.first_stage = self.calculate_ignition_delay()

            try:
                self.T_EOC = self.calculate_EOC_temperature()
            except ct.CanteraError as e:
                self.T_EOC = 0
                print('Exception in computing the temperature at EOC', e)
        else:
            self.ignition_delay, self.first_stage, self.T_EOC = 0, 0, 0

    def change_EOC_time(self, time, is_reactive=True):
        """Change the EOC time for an experiment

        Parameters
        ----------
        time : `float`
            The new value of the EOC time
        is_reactive : `boolean`
            The experiment is reactive or not
        """
        self.pressure_trace.change_EOC_time(time, is_reactive)
        self.process_pressure_trace()
        self.copy_to_clipboard()
        if hasattr(self, 'p_axis'):
            # Create a list of the lines in the axis so that
            # removing one doesn't affect the list as its
            # looping. Use the loop rather than axis.clear()
            # so that the axis labels and legend are preserved.
            for l in list(self.p_axis.lines):
                l.remove()

            self.p_axis.plot(self.pressure_trace.zeroed_time*1000,
                             self.pressure_trace.raw_pressure,
                             'b', label="Raw Pressure")
            self.p_axis.plot(self.pressure_trace.zeroed_time*1000.0,
                             self.pressure_trace.pressure,
                             'g', label="Pressure")

        if hasattr(self, 'dpdt_axis'):
            for l in list(self.dpdt_axis.lines):
                l.remove()

            self.dpdt_axis.plot(self.pressure_trace.zeroed_time*1000.0,
                                self.pressure_trace.derivative/1000.0,
                                'm', label="Derivative")


    def change_filter_freq(self, value):
        """Change the cutoff frequency of the filter for the voltage trace

        Parameters
        ----------
        value : `float`
            The new value of the cutoff frequency
        """
        self.voltage_trace.change_filter_freq(value)
        self.pressure_trace = ExperimentalPressureTrace(self.voltage_trace,
                                                        self.experiment_parameters['pin'],
                                                        self.experiment_parameters['factor'],
                                                        )
        self.process_pressure_trace()
        self.copy_to_clipboard()
        if hasattr(self, 'p_axis'):
            # Create a list of the lines in the axis so that
            # removing one doesn't affect the list as its
            # looping. Use the loop rather than axis.clear()
            # so that the axis labels and legend are preserved.
            for l in list(self.p_axis.lines):
                l.remove()

            self.p_axis.plot(self.pressure_trace.zeroed_time*1000,
                             self.pressure_trace.raw_pressure,
                             'b', label="Raw Pressure")
            self.p_axis.plot(self.pressure_trace.zeroed_time*1000.0,
                             self.pressure_trace.pressure,
                             'g', label="Pressure")

        if hasattr(self, 'dpdt_axis'):
            for l in list(self.dpdt_axis.lines):
                l.remove()

            self.dpdt_axis.plot(self.pressure_trace.zeroed_time*1000.0,
                                self.pressure_trace.derivative/1000.0,
                                'm', label="Derivative")

    def copy_to_clipboard(self):
        """Copy experimental information to the clipboard

        The information is intended to be pasted to a spreadsheet.
        The information is ordered by columns in the following way:

        A. 4-digit time of day
        B. The initial pressure of the experiment, in Torr
        C. The initial temperature of the experiment, in K
        D. The end of compression pressure
        E. The overall ignition delay
        F. The first stage ignition delay
        G. The estimated end of compression temperature
        H. The inches of spacers for the experiment
        I. The millimeters of shims for the experiment
        J. The cutoff frequency that was used to filter the voltage trace
        """
        copy('\t'.join(map(str, [
            self.experiment_parameters['time_of_day'], self.experiment_parameters['pin'],
            self.experiment_parameters['Tin'], self.pressure_trace.p_EOC, self.ignition_delay,
            self.first_stage, self.T_EOC, self.experiment_parameters['spacers'],
            self.experiment_parameters['shims'], self.voltage_trace.filter_frequency])))

    def calculate_ignition_delay(self):
        """Calculate the ignition delay from the pressure trace.
        """
        # tau_points is an offset from the EOC to ensure that if
        # ignition is weak, the peak in dP/dt from the compression
        # stroke is not treated as the ignition event. Define points
        # to start looking for and stop looking for ignition.
        tau_points = int(0.002*self.pressure_trace.frequency)
        start_point = self.pressure_trace.EOC_idx + tau_points
        end_point = self.pressure_trace.EOC_idx + tau_points + 100000

        # The index of the maximum of the derivative trace
        # is taken as the point of ignition
        idx_of_ig = np.argmax(self.pressure_trace.derivative[start_point:end_point])

        # Take the offset into account when calculating the ignition
        # delay. The index of ignition is already relative to zero,
        # so we don't need to subtract any times. Stored in milliseconds
        ignition_delay = self.pressure_trace.time[idx_of_ig + tau_points]*1000

        try:
            end_first_stage = start_point + idx_of_ig - tau_points
            idx_of_first_stage = np.argmax(
                self.pressure_trace.derivative[start_point:end_first_stage]
            )
            first_stage = self.pressure_trace.time[idx_of_first_stage + tau_points]*1000
        except ValueError:
            first_stage = 0.0

        return ignition_delay, first_stage

    def calculate_EOC_temperature(self):
        est_comp_time = 0.03
        comp_time_length = int(est_comp_time*self.pressure_trace.frequency)
        pres_to_temp_start_idx = self.pressure_trace.EOC_idx - comp_time_length
        tempp = self.pressure_trace.pressure[pres_to_temp_start_idx:self.pressure_trace.EOC_idx]
        temperature_trace = TemperatureFromPressure(
            tempp,
            self.experiment_parameters['Tin'],
            chem_file=str(self.cti_file),
        )
        return np.amax(temperature_trace.temperature)

    def plot_pressure_trace(self):
        # Plot the filtered pressure and the derivative
        # on a new figure every time
        self.exp_fig = plt.figure(self.experiment_parameters['date'])
        self.p_axis = self.exp_fig.add_subplot(1, 1, 1)
        raw_line, = self.p_axis.plot(self.pressure_trace.zeroed_time*1000,
                                     self.pressure_trace.raw_pressure,
                                     'b', label="Raw Pressure")
        p_line, = self.p_axis.plot(self.pressure_trace.zeroed_time*1000.0,
                                   self.pressure_trace.pressure,
                                   'g', label="Pressure")
        self.dpdt_axis = self.p_axis.twinx()
        dpdt_line, = self.dpdt_axis.plot(self.pressure_trace.zeroed_time*1000.0,
                                         self.pressure_trace.derivative/1000.0,
                                         'm', label="Derivative")
        self.p_axis.legend([raw_line, p_line, dpdt_line],
                           [l.get_label() for l in [raw_line, p_line, dpdt_line]],
                           loc='best')
        self.p_axis.set_xlabel('Time [ms]')
        self.p_axis.set_ylabel('Pressure [bar]')
        self.dpdt_axis.set_ylabel('Time Derivative of Pressure [bar/ms]')
        if platform.system() == 'Windows':
            m = plt.get_current_fig_manager()
            m.window.showMaximized()


class AltExperiment(Experiment):
    """Contains all the information of a single alternate RCM experiment.
    See the documentation for `Experiment` for attribute descriptions.
    """
    def __init__(self, file_path=None, cti_file=None, cti_source=None, copy=True):
        self.resolve_file_path(file_path)
        self.experiment_parameters = self.parse_file_name(self.file_path)
        self.pressure_trace = AltExperimentalPressureTrace(self.file_path,
                                                           self.experiment_parameters['pin'])
        self.compression_time = None
        self.output_end_time = None
        self.offset_points = None
        if cti_source is None and cti_file is None:
            raise ValueError('One of cti_file or cti_source must be specified')
        elif cti_source is not None and cti_file is not None:
            raise ValueError('Only one of cti_file or cti_source can be specified')
        elif cti_source is None:
            self.cti_file = Path(cti_file).resolve()
            with open(str(cti_file), 'r') as in_file:
                self.cti_source = in_file.read()
        else:
            self.cti_source = cti_source
        self.process_pressure_trace()
        if copy:
            self.copy_to_clipboard()

    def __repr__(self):
        return 'AltExperiment(file_path={self.file_path!r})'.format(self=self)

    def copy_to_clipboard(self):
        """Copy experimental information to the clipboard

        The information is intended to be pasted to a spreadsheet.
        The information is ordered by columns in the following way:

        A. 4-digit time of day
        B. The initial pressure of the experiment, in Torr
        C. The initial temperature of the experiment, in K
        D. The end of compression pressure
        E. The overall ignition delay
        F. The first stage ignition delay
        G. The estimated end of compression temperature
        H. The inches of spacers for the experiment
        I. The millimeters of shims for the experiment
        J. The cutoff frequency that was used to filter the voltage trace
        """
        copy('\t'.join(map(str, [
            self.experiment_parameters['time_of_day'], self.experiment_parameters['pin'],
            self.experiment_parameters['Tin'], self.pressure_trace.p_EOC, self.ignition_delay,
            self.first_stage, self.T_EOC, self.experiment_parameters['spacers'],
            self.experiment_parameters['shims'], self.pressure_trace.filter_frequency])))

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
