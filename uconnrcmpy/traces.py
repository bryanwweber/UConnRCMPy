"""All of the kinds of traces in UConnRCMPy"""

# System imports

# Third-party imports
import numpy as np
import cantera as ct
from scipy import signal as sig
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress

# Local imports
from .constants import (one_atm_in_bar,
                        one_atm_in_torr,
                        one_bar_in_pa,
                        )


class VoltageTrace(object):
    """Voltage signal from a single experiment.

    Parameters
    ----------
    file_path : `pathlib.Path`
        `~pathlib.Path` object associated with the particular experiment

    Attributes
    ----------
    signal : `numpy.ndarray`
        2-D array containing the raw signal from the experimental
        text file. First column is the time, second column is the
        voltage.
    time : `numpy.ndarray`
        The time loaded from the signal trace
    frequency : `int`
        The sampling frequency of the pressure trace
    filtered_voltage : `numpy.ndarray`
        The voltage trace after filtering

    Note
    ----
    The first sample of the voltage is set equal to the
    mean of the first 200 points to eliminate DAQ startup
    effects seen in some data.
    """
    def __init__(self, file_path):
        self.file_path = file_path
        self.signal = np.genfromtxt(str(self.file_path))
        self.time = self.signal[:, 0]
        self.frequency = np.rint(1/self.time[1])

        self.filter_frequency = None

        self.signal[0, 1] = np.mean(self.signal[:200, 1])

        self.filtered_voltage = self.filtering(self.signal[:, 1])

    def __repr__(self):
        return 'VoltageTrace(file_path={self.file_path!r})'.format(self=self)

    @property
    def filter_frequency(self):
        """The cutoff frequency for the low-pass filter

        When setting the frequency, if the ``value`` is `None`,
        determines the optimal cutoff frequency for a first-order
        Butterworth low-pass filter by analyzing the root-mean-squared
        residuals for a sequence of cutoff frequencies. The residuals
        plotted as a function of the cutoff frequency tend to have a
        linear portion for a range of cutoff frequencies. Analysis of
        typical data files from our RCM has shown this range to start
        near ``nyquist_freq*0.05``. The end point is determined by
        looping through values from ``nyquist_freq*0.5`` to
        ``nyquist_freq*0.1`` and finding the location where the
        coefficient of determination of a linear fit is maximized.
        A line is fit to this portion of the residuals curve and
        the intersection point of a horizontal line through the
        y-intercept of the fit and the residuals curve is used to
        determine the optimal cutoff frequency (see Figure 2 in Yu
        et al. [1]_). The methodology is described by Yu et al. [1]_,
        and the code is modifed from Duarte [2]_.

        References
        ----------
        .. [1] B. Yu, D. Gabriel, L. Noble, and K.N. An, "Estimate of
               the Optimum Cutoff Frequency for the Butterworth Low-Pass
               Digital Filter", Journal of Applied Biomechanics, Vol. 15,
               pp. 318-329, 1999.
               DOI: `10.1123/jab.15.3.318 <http://dx.doi.org/10.1123/jab.15.3.318>`_

        .. [2] M. Duarte, "Residual Analysis", v.3 2014/06/13,
               http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/ResidualAnalysis.ipynb
        """

        return self._filter_frequency

    @filter_frequency.setter
    def filter_frequency(self, value):
        if value is None:
            nyquist_freq = self.frequency/2.0
            n_freqs = 101
            freqs = np.linspace(nyquist_freq/n_freqs, nyquist_freq, n_freqs)
            resid = np.zeros(n_freqs)
            for i, fc in enumerate(freqs):
                b, a = sig.butter(1, fc/nyquist_freq)
                yf = sig.filtfilt(b, a, self.signal[:, 1])
                resid[i] = np.sqrt(np.mean((yf - self.signal[:, 1])**2))

            end_points = np.linspace(0.5, 0.1, 9)
            r_sq = np.zeros(len(end_points))
            intercepts = np.zeros(len(end_points))
            for i, end_point in enumerate(end_points):
                # The indices of the frequencies used for fitting the straight line
                fit_freqs = np.arange(np.nonzero(freqs >= nyquist_freq*0.05)[0][0],
                                      np.nonzero(freqs >= nyquist_freq*end_point)[0][0] + 1)
                _, intercepts[i], r, _, _ = linregress(freqs[fit_freqs], resid[fit_freqs])
                r_sq[i] = r**2

            intercept = intercepts[np.argmax(r_sq)]

            # The UnivariateSpline with s=0 forces the spline fit through every
            # data point in the array. The residuals are shifted down by the
            # intercept so that the root of the spline is the optimum cutoff
            # frequency
            try:
                self._filter_frequency = UnivariateSpline(freqs, resid - intercept, s=0).roots()[0]
            except IndexError:
                self._filter_frequency = float(input(
                    'Automatic setting of the filter frequency failed. Please input a frequency; '
                    'typical values are between 1000-5000 Hz: '))
        else:
            self._filter_frequency = value

    def change_filter_freq(self, value):
        """Change the filter frequency

        This method is intended to be used after the `VoltageTrace` has
        already been created. It sets the `~VoltageTrace.filter_frequency`
        attribute and then runs the filtering.

        Parameters
        ----------
        value : `float` or `None`
            Value for the filter frequency. If `None`, an optimal frequency
            will be estimated by the procedure detailed in the documentation
            for the `~VoltageTrace.filter_frequency` attribute.
        """
        self.filter_frequency = value
        self.filtered_voltage = self.filtering(self.signal[:, 1])

    def savetxt(self, filename, **kwargs):
        """Save a text file output of the voltage trace.

        Save a text file with the time in the first column and the filtered
        voltage in the second column. The keyword arguments are the same as
        `numpy.savetxt`.

        Parameters
        ----------
        filename : `str`
            Filename of the output file
        """
        np.savetxt(fname=filename, X=np.vstack(self.time, self.filtered_voltage).T, **kwargs)

    def filtering(self, data):
        """Filter the input using a low-pass filter.

        The filter is a first-order Butterworth filter, with the
        cutoff frequency taken from the instance attribute
        `~VoltageTrace.filter_frequency`.

        Parameters
        ----------
        data : `numpy.ndarray`
            The data that should be filtered

        Returns
        -------
        `numpy.ndarray`
            1-D array of the same length as the input data
        """
        nyquist_freq = self.frequency/2.0
        b, a = sig.butter(1, (self.filter_frequency)/nyquist_freq)
        return sig.filtfilt(b, a, data, padtype='odd', padlen=101, method='pad')


class ExperimentalPressureTrace(object):
    """Pressure trace from a single experiment.

    Parameters
    ----------
    voltage_trace : `VoltageTrace`
        Instance of class containing the voltage trace of the
        experiment.
    initial_pressure_in_torr : `float`
        The initial pressure of the experiment, in units of Torr
    factor : `float`
        The factor set on the charge amplifier

    Attributes
    ----------
    pressure : `numpy.ndarray`
        The pressure trace computed from the filtered voltage trace
    time : `numpy.ndarray`
        A 1-D array containting the time. Copied from
        `VoltageTrace.time`
    frequency : `int`
        Integer sampling frequency of the experiment. Copied from
        `VoltageTrace.frequency`
    p_EOC : `float`
        Pressure at the end of compression
    EOC_idx : `int`
        Integer index in the `pressure` and `time` arrays
        of the end of compression.
    is_reactive : `bool`
        Boolean if the pressure trace represents a reactive or
        or non-reactive experiment
    derivative : `numpy.ndarray`
        1-D array containing the raw derivative computed from the
        `pressure` trace.
    zeroed_time : `numpy.ndarray`
        1-D array containing the time, with the zero point set at
        the end of compression.
    """
    def __init__(self, voltage_trace, initial_pressure_in_torr, factor):
        initial_pressure_in_bar = initial_pressure_in_torr*one_atm_in_bar/one_atm_in_torr
        self.pressure = (voltage_trace.filtered_voltage - voltage_trace.filtered_voltage[0])
        self.pressure *= factor
        self.pressure += initial_pressure_in_bar

        self.raw_pressure = (voltage_trace.signal[:, 1] - voltage_trace.signal[0, 1])
        self.raw_pressure *= factor
        self.raw_pressure += initial_pressure_in_bar

        self.time = voltage_trace.time
        self.frequency = voltage_trace.frequency

        self.p_EOC, self.EOC_idx, self.is_reactive = self.find_EOC()
        self.derivative = self.calculate_derivative(self.pressure, self.time)
        # Smooth the derivative with a moving average 151 points wide
        self.derivative = sig.fftconvolve(self.derivative, np.ones(151)/151, mode='same')
        self.zeroed_time = self.time - self.time[self.EOC_idx]

    def __repr__(self):
        return ('ExperimentalPressureTrace(p_EOC={self.p_EOC!r}, '
                'is_reactive={self.is_reactive!r})').format(self=self)

    def savetxt(self, filename, **kwargs):
        """Save a text file output of the pressure trace.

        Save a text file with the time in the first column and the filtered
        pressure in the second column. The keyword arguments are the same as
        `numpy.savetxt`.

        Parameters
        ----------
        filename : `str`
            Filename of the output file
        """
        np.savetxt(fname=filename, X=np.vstack(self.time, self.pressure).T, **kwargs)

    def pressure_fit(self, comptime=0.08):
        """Fit a line to the pressure trace before compression starts.

        Parameters
        ----------
        comptime : `float`, optional
            Desired compression time, computed from the EOC, to when
            the pressure fit should start

        Returns
        -------
        `numpy.polyfit`
            Numpy object containing the parameters of the fit
        """
        beg_compress = int(np.floor(self.EOC_idx - comptime*self.frequency))
        time = np.linspace(0, (beg_compress - 1)/self.frequency, beg_compress)
        fit_pres = self.pressure[:beg_compress]
        fit_pres[0:9] = fit_pres[10]
        linear_fit = np.polyfit(time, fit_pres, 1)
        return linear_fit

    def change_EOC_time(self, time, is_reactive=True):
        """Change the EOC time for an experiment

        Parameters
        ----------
        time : `float`
            The new value of the EOC time
        is_reactive : `boolean`
            The experiment is reactive or not
        """
        offset = int(round(time/1000*self.frequency, 0))
        self.EOC_idx += offset
        self.p_EOC = self.pressure[self.EOC_idx]
        self.is_reactive = is_reactive
        self.zeroed_time = self.time - self.time[self.EOC_idx]

    def find_EOC(self):
        """Find the index and pressure at the end of compression.

        Returns
        -------
        `tuple`
            Returns a tuple with types (`float`, `int`,
            `bool`) representing the pressure at the end of
            compression, the index of the end of compression relative
            to the start of the pressure trace, and a boolean that is
            True if the case is reactive and False otherwise,
            respectively

        Notes
        -----
        The EOC is found by moving backwards from the maximum pressure
        point and testing the values of the pressure. When the test
        value becomes greater than the previous pressure, we have reached
        the minimum pressure before ignition, in the case of a reactive
        experiment. Then, the EOC is the maximum of the pressure before
        this minimum point. If the pressure at the minimum is close to
        the initial pressure, assume the case is non-reactive and set
        the EOC pressure and the index to the max pressure point.
        """
        is_reactive = True
        max_p = np.amax(self.pressure)
        max_p_idx = np.argmax(self.pressure)
        min_p_idx = max_p_idx - 100
        while self.pressure[min_p_idx] >= self.pressure[min_p_idx - 50]:
            min_p_idx -= 1

        p_EOC = np.amax(self.pressure[0:min_p_idx])
        p_EOC_idx = np.argmax(self.pressure[0:min_p_idx])
        diff = abs(self.pressure[p_EOC_idx] - self.pressure[15])
        if diff < 5.0:
            p_EOC, p_EOC_idx = max_p, max_p_idx
            is_reactive = False

        return p_EOC, p_EOC_idx, is_reactive

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
        The derivative is calculated by a second-order forward method
        and any places where the derivative is infinite are set to
        zero.
        """
        m = len(dep_var)
        ddt = np.zeros(m)
        ddt[:m-2] = (-3*dep_var[:m-2] + 4*(dep_var[1:m-1]) - dep_var[2:m])/(2*np.diff(indep_var[:m-1]))  # NOQA
        ddt[np.isinf(ddt)] = 0
        return ddt


class AltExperimentalPressureTrace(ExperimentalPressureTrace):
    """Process an alternate experimental pressure trace.

    These pressure traces do not have an associated voltage trace,
    but the machinery in the VoltageTrace class is useful for
    filtering.
    """
    def __init__(self, file_path, initial_pressure_in_torr):
        # This is not a real voltage trace
        pressure_trace = VoltageTrace(file_path)

        self.time = pressure_trace.time
        self.frequency = pressure_trace.frequency
        self.filter_frequency = pressure_trace.filter_frequency

        self.pressure = pressure_trace.filtered_voltage
        pressure_start = np.mean(self.pressure[20:500])
        self.pressure -= pressure_start
        self.pressure += initial_pressure_in_torr*one_atm_in_bar/one_atm_in_torr

        self.raw_pressure = pressure_trace.signal[:, 1]
        raw_pressure_start = np.mean(self.raw_pressure[20:500])
        self.raw_pressure -= raw_pressure_start
        self.raw_pressure += initial_pressure_in_torr*one_atm_in_bar/one_atm_in_torr

        self.p_EOC, self.EOC_idx, self.is_reactive = self.find_EOC()
        self.derivative = self.calculate_derivative(self.pressure, self.time)
        self.zeroed_time = self.time - self.time[self.EOC_idx]

    def __repr__(self):
        return ('AltExperimentalPressureTrace(p_EOC={self.p_EOC!r}, '
                'is_reactive={self.is_reactive!r})').format(self=self)


class PressureFromVolume(object):
    """Create a pressure trace given a volume trace.

    Using Cantera to evaluate the thermodynamic properties, compute a
    pressure trace from a volume trace.

    Parameters
    ----------
    volume : `numpy.ndarray`
        1-D array containing the reactor volume
    p_initial : `float`
        Initial pressure of the experiment, in bar
    T_initial : `float`, optional
        Initial temperature of the experiment, in Kelvin.
        Optional for Cantera versions greater than 2.2.0.
    chem_file : `str`, optional
        Filename of the chemistry file to be used

    Attributes
    ----------
    pressure : `numpy.ndarray`
        The pressure trace

    Notes
    -----
    The pressure is computed in a `cantera.Solution` object by
    setting the volume and the entropy according to an isentropic
    process using the given volume trace.
    """
    def __init__(self, volume, p_initial, T_initial, chem_file='species.cti', cti_source=None):
        if cti_source is None:
            gas = ct.Solution(chem_file)
        else:
            gas = ct.Solution(source=cti_source)
        gas.TP = T_initial, p_initial
        initial_volume = gas.volume_mass
        initial_entropy = gas.entropy_mass
        self.pressure = np.zeros((len(volume)))
        for i, v in enumerate(volume):
            gas.SV = initial_entropy, v*initial_volume
            self.pressure[i] = gas.P/one_bar_in_pa

    def __repr__(self):
        return 'PressureFromVolume(pressure={self.pressure!r})'.format(self=self)


class VolumeFromPressure(object):
    r"""Create a volume trace given a pressure trace.

    Using Cantera to evaluate the thermodynamic properties, compute a
    volume trace from a pressure trace.

    Parameters
    ----------
    pressure : `numpy.ndarray`
        1-D array containing the reactor pressure
    v_initial : `float`
        Initial volume of the experiment, in m**3
    T_initial : `float`, optional
        Initial temperature of the experiment, in Kelvin. Optional for
        Cantera versions greater than 2.2.0.
    chem_file : `str`, optional
        Filename of the chemistry file to be used

    Attributes
    ----------
    volume : `numpy.ndarray`
        The volume trace

    Notes
    -----
    The volume is computed according to the formula

    .. math:: v_i = v_{initial}*\rho_{initial}/\rho_i

    where the index :math:`i` indicates the current point. The state
    is set at each point by setting the pressure from the input array
    and the entropy to be constant. The volume is computed by the
    isentropic relationship described above.
    """
    def __init__(self, pressure, v_initial, T_initial, chem_file='species.cti', cti_source=None):
        if cti_source is None:
            gas = ct.Solution(chem_file)
        else:
            gas = ct.Solution(source=cti_source)
        gas.TP = T_initial, pressure[0]*one_bar_in_pa
        initial_entropy = gas.entropy_mass
        initial_density = gas.density
        self.volume = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p*one_bar_in_pa
            self.volume[i] = v_initial*initial_density/gas.density

    def __repr__(self):
        return 'VolumeFromPressure(volume={self.volume!r})'.format(self=self)


class TemperatureFromPressure(object):
    """Create a temperature trace given a pressure trace.

    Using Cantera to evaluate the thermodynamic properties, compute a
    pressure trace from a volume trace.

    Parameters
    ----------
    pressure : `numpy.ndarray`
        1-D array containing the pressure
    T_initial : `float`
        Initial temperature of the experiment, in Kelvin.
        Optional for Cantera versions greater than 2.2.0.
    chem_file : `str`, optional
        Filename of the chemistry file to be used

    Attributes
    ----------
    temperature : `numpy.ndarray`
        The temperature trace

    Notes
    -----
    The temperature is computed in a `cantera.Solution` object by
    setting the pressure and the entropy according to an isentropic
    process using the given pressure trace.
    """
    def __init__(self, pressure, T_initial, chem_file='species.cti', cti_source=None):
        if cti_source is None:
            gas = ct.Solution(chem_file)
        else:
            gas = ct.Solution(source=cti_source)
        gas.TP = T_initial, pressure[0]*one_bar_in_pa
        initial_entropy = gas.entropy_mass
        self.temperature = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p*one_bar_in_pa
            self.temperature[i] = gas.T

    def __repr__(self):
        return 'TemperatureFromPressure(temperature={self.temperature!r})'.format(self=self)
