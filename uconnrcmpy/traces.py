"""All of the kinds of traces in UConnRCMPy"""

# System imports

# Third-party imports
import numpy as np
import cantera as ct
from scipy import signal as sig

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
    smoothed_voltage : `numpy.ndarray`
        The voltage trace after filtering and smoothing
    """
    def __init__(self, file_path):
        self.signal = np.genfromtxt(str(file_path))

        self.time = self.signal[:, 0]
        self.frequency = np.rint(1/self.time[1])

        self.filtered_voltage = self.filtering(self.signal[:, 1])
        self.smoothed_voltage = self.smoothing(self.filtered_voltage)

    def smoothing(self, data, span=21):
        """Smooth the input using a moving average.

        Parameters
        ----------
        data : `numpy.ndarray`
            The data that should be smoothed
        span : `int`, optional
            The width of the moving average. Should be an odd integer.
            The number of points included in the average on either side
            of the current point is given by ``(span-1)/2``.

        Returns
        -------
        `numpy.ndarray`
            A 1-D array of the same length as the input data.

        Notes
        -----
        This function effects the smoothing by convolving the input
        array with a uniform window array whose values are equal to
        ``1.0/span`` and whose length is equal to ``span``. The
        `~scipy.signal.fftconvolve` function from SciPy is used
        to do the convolution for speed. Since we desire an output
        array of the same length as the input, the first ``(span-1)/2``
        points will have improper values, so these are set equal to the
        value of the average at the point ``(span-1)/2``.
        """
        window = np.ones(span)/span
        output = sig.fftconvolve(data, window, mode='same')
        midpoint = int((span - 1)/2)
        output[:midpoint] = output[midpoint]
        return output

    def filtering(self, data, cutoff_hz=10000):
        """Filter the input using a low-pass filter.

        Parameters
        ----------
        data : `numpy.ndarray`
            The data that should be filtered
        cutoff_hz : `int`, optional
            The cutoff frequency for the filter, in Hz. The default
            value was chosen empirically for a particular set of data
            and may need to be adjusted.

        Returns
        -------
        `numpy.ndarray`
            1-D array of the same length as the input data

        Notes
        -----
        Creates a low-pass filter using the window construction funtion
        `~scipy.signal.firwin`. Applies the filter using the
        `~scipy.signal.fftconvolve` function from the same module
        for speed. Defaults to the `~scipy.signal.blackman` window
        for the filter.
        """
        nyquist_freq = self.frequency/2.0
        n_taps = 2**14
        low_pass_filter = sig.firwin(
            n_taps,
            cutoff_hz/nyquist_freq,
            window='blackman',
        )
        return sig.fftconvolve(data, low_pass_filter, mode='same')


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
        The pressure trace computed from the smoothed and filtered
        voltage trace
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
    smoothed_derivative : `numpy.ndarray`
        1-D array containing the smoothed derivative computed from
        the `derivative`
    zeroed_time : `numpy.ndarray`
        1-D array containing the time, with the zero point set at
        the end of compression.
    """
    def __init__(self, voltage_trace, initial_pressure_in_torr, factor):
        initial_pressure_in_bar = initial_pressure_in_torr*one_atm_in_bar/one_atm_in_torr
        self.pressure = (voltage_trace.smoothed_voltage - voltage_trace.smoothed_voltage[0])
        self.pressure *= factor
        self.pressure += initial_pressure_in_bar

        self.time = voltage_trace.time
        self.frequency = voltage_trace.frequency

        self.p_EOC, self.EOC_idx, self.is_reactive = self.find_EOC()
        self.derivative = self.calculate_derivative(self.pressure, self.time)
        self.smoothed_derivative = voltage_trace.smoothing(self.derivative, span=151)
        self.zeroed_time = self.time - self.time[self.EOC_idx]

    def pressure_fit(self, comptime=0.08):
        """Fit a line to the pressure trace before EOC.

        Fit a line to the part of the pressure trace before compression
        starts.

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
        for i in range(m-2):
            ddt[i] = (-dep_var[i+2] + 4*dep_var[i+1] -
                      3*dep_var[i])/(2*(indep_var[i+1] -
                                        indep_var[i]))
        ddt[np.isinf(ddt)] = 0
        return ddt


class AltExperimentalPressureTrace(ExperimentalPressureTrace):
    """Process an alternate experimental pressure trace.

    These pressure traces do not have an associated voltage trace.
    """
    def __init__(self, file_path, initial_pressure_in_torr):
        self.signal = np.genfromtxt(str(file_path))

        self.time = self.signal[:, 0]
        self.frequency = np.rint(1/self.time[1])

        self.filtered_pressure = VoltageTrace.filtering(self, self.signal[:, 1])
        self.pressure = VoltageTrace.smoothing(self, self.filtered_pressure)
        pressure_start = np.average(self.pressure[20:500])
        self.pressure -= pressure_start
        self.pressure += initial_pressure_in_torr*one_atm_in_bar/one_atm_in_torr

        self.p_EOC, self.EOC_idx, self.is_reactive = self.find_EOC()
        self.derivative = self.calculate_derivative(self.pressure, self.time)
        self.smoothed_derivative = VoltageTrace.smoothing(self, self.derivative, span=21)
        self.zeroed_time = self.time - self.time[self.EOC_idx]


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
    def __init__(self, volume, p_initial, T_initial, chem_file='species.cti'):
        gas = ct.Solution(chem_file)
        gas.TP = T_initial, p_initial
        initial_volume = gas.volume_mass
        initial_entropy = gas.entropy_mass
        self.pressure = np.zeros((len(volume)))
        for i, v in enumerate(volume):
            gas.SV = initial_entropy, v*initial_volume
            self.pressure[i] = gas.P/one_bar_in_pa


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
    def __init__(self, pressure, v_initial, T_initial, chem_file='species.cti'):
        gas = ct.Solution(chem_file)
        gas.TP = T_initial, pressure[0]*one_bar_in_pa
        initial_entropy = gas.entropy_mass
        initial_density = gas.density
        self.volume = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p*one_bar_in_pa
            self.volume[i] = v_initial*initial_density/gas.density


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
    def __init__(self, pressure, T_initial, chem_file='species.cti'):
        gas = ct.Solution(chem_file)
        gas.TP = T_initial, pressure[0]*one_bar_in_pa
        initial_entropy = gas.entropy_mass
        self.temperature = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p*one_bar_in_pa
            self.temperature[i] = gas.T
