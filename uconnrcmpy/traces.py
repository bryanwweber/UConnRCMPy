"""All of the kinds of traces in UConnRCMPy"""

# System imports

# Third-party imports
import numpy as np
import cantera as ct
from scipy import signal as sig

# Local imports
from .constants import (cantera_version,
                        one_atm_in_bar,
                        one_atm_in_torr,
                        one_bar_in_pa,
                        )


class VoltageTrace(object):
    """Voltage signal from a single experiment.

    Parameters
    ----------
    file_path : :class:`pathlib.Path`
        :class:`pathlib.Path` object associated with the particular experiment

    Attributes
    ----------
    signal : :class:`numpy.ndarray`
        2-D array containing the raw signal from the experimental
        text file. First column is the time, second column is the
        voltage.
    time : :class:`numpy.ndarray`
        The time loaded from the signal trace
    frequency : :class:`int`
        The sampling frequency of the pressure trace
    filtered_voltage : :class:`numpy.ndarray`
        The voltage trace after filtering
    smoothed_voltage : :class:`numpy.ndarray`
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
        data : :class:`numpy.ndarray`
            The data that should be smoothed
        span : :class:`int`, optional
            The width of the moving average. Should be an odd integer.
            The number of points included in the average on either side
            of the current point is given by ``(span-1)/2``.

        Returns
        -------
        :class:`numpy.ndarray`
            A 1-D array of the same length as the input data.

        Notes
        -----
        This function effects the smoothing by convolving the input
        array with a uniform window array whose values are equal to
        ``1.0/span`` and whose length is equal to ``span``. The
        :func:`~scipy.signal.fftconvolve` function from SciPy is used
        to do the convolution for speed. Since we desire an output
        array of the same length as the input, the first ``(span-1)/2``
        points will have improper values, so these are set equal to the
        value of the average at the point ``(span-1)/2``.
        """
        window = np.ones(span)/span
        output = sig.fftconvolve(data, window, mode='same')
        midpoint = (span - 1)/2
        output[:midpoint] = output[midpoint]
        return output

    def filtering(self, data, cutoff_hz=10000):
        """Filter the input using a low-pass filter.

        Parameters
        ----------
        data : :class:`numpy.ndarray`
            The data that should be filtered
        cutoff_hz : :class:`int`, optional
            The cutoff frequency for the filter, in Hz. The default
            value was chosen empirically for a particular set of data
            and may need to be adjusted.

        Returns
        -------
        :class:`numpy.ndarray`
            1-D array of the same length as the input data

        Notes
        -----
        Creates a low-pass filter using the window construction funtion
        :func:`~scipy.signal.firwin`. Applies the filter using the
        :func:`~scipy.signal.fftconvolve` function from the same module
        for speed. Defaults to the :func:`~scipy.signal.blackman` window
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
    voltage_trace : :class:`VoltageTrace`
        Instance of class containing the voltage trace of the
        experiment.
    initial_pressure_in_torr : :class:`float`
        The initial pressure of the experiment, in units of Torr
    factor : :class:`float`
        The factor set on the charge amplifier

    Attributes
    ----------
    pressure : :class:`numpy.ndarray`
        The pressure trace computed from the smoothed and filtered
        voltage trace
    time : :class:`numpy.ndarray`
        A 1-D array containting the time. Copied from
        :class:`VoltageTrace.time`
    frequency : :class:`int`
        Integer sampling frequency of the experiment. Copied from
        :class:`VoltageTrace.frequency`
    p_EOC : :class:`float`
        Pressure at the end of compression
    EOC_idx : :class:`int`
        Integer index in the :class:`pressure` and :class:`time` arrays
        of the end of compression.
    is_reactive : :class:`bool`
        Boolean if the pressure trace represents a reactive or
        or non-reactive experiment
    derivative : :class:`numpy.ndarray`
        1-D array containing the raw derivative computed from the
        :class:`pressure` trace.
    smoothed_derivative : :class:`numpy.ndarray`
        1-D array containing the smoothed derivative computed from
        the :class:`derivative`
    zeroed_time : :class:`numpy.ndarray`
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
        comptime : :class:`float`, optional
            Desired compression time, computed from the EOC, to when
            the pressure fit should start

        Returns
        -------
        :class:`numpy.polyfit`
            Numpy object containing the parameters of the fit
        """
        beg_compress = np.floor(self.EOC_idx - comptime*self.frequency)
        time = np.linspace(0, (beg_compress - 1)/self.frequency, beg_compress)
        fit_pres = self.pressure[:beg_compress]
        fit_pres[0:9] = fit_pres[10]
        linear_fit = np.polyfit(time, fit_pres, 1)
        return linear_fit

    def find_EOC(self):
        """Find the index and pressure at the end of compression.

        Returns
        -------
        :class:`tuple`
            Returns a tuple with types (:class:`float`, :class:`int`,
            :class:`bool`) representing the pressure at the end of
            compression, the index of the end of compression relative
            to the start of the pressure trace, and a boolean that is
            True if the case is reactive and False otherwise,
            respectively

        Notes
        -----
        The EOC is found by moving backwards from the maximum pressure
        point and testing the values of the pressure. When the test
        value becomes less than the previous pressure, we have reached
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
        while self.pressure[min_p_idx] >= self.pressure[min_p_idx - 100]:
            min_p_idx -= 1

        p_EOC = np.amax(self.pressure[0:min_p_idx])
        p_EOC_idx = np.argmax(self.pressure[0:min_p_idx])
        diff = abs(self.pressure[p_EOC_idx] - self.pressure[15])
        if diff < 5:
            p_EOC, p_EOC_idx = max_p, max_p_idx
            is_reactive = False

        return p_EOC, p_EOC_idx, is_reactive

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


class SimulatedPressureTrace(object):
    """The pressure trace from a simulation.

    Parameters
    ----------
    filename : :class:`str`, optional
        The filename to load the data from.
    data : :class:`numpy.recarray`, optional
        The array to load the data from.

    Attributes
    ----------
    data : :class:`numpy.recarray`
        The data as loaded from the file or the input array
    pres : :class:`numpy.ndarray`
        The pressure
    time : :class:`numpy.ndarray`
        The time
    dpdt : :class:`numpy.ndarray`
        The derivative

    Notes
    -----
    By default, the pressure trace is loaded from the file with the
    name ``export.csv``. If a :class:`numpy.recarray` is passed to the
    ``data`` argument, it will be used instead of loading from a file.
    The header in the csv file, or the names of the columns in the
    record array are expected to be in the format

    * ``Pressure_(bar)`` for the pressure
    * ``Time_(sec)`` for the time
    """

    def __init__(self, filename='export.csv', data=None):
        if data is None:
            self.data = np.genfromtxt(filename, delimiter=',', names=True)
        else:
            self.data = data

        self.pres = self.data['Pressure_(bar)']
        self.time = self.data['Time_(sec)']

        self.dpdt = self.derivative(self.pres, self.time)

    def derivative(self, dep_var, indep_var):
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


class PressureFromVolume(object):
    """Create a pressure trace given a volume trace.

    Using Cantera to evaluate the thermodynamic properties, compute a
    pressure trace from a volume trace.

    Parameters
    ----------
    volume : :class:`numpy.ndarray`
        1-D array containing the reactor volume
    p_initial : :class:`float`
        Initial pressure of the experiment, in bar
    T_initial : :class:`float`, optional
        Initial temperature of the experiment, in Kelvin.
        Optional for Cantera versions greater than 2.2.0.
    chem_file : :class:`str`, optional
        Filename of the chemistry file to be used

    Attributes
    ----------
    pressure : :class:`numpy.ndarray`
        The pressure trace

    Notes
    -----
    The pressure is computed in a :class:`cantera.Solution` object by
    setting the volume and the entropy according to an isentropic
    process using the given volume trace.
    """
    def __init__(self, volume, p_initial, T_initial=None, chem_file='species.cti'):
        gas = ct.Solution(chem_file)
        if cantera_version[1] > 2:
            gas.DP = 1.0/volume[0], p_initial*one_bar_in_pa
        elif T_initial is None:
            raise RuntimeError("T_initial must be provided for this version of Cantera.")
        else:
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
    pressure : :class:`numpy.ndarray`
        1-D array containing the reactor pressure
    v_initial : :class:`float`
        Initial volume of the experiment, in m**3
    T_initial : :class:`float`, optional
        Initial temperature of the experiment, in Kelvin. Optional for
        Cantera versions greater than 2.2.0.
    chem_file : :class:`str`, optional
        Filename of the chemistry file to be used

    Attributes
    ----------
    volume : :class:`numpy.ndarray`
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
    def __init__(self, pressure, v_initial, T_initial=None, chem_file='species.cti'):
        gas = ct.Solution(chem_file)
        if cantera_version[1] > 2:
            gas.DP = 1.0/v_initial, pressure[0]*one_bar_in_pa
        elif T_initial is None:
            raise RuntimeError("T_initial must be provided for this version of Cantera.")
        else:
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
    pressure : :class:`numpy.ndarray`
        1-D array containing the pressure
    T_initial : :class:`float`
        Initial temperature of the experiment, in Kelvin.
        Optional for Cantera versions greater than 2.2.0.
    chem_file : :class:`str`, optional
        Filename of the chemistry file to be used

    Attributes
    ----------
    temperature : :class:`numpy.ndarray`
        The temperature trace

    Notes
    -----
    The temperature is computed in a :class:`cantera.Solution` object by
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
