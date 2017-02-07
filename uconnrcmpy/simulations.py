"""Simulations Module"""

# Third party imports
import numpy as np
import cantera as ct


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
    end_temp : `float`
        Reactor temperature at which the simulation will be ended
    end_time : `float`
        Time at which the simulation will be ended
    chem_file : `str`
        String filename of the chemistry file to use
    initial_temperature : `float`
        The initial temperature of the simulation
    initial_pressure : `float`
        The initial pressure of the simulation
    """

    def __init__(self, initial_temperature, initial_pressure, volume, is_reactive,
                 end_temp=2500., end_time=0.2, chem_file='species.cti', cti_source=None):

        if volume is None:
            volume = np.genfromtxt('volume.csv', delimiter=',')

        inp_time = volume[:, 0]
        inp_vol = volume[:, 1]

        self.time = []
        self.temperature = []
        self.pressure = []
        self.input_volume = volume
        self.simulated_volume = []
        self.end_temp = end_temp
        self.end_time = end_time
        self.is_reactive = is_reactive
        self.chem_file = chem_file
        self.initial_temperature = initial_temperature
        self.initial_pressure = initial_pressure

        if cti_source is None:
            gas = ct.Solution(chem_file)
        else:
            gas = ct.Solution(source=cti_source)
        gas.TP = self.initial_temperature, self.initial_pressure
        if not self.is_reactive:
            gas.set_multiplier(0)
        reac = ct.IdealGasReactor(gas)
        env = ct.Reservoir(ct.Solution('air.xml'))
        ct.Wall(reac, env, A=1.0, velocity=VolumeProfile(inp_time, inp_vol))
        netw = ct.ReactorNet([reac])
        netw.set_max_time_step(inp_time[1])
        self.time.append(netw.time)
        self.temperature.append(reac.T)
        self.pressure.append(gas.P/1E5)
        self.simulated_volume.append(reac.volume)

        while reac.T < self.end_temp and netw.time < self.end_time:
            netw.step()
            self.time.append(netw.time)
            self.temperature.append(reac.T)
            self.pressure.append(gas.P/1E5)
            self.simulated_volume.append(reac.volume)

        self.time = np.array(self.time)
        self.pressure = np.array(self.pressure)
        self.temperature = np.array(self.temperature)
        self.simulated_volume = np.array(self.simulated_volume)
        self.derivative = self.calculate_derivative(self.pressure, self.time)

    def __repr__(self):
        return ('Simulation(initial_temperature={self.initial_temperature!r}, '
                'initial_pressure={self.initial_pressure!r}, volume={self.input_volume!r}, '
                'is_reactive={self.is_reactive!r}, end_temp={self.end_temp!r}, '
                'end_time={self.end_time!r}, chem_file={self.chem_file!r})').format(
                    self=self,
                )

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


class VolumeProfile(object):
    """
    Set the velocity of the piston by using a user specified volume
    profile. The initialization and calling of this class are handled
    by the `cantera.Func1` interface of Cantera.

    The velocity is calculated by assuming a unit area and using the
    forward difference, calculated by `numpy.diff`. This function is
    only called once when the class is initialized at the beginning of
    a problem so it is efficient.

    Parameters
    ----------
    time: `numpy.ndarray`
        Array or list of time values
    volume: `numpy.ndarray`
        Array or list of volume values

    Attributes
    ----------
    time: `numpy.ndarray`
        Array of time values
    volume: `numpy.ndarray`
        Array of volume values
    velocity: `numpy.ndarray`
        Array of velocity values
    """

    def __init__(self, time, volume):
        # The time and volume are stored as lists in the keywords
        # dictionary. The volume is normalized by the first volume
        # element so that a unit area can be used to calculate the
        # velocity.
        self.time = np.array(time)
        self.volume = np.array(volume)/volume[0]

        # The velocity is calculated by the forward difference.
        # numpy.diff returns an array one element smaller than the
        # input array, so we append a zero to match the length of the
        # self.time array.
        self.velocity = np.diff(self.volume)/np.diff(self.time)
        self.velocity = np.append(self.velocity, 0)

    def __repr__(self):
        return ('VolumeProfile(time={self.time!r}, volume={self.volume!r}, '
                'velocity={self.velocity!r}').format(
                    self=self,
                )

    def __call__(self, t):
        """Return the velocity when called during a time step.

        Parameters
        ----------
        t : `float`
            Current simulation time.
        """

        if t < self.time[-1]:
            # prev_time_point is the previous value in the time array
            # after the current simulation time
            prev_time_point = self.time[self.time <= t][-1]
            # index is the index of the time array where
            # prev_time_point occurs
            index = np.where(self.time == prev_time_point)[0][0]
            return self.velocity[index]
        else:
            return 0
