"""All of the kinds of traces in UConnRCMPy"""

# System imports

# Third-party imports
import numpy as np
from scipy import signal as sig

# Local imports


class VoltageTrace(object):
    """Class for the voltage trace of an experiment"""

    def __init__(self, file_path):
        self.signal = np.genfromtxt(str(self.file_path))

        self.time = self.signal[:, 0]
        """The time loaded from the signal trace."""
        self.frequency = np.rint(1/self.time[1])
        """The sampling frequency of the pressure trace."""

        self.filtered_voltage = self.filtering(self.signal[:, 1])
        self.smoothed_voltage = self.smoothing(self.filtered_voltage)

    def smoothing(self, data, span=21):
        """
        Smooth the input `data` using a moving average of width `span`.
        """
        window = np.ones(span)/span
        output = sig.fftconvolve(data, window, mode='same')
        midpoint = (span - 1)/2
        output[:midpoint] = output[midpoint]
        return output

    def filtering(self, data, cutoff_hz=10000):
        """
        Filter the input `data` using a low-pass filter with cutoff at 10 kHz
        """
        nyquist_freq = self.frequency/2.0
        n_taps = 2**14
        low_pass_filter = sig.firwin(
            n_taps,
            cutoff_hz/nyquist_freq,
            window='blackman',
        )
        return sig.fftconvolve(data, low_pass_filter, mode='same')
