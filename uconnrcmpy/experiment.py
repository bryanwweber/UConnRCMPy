"""Experiment Module"""

# System imports
from pathlib import Path

# Third-party imports
import numpy as np

# Local imports
from .utilities import parse_file_name
from .traces import (VoltageTrace,)


class Experiment(object):
    """Class containing a single RCM experiment"""

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
