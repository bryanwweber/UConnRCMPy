"""
Utilities for UConnRCMPy
"""
# System imports
from datetime import datetime

# Third Party imports
import win32clipboard


def parse_file_name(file_path):
    """Parse the file name of an experimental trace.

    Parameters
    ----------
    file_path : Path
        The Path object that contains the path of the current experimental data file.
        The filename associated with the path should be in the following format::

            [NR_]XX_in_YY_mm_ZZZK-AAAAt-BBBx-DD-Mon-YY-Time.txt

        where:

        - `[NR_]`: Optional non-reactive indicator
        - `XX`: inches of spacers
        - `YY`: millimeters of shims
        - `ZZZ`: Initial temperature in Kelvin
        - `AAAA`: Initial pressure in Torr
        - `BBB`: Multiplication factor set on the charge amplifier
        - `DD-Mon-YY-Time`: Day, Month, Year, and Time of experiment

    Returns
    -------
    dict
        Dictionary containing the parameters of the experiment with the
        following names::

        - spacers: Inches of spacers
        - shims: Millimeters of shims
        - Tin: Initial temperature in Kelvin
        - pin: Initial pressure in Torr
        - factor: Multiplication factor set on the charge amplifier
        - time_of_day: Time of day of the experiment
        - date: Date of the experiment

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
    return name_parts


class ParsedFilename(object):
    """
    Parse a filename for experimental conditions.
    """

    def __init__(self, file_path):
        """
        Given a filename of the form
        `[NR_]XX_in_YY_mm_ZZZK-AAAAt-BBBx-DD-Mon-YY-Time.txt,`
        parse the filename for the experimental conditions.

        * `[NR_] = Optional non-reactive indicator
        * `XX` = inches of spacers
        * `YY` = millimeters of shims
        * `ZZZ` = Initial temperature in Kelvin
        * `AAAA` = Initial pressure in Torr
        * `BBB` = Multiplication factor set on the charge amplifier
        * `DD-Mon-YY-Time` = Day, Month, year, and time of experiment
        """
        fname = file_path.name.lstrip('NR_')
        """Processed filename."""
        fname = fname.rstrip('.txt')
        self.name_split = fname.split('_')
        """Split list of filename parts."""
        self.spacers = int(self.name_split[0])/10
        """Inches of spacers."""
        self.shims = int(self.name_split[2])
        """Millimeters of shims."""
        self.name_split_end = self.name_split[4].split('-', maxsplit=3)
        """Split list of filename parts."""
        self.Tin = int(self.name_split_end[0][:-1])
        """Initial temperature in Kelvin."""
        self.pin = int(self.name_split_end[1][:-1])
        """Initial pressure in torr."""
        self.factor = int(self.name_split_end[2][:-1])
        """Multiplication factor from charge amplifier."""
        self.data_date = datetime.strptime(self.name_split_end[3],
                                           '%d-%b-%y-%H%M')
        """Datetime object storing the processed date."""
        self.time_of_day = self.data_date.strftime('%H%M')
        """Time of day of the experiment."""
        self.date = self.data_date.strftime('%d-%b-%H%M')
        """Date of the experiment."""


def copy(text):
    """Copy the input `text` to the Windows clipboard."""
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardText(text)
    win32clipboard.CloseClipboard()
