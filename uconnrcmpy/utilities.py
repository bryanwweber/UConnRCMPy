"""
Utilities for UConnRCMPy
"""
# System imports
from datetime import datetime

# Third Party imports
import win32clipboard


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
