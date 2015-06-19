"""
Utilities for UConnRCMPy
"""
# System imports
from datetime import datetime

# Library imports
import win32clipboard


class ParsedFilename(object):
    """
    Parse a filename for experimental conditions.

    Given a filename of the form
    [NR_]XX_in_YY_mm_ZZZK-AAAAt-BBBx-Mon-DD-YY-Time.txt,
    parse the filename for the experimental conditions.
    XX = inches of spacers
    YY = millimeters of shims
    ZZZ = Initial temperature in Kelvin
    AAAA = Initial pressure in Torr
    BBB = Multiplication factor set on the charge amplifier
    Mon-DD-YY-Time = Month, day, year, and time of experiment
    """

    def __init__(self, filename):
        self.fname = filename.lstrip('NR_')
        self.fname = self.fname.rstrip('.txt')
        self.name_split = self.fname.split('_')
        self.spacers = int(self.name_split[0])/10
        self.shims = int(self.name_split[2])
        self.name_split_end = self.name_split[4].split('-', maxsplit=3)
        self.Tin = int(self.name_split_end[0][:-1])
        self.pin = int(self.name_split_end[1][:-1])
        self.factor = int(self.name_split_end[2][:-1])
        self.data_date = datetime.strptime(self.name_split_end[3],
                                           '%d-%b-%y-%H%M')
        self.time_of_day = self.data_date.strftime('%H%M')
        self.date = self.data_date.strftime('%d-%b-%H%M')


def copy(text):
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardText(text)
    win32clipboard.CloseClipboard()
