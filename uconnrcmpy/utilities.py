"""
Utilities for UConnRCMPy
"""
# System imports
from datetime import datetime
import platform

# Third Party imports
if platform.system() == 'Windows':
    import win32clipboard
elif platform.system() == 'Darwin':
    import subprocess


def parse_file_name(file_path):
    """Parse the file name of an experimental trace.

    Parameters
    ----------
    file_path : Path
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
    dict
        Dictionary containing the parameters of the experiment with the
        following names:

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


def copy(text):
    """Copy the input ``text`` to the system clipboard."""
    if platform.system() == 'Windows':
        win32clipboard.OpenClipboard()
        win32clipboard.EmptyClipboard()
        win32clipboard.SetClipboardText(text)
        win32clipboard.CloseClipboard()
    elif platform.system() == 'Darwin':
        process = subprocess.Popen(
            'pbcopy', env={'LANG': 'en_US.UTF-8'}, stdin=subprocess.PIPE)
        process.communicate(text.encode('utf-8'))
