import sys
from .conditions import Condition, AltCondition, process_folder, process_alt_folder
from ._version import __version__

if sys.version_info[0] < 3 and sys.version_info[1] < 4:
    raise Exception('Python 3.4 or greater is required to use this package.')
