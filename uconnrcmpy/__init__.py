import sys
from .dataprocessing import Condition, AltCondition, process_folder, process_alt_folder

if sys.version_info[0] < 3 and sys.version_info[1] < 4:
    raise Exception('Python 3.4 or greater is required to use this package.')
