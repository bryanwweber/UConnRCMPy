import sys
if sys.version_info[0] < 3 and sys.version_info[1] < 4:
    raise Exception('Python 3.4 or greater is required to use this package.')
