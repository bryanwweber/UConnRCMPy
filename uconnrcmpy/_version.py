__version_info__ = ('2', '1', '0')
__version__ = '.'.join(__version_info__[:3])
if len(__version_info__) == 4:
    __version__ += __version_info__[-1]
