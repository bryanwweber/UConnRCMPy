"""
Utilities for UConnRCMPy
"""
# System imports
import platform

# Third Party imports
if platform.system() == 'Windows':
    import win32clipboard
elif platform.system() == 'Darwin':
    import subprocess


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
