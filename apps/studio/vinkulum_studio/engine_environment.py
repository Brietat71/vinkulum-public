"""Keep external engines independent of the desktop's bundled libraries."""

import sys

from PySide6.QtCore import QProcessEnvironment


def external_engine_environment():
    environment = QProcessEnvironment.systemEnvironment()
    for key in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV"):
        environment.remove(key)
    if getattr(sys, "frozen", False):
        # PyInstaller prepends the desktop's library directory on Linux.
        # Restore the caller's loader path for independently installed engines.
        original = environment.value("LD_LIBRARY_PATH_ORIG")
        if original:
            environment.insert("LD_LIBRARY_PATH", original)
        else:
            environment.remove("LD_LIBRARY_PATH")
    return environment
