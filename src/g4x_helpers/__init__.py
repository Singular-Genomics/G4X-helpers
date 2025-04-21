from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("g4x_helpers")
except PackageNotFoundError:
    __version__ = "unknown"

import g4x_helpers.plotting
from g4x_helpers.models import G4Xoutput
