from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("g4x_helpers")
except PackageNotFoundError:
    __version__ = "unknown"

import g4x_helpers.plotting as pl
import g4x_helpers.utils as ut
from g4x_helpers.models import G4Xoutput
