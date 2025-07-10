from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('g4x_helpers')
except PackageNotFoundError:
    __version__ = 'unknown'

from g4x_helpers.models import G4Xoutput as G4Xoutput
