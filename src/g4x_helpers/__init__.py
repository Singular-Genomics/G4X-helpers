from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('g4x_helpers')
except PackageNotFoundError:
    __version__ = 'unknown'

from .main_features import new_bin as new_bin
from .main_features import redemux as redemux
from .main_features import resegment as resegment
from .main_features import tar_viewer as tar_viewer
from .main_features import update_bin as update_bin
from .models import G4Xoutput as G4Xoutput
