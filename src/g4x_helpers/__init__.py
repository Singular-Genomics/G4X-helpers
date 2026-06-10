# from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('g4x_helpers')
except PackageNotFoundError:
    __version__ = 'unknown'

from . import constants as c
from . import io, schema
from . import sample_ops as ops
from . import utils as ut
from .g4x_output import G4Xoutput
from .modules import aggregate, demux, migrate, single_cell, viewer

__all__ = [
    '__version__',
    'G4Xoutput',
    'c',
    'io',
    'schema',
    'ut',
    'ops',
    'aggregate',
    'demux',
    'migrate',
    'single_cell',
    'viewer',
]
