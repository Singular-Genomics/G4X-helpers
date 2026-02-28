from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

from . import constants as c
from . import io

try:
    __version__ = version('g4x_helpers')
except PackageNotFoundError:
    __version__ = 'unknown'

_LAZY_ATTRS = {
    'redemux': ('.main_features', 'redemux'),
    'resegment': ('.main_features', 'resegment'),
    'migrate': ('.main_features', 'migrate'),
    'G4Xoutput': ('.g4x_output', 'G4Xoutput'),
}

__all__ = ['__version__', 'c', 'io', *sorted(_LAZY_ATTRS)]


def __getattr__(name):
    if name not in _LAZY_ATTRS:
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

    module_name, attr_name = _LAZY_ATTRS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(globals().keys() | set(__all__))
