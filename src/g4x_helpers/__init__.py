from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('g4x_helpers')
except PackageNotFoundError:
    __version__ = 'unknown'

_LAZY_ATTRS = {
    'new_bin': ('.main_features', 'new_bin'),
    'redemux': ('.main_features', 'redemux'),
    'resegment': ('.main_features', 'resegment'),
    'tar_viewer': ('.main_features', 'tar_viewer'),
    'update_bin': ('.main_features', 'update_bin'),
    'G4Xoutput': ('.models', 'G4Xoutput'),
}

__all__ = ['__version__', *sorted(_LAZY_ATTRS)]


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
