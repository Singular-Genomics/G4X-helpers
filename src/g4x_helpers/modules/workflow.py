import functools
import logging
from pathlib import Path
from typing import TYPE_CHECKING

from .. import logging_utils as logut
from ..io import pathval

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput
    from ..schema.validator import BaseValidator

LOGGER = logging.getLogger(__name__)
PRESET_SOURCE = '__g4x_out_tree__'


def collect_input(
    smp: 'G4Xoutput',
    path: str,
    validator: 'BaseValidator',
    validate: bool = True,
    logger: logging.Logger | None = None,
) -> 'BaseValidator':
    log = logger or LOGGER

    if path == PRESET_SOURCE:
        in_obj = getattr(smp.out, validator.__name__)
        prefix = 'pre-set'
    else:
        path_valid = pathval.validate_file_path(path)
        in_obj = validator(target_path=path_valid)
        prefix = 'provided'

    if validate and not in_obj.is_valid:
        raise ValueError(f'Provided {validator.__name__} is not valid!\nreason: {in_obj.report_validation()}')

    logut.log_with_path(f'Using {prefix} {validator.__name__} as input:', in_obj.p, logger=log)
    return in_obj


def reroute_source(
    smp: 'G4Xoutput',
    out_dir: str,
    validator: 'BaseValidator',
    overwrite: bool = False,
    logger: logging.Logger | None = None,
) -> bool:
    log = logger or LOGGER

    out_obj = getattr(smp.out, validator.__name__)
    out_obj.root = Path(out_dir)
    pathval.ensure_parent_dir(out_obj.p)

    if out_obj.path_exists() and not overwrite:
        raise RuntimeError(
            f'Operation aborted! {validator.__name__} already exists at:\n{logut.PGAP}{out_obj.p}\nUse overwrite=True to ignore this.',
        )

    suffix = 'overriding existing file' if out_obj.path_exists() else 'creating new file'
    logut.log_with_path(f'Using the following path for {validator.__name__} output ({suffix}):', out_obj.p, logger=log)

    return True


def g4x_workflow(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wflow_name = func.__name__

        log = kwargs.pop('logger', None)
        if log is None:
            log = logut.configure_g4x_logging(
                logger_name=f'{func.__module__}.{func.__name__}',
                level='INFO',
            )
            # log.info('No logger provided, using default (stream only).')

        log.info('-' * 10)
        log.info('Initializing workflow: %s', wflow_name)

        try:
            result = func(*args, logger=log, **kwargs)
        except Exception:
            log.exception('Exception occurred during workflow: %s', wflow_name)
            raise

        log.info('Completed workflow: %s', wflow_name)
        return result

    return wrapper
