import functools
import logging
from pathlib import Path

# import sys
from typing import TYPE_CHECKING

from .. import logging_utils
from ..io import pathval

if TYPE_CHECKING:
    from ..schema.validator import BaseValidator

LOGGER = logging.getLogger(__name__)
DEFAULT_INPUT = '__g4x_default__'


def collect_input(
    smp,
    path: str,
    validator: 'BaseValidator',
    validate: bool = True,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    if path == DEFAULT_INPUT:
        in_obj = getattr(smp.src, validator.__name__)
        prefix = 'default'
    else:
        path_valid = pathval.validate_file_path(path)
        in_obj = validator(target_path=path_valid)
        prefix = 'provided'

    if validate and not in_obj.is_valid:
        raise ValueError(f'Provided {validator.__name__} is not valid!\nreason: {in_obj.report_validation()}')

    log.debug(f'Using {prefix} {validator.__name__} as input:\n%s%s', logging_utils.PGAP, in_obj.p)
    return in_obj


def prepare_output(
    smp,
    out_dir: str,
    validator: 'BaseValidator',
    overwrite: bool = False,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    out_obj = getattr(smp.out, validator.__name__)
    out_obj.root = Path(out_dir)
    pathval.ensure_parent_dir(out_obj.p)

    if out_obj.path_exists and not overwrite:
        raise RuntimeError(
            f'Operation aborted! {validator.__name__} already exists at:\n{logging_utils.PGAP}{out_obj.p}\nUse overwrite=True to ignore this.',
        )

    suffix = 'overriding existing file' if out_obj.path_exists else 'creating new file'
    log.debug(
        f'Using the following path for {validator.__name__} output ({suffix}):\n%s%s', logging_utils.PGAP, out_obj.p
    )
    return True


def g4x_workflow(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wflow_name = func.__name__.removesuffix('_core')
        logger = kwargs.pop('logger', None)
        if logger is None:
            logger = logging_utils.configure_g4x_logging(level='DEBUG')
            logger.info('No logger provided, using default package logger (stream only).')

        workflow_logger = logger.getChild(f'workflow.{wflow_name}')
        workflow_logger.info('-' * 10)
        workflow_logger.info('Initializing %s workflow.', wflow_name)

        try:
            result = func(*args, logger=workflow_logger, **kwargs)
        except Exception as e:
            msg = f'Exception occurred during {wflow_name} workflow: \n\n{e}'
            workflow_logger.error(msg)
            raise RuntimeError(msg) from e
        else:
            workflow_logger.info('Completed %s workflow.', wflow_name)
            return result

    return wrapper


# def workflow(func):
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         wflow_name = func.__name__.removesuffix('_core')
#         logger = kwargs.pop('logger', None)
#         if logger is None:
#             logger = setup_logger(logger_name=wflow_name, file_logger=False)
#             logger.info('No logger provided, using default logger (stream only).')

#         logger.info('-' * 10)
#         logger.info(f'Initializing {wflow_name} workflow.')

#         # Save the original stdout so we can restore it later
#         original_stdout = sys.stdout
#         sys.stdout = LoggerWriter(logger, logging.DEBUG)

#         try:
#             result = func(*args, logger=logger, **kwargs)
#         except Exception as e:
#             msg = f'Exception occurred during {wflow_name} workflow: \n\n{e}'
#             logger.error(msg)
#             raise RuntimeError(msg) from e
#         else:
#             logger.info(f'Completed {wflow_name} workflow.')
#             return result
#         finally:
#             # Always restore stdout
#             sys.stdout = original_stdout

#     return wrapper


# class OutSchema:
#     def __init__(self, out_dir, subdirs=[]):
#         self.out_dir = validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False, resolve_path=False)
#         self.subdirs = subdirs

#         for subdir in self.subdirs:
#             p = self.out_dir / subdir
#             p.mkdir(exist_ok=True, parents=True)
#             setattr(self, subdir, p)
