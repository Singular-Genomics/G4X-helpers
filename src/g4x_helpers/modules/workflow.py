import functools
import logging
import sys

from .. import logging_utils

# from ..utils import LoggerWriter, setup_logger, validate_path


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
