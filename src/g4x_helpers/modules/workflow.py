import functools
import logging
import sys

from ..utils import setup_logger


def workflow(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger = kwargs.pop('logger', None)
        if logger is None:
            logger = setup_logger(logger_name=func.__name__, file_logger=False)
            logger.info('No logger provided, using default logger (stream only).')

        logger.info('-' * 10)
        logger.info(f'Initializing {func.__name__} workflow.')

        # Redirect prints → INFO logs
        sys.stdout = LoggerWriter(logger, logging.DEBUG)
        result = func(*args, logger=logger, **kwargs)

        logger.info(f'Completed {func.__name__} workflow.')
        return result

    return wrapper


class LoggerWriter:
    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self._buffer = ''

    def write(self, message):
        message = message.strip()
        if message:  # ignore blank writes
            self.logger.log(self.level, message)

    def flush(self):
        pass  # required for Python's IO interface
