import functools
import logging
import shutil
import sys

from ..utils import setup_logger


def workflow(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        func_name = func.__name__.removesuffix('_core')
        logger = kwargs.pop('logger', None)
        if logger is None:
            logger = setup_logger(logger_name=func_name, file_logger=False)
            logger.info('No logger provided, using default logger (stream only).')

        logger.info('-' * 10)
        logger.info(f'Initializing {func_name} workflow.')

        # Save the original stdout so we can restore it later
        original_stdout = sys.stdout
        sys.stdout = LoggerWriter(logger, logging.DEBUG)

        try:
            result = func(*args, logger=logger, **kwargs)
        except Exception as e:
            # Delete out_dir on failure, but only if present
            out_dir = kwargs.get('out_dir')
            if out_dir is not None:
                logger.error(f"Exception occurred — deleting out_dir '{out_dir}'.")
                shutil.rmtree(out_dir, ignore_errors=True)
            # Re-raise so caller sees the exception
            raise RuntimeError(f'Error during {func_name} workflow: {e}') from e
        else:
            logger.info(f'Completed {func_name} workflow.')
            return result
        finally:
            # Always restore stdout
            sys.stdout = original_stdout

    return wrapper


class LoggerWriter:
    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self._buffer = ''

    def write(self, message):
        # stdout writes can be chunked, so we buffer and split on newlines
        if not isinstance(message, str):
            message = str(message)

        self._buffer += message
        while '\n' in self._buffer:
            line, self._buffer = self._buffer.split('\n', 1)
            line = line.strip()
            if line:
                self.logger.log(self.level, line)

    def flush(self):
        # Flush any remaining text in the buffer
        if self._buffer.strip():
            self.logger.log(self.level, self._buffer.strip())
        self._buffer = ''
