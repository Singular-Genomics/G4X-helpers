import functools
import logging

from ..utils import setup_logger


def workflow(func):
    """Decorator to apply workflow initialization logic."""

    @functools.wraps(func)
    def wrapper(
        logger: logging.Logger | None = None,
        **kwargs,
    ):
        func_name = func.__name__
        if logger is None:
            logger = setup_logger(logger_name=func_name, file_logger=False)
            logger.info('No logger provided, using default logger (stream only).')

        logger.info(f'{"-" * 10}')
        logger.info(f'Initializing {func_name} workflow.')

        result = func(
            logger=logger,
            **kwargs,
        )

        logger.info(f'Completed {func_name} workflow.')
        return result

    return wrapper
