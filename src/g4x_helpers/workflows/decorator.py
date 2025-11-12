import functools

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

        result = func(*args, logger=logger, **kwargs)

        logger.info(f'Completed {func.__name__} workflow.')
        return result

    return wrapper
