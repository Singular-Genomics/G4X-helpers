import os
import logging
from pathlib import Path
from typing import List, Tuple, Union, Iterator, Iterable, Any, Generator, Optional, Literal, Callable


def setup_logger(
    logger_name: str,
    *,
    stream_logger: Optional[bool] = True,
    stream_level: Optional[int] = logging.DEBUG,
    file_logger: Optional[bool] = False,
    file_level: Optional[int] = logging.DEBUG,
    file_mode: Optional[str] = 'a',
    out_dir: Optional[Union[Path, str]] = None,
    format: Optional[str] = None,
    clear_handlers: Optional[bool] = False
) -> logging.Logger:
    """
    Sets up a logger with configurable stream and file handlers.

    Parameters
    ----------
    logger_name : str
        Name of the logger.
    stream_logger : bool, optional
        Whether to enable logging to the console (default is True).
    stream_level : int, optional
        Logging level for the stream handler (default is logging.DEBUG).
    file_logger : bool, optional
        Whether to enable logging to a file (default is False).
    file_level : int, optional
        Logging level for the file handler (default is logging.DEBUG).
    file_mode : str, optional
        File mode for the log file. Common options: 'a' for append, 'w' for overwrite (default is 'a').
    out_dir : Path or str, optional
        Directory where the log file will be saved. Required if `file_logger` is True.
    format : str, optional
        Custom log message format. If not provided, a default format will be used.
    clear_handlers : bool, optional
        Whether to clear existing handlers from the logger before adding new ones (default is False).

    Returns
    -------
    logging.Logger
        Configured logger instance.
    
    """
    
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    if format is None:
        format = '[%(asctime)s %(levelname)s %(funcName)s %(name)s] %(message)s'
    formatter = logging.Formatter(format)

    ## optionally clear existing handlers
    if clear_handlers:
        logger.handlers.clear()

    if stream_logger:
        h = logging.StreamHandler()
        h.setLevel(stream_level)
        h.setFormatter(formatter)
        logger.addHandler(h)

    if file_logger:
        assert out_dir is not None, "out_dir must be provided if file_logger is True"
        os.makedirs(out_dir, exist_ok= True)
        fh = logging.FileHandler(f'{out_dir}/{logger_name}.log', mode= file_mode)
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger