from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # This import is only for type checkers (mypy, PyCharm, etc.), not at runtime
    from g4x_helpers.models import G4Xoutput

import logging
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def verbose_to_log_level(verbose: int) -> int:
    """
    returns a logging level based on verbose integer:

    0 == WARNING

    1 == REPORT

    2 == INFO

    any other integer == DEBUG
    """

    if verbose == 0:
        log_level = logging.WARNING
    elif verbose == 1:
        log_level = logging.REPORT
    elif verbose == 2:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG

    return log_level


def setup_logger(
    logger_name: str,
    *,
    stream_logger: bool = True,
    stream_level: int = logging.DEBUG,
    file_logger: bool = False,
    file_level: int = logging.DEBUG,
    file_mode: str = 'a',
    out_dir: Path | str | None = None,
    format: str | None = None,
    clear_handlers: bool = False,
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
        assert out_dir is not None, 'out_dir must be provided if file_logger is True'
        os.makedirs(out_dir, exist_ok=True)
        fh = logging.FileHandler(f'{out_dir}/{logger_name}.log', mode=file_mode)
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger