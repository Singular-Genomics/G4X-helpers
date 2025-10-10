from __future__ import annotations

import logging
import os
from pathlib import Path

import rich_click as click

from .models import G4Xoutput


def verbose_to_log_level(verbose: int) -> int:
    """
    returns a logging level based on verbose integer:
    0 == WARNING
    1 == REPORT
    2 == INFO
    any other integer == DEBUG
    """
    mapping = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG,
    }
    return mapping.get(verbose, logging.DEBUG)


def setup_logger(
    logger_name: str,
    *,
    stream_logger: bool = True,
    stream_level: int = 2,
    file_logger: bool = True,
    file_level: int = 2,
    file_mode: str = 'a',
    out_dir: Path | str | None = None,
    format: str | None = None,
    clear_handlers: bool = True,
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

    stream_level = verbose_to_log_level(stream_level)
    file_level = verbose_to_log_level(file_level)

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    if format is None:
        format = '[%(asctime)s | %(levelname)s | %(name)s] %(message)s'
    formatter = logging.Formatter(format, datefmt='%Y-%m-%d %H:%M:%S')

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


def initialize_sample(sample_dir: Path | str, sample_id: str | None) -> None:
    try:
        sample = G4Xoutput(run_base=sample_dir, sample_id=sample_id)
    except Exception as e:
        click.echo('')
        click.secho('Failed to load G4Xoutput:', fg='red', err=True, bold=True)
        raise click.ClickException(f'{type(e).__name__}: {e}')

    click.echo(sample)
