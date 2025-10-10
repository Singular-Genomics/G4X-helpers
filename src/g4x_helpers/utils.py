from __future__ import annotations

import gzip
import logging
import os
import shutil
import zipfile
from pathlib import Path

import glymur
import numpy as np
import rich_click as click

from .models import G4Xoutput

DEFAULT_THREADS = max(1, (os.cpu_count() // 2 or 4))


def validate_path(path_str, must_exist=True, is_dir_ok=True, is_file_ok=True):
    path = Path(path_str)  # .expanduser().resolve()

    if must_exist and not path.exists():
        raise FileNotFoundError(f'Path does not exist: {path}')

    if path.exists():
        if path.is_dir() and not is_dir_ok:
            raise ValueError(f'Expected a file but got a directory: {path}')
        if path.is_file() and not is_file_ok:
            raise ValueError(f'Expected a directory but got a file: {path}')

    return path


def initialize_sample(sample_dir: Path | str, sample_id: str | None, n_threads: int = DEFAULT_THREADS) -> None:
    glymur.set_option('lib.num_threads', n_threads)
    try:
        sample_dir = validate_path(sample_dir, must_exist=True, is_dir_ok=True, is_file_ok=False)
        sample = G4Xoutput(run_base=sample_dir, sample_id=sample_id, n_threads=n_threads)
    except Exception as e:
        click.echo('')
        click.secho('Failed to load G4Xoutput:', fg='red', err=True, bold=True)
        raise click.ClickException(f'{type(e).__name__}: {e}')

    click.echo(sample)


# region file operations
def npzGetShape(npz, key):
    """Takes a path to an .npz file and key to stored array and returns array shape
    without loading the array into memory

    Parameters
        npz: str
        key: str

    Raises: KeyError

    Returns: tuple

    Reference: http://bit.ly/2qsSxy8
    """
    with zipfile.ZipFile(npz) as archive:
        name = '{}.npy'.format(key)
        if name in archive.namelist():
            npy = archive.open(name)
            version = np.lib.format.read_magic(npy)
            shape, _, _ = np.lib.format._read_array_header(npy, version)
            return shape
        else:
            raise KeyError('{} not in archive'.format(key))


def delete_existing(outfile: str | Path) -> None:
    outfile = Path(outfile)
    if outfile.exists() or outfile.is_symlink():
        outfile.unlink()


def gzip_file(outfile: str | Path, remove_original: bool = False) -> None:
    outfile = Path(outfile)
    delete_existing(outfile.with_suffix('.csv.gz'))
    with open(outfile, 'rb') as f_in:
        with gzip.open(outfile.with_suffix('.csv.gz'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if remove_original:
        os.remove(outfile)


def create_custom_out(
    out_dir: Path | str,
    parent_folder: Path | str | None = None,
    file_name: str | None = None,
) -> Path:
    custom_out = Path(out_dir) / parent_folder

    # Ensure the directory exists
    custom_out.mkdir(parents=True, exist_ok=True)

    # Prepare output file path
    outfile = custom_out / file_name

    delete_existing(outfile)

    return outfile


# region logger
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
