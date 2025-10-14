from __future__ import annotations

import gzip

# import inspect
import logging
import os
import shutil
import traceback
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING

import glymur
import numpy as np
import rich_click as click

if TYPE_CHECKING:
    from .models import G4Xoutput

DEFAULT_THREADS = max(1, (os.cpu_count() // 2 or 4))


def validate_run_base(run_base):
    """check that all expected outputs are present."""

    run_base = validate_path(path_str=run_base, must_exist=True, is_file_ok=False)

    # TDOD: add more required files to check for
    required_paths = ['run_meta.json', 'single_cell_data/feature_matrix.h5', 'segmentation/segmentation_mask.npz']

    for p in required_paths:
        validate_path(run_base / p, must_exist=True, is_dir_ok=False, is_file_ok=True)

    return run_base


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


def _fail_message(func_name, e, trace_back=False):
    click.echo('')
    click.secho(f'Failed {func_name}:', fg='red', err=True, bold=True)
    if trace_back:
        traceback.print_exc()
    raise click.ClickException(f'{type(e).__name__}: {e}')


def initialize_sample(g4x_dir: Path | str, sample_id: str | None, n_threads: int = DEFAULT_THREADS) -> None:
    from .models import G4Xoutput

    glymur.set_option('lib.num_threads', n_threads)
    try:
        sample = G4Xoutput(run_base=g4x_dir, sample_id=sample_id)
    except Exception as e:
        click.echo('')
        click.secho('Failed to load G4Xoutput:', fg='red', err=True, bold=True)
        raise click.ClickException(f'{type(e).__name__}: {e}')

    return sample


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
        format = '[%(asctime)s | %(name)s | %(levelname)s] %(message)s'
    
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
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        log_path = out_dir / f'{logger_name}.log'

        prior_size = log_path.stat().st_size if log_path.exists() else 0
        if file_mode == 'w':
            prior_size = 0

        fh = logging.FileHandler(log_path, mode=file_mode, encoding='utf-8')
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        session_lines = ['log session created']
        if prior_size > 0:
            session_lines.insert(0, '')
        fh.stream.write('\n'.join(session_lines) + '\n')
        fh.flush()

    return logger


def print_k_v(item, value, gap=2):
    value = '<undefined>' if not value else value
    click.secho(f'{item:<{gap}}', dim=True, nl=False)
    click.secho('- ', dim=True, nl=False)
    click.secho(f'{value}', fg='blue', bold=True)


def symlink_original_files(g4x_out: 'G4Xoutput', out_dir: Path | str) -> None:
    ignore_file_list = ['clustering_umap.csv.gz', 'dgex.csv.gz', 'transcript_panel.csv']
    run_base = Path(g4x_out.run_base).resolve()

    for root, dirs, files in os.walk(run_base):
        rel_root = Path(root).relative_to(run_base)
        if str(rel_root) == 'metrics':
            continue
        dst_root = out_dir / rel_root
        dst_root.mkdir(parents=True, exist_ok=True)

        for f in files:
            if f in ignore_file_list:
                continue
            # src_file = Path(root) / f
            dst_file = dst_root / f
            if dst_file.exists():
                dst_file.unlink()
            dst_file.symlink_to((Path(root) / f).resolve())
