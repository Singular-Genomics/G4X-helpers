from __future__ import annotations

import gzip
import os
import shutil
import zipfile
from pathlib import Path

import numpy as np


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
