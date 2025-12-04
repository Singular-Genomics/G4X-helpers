from __future__ import annotations

import gzip
import logging
import os
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd
import polars as pl

DEFAULT_THREADS = max(1, (os.cpu_count() // 2 or 4))
PROBE_PATTERN = r'^(.*?)-([ACGT]{2,30})-([^-]+)$'

primer_read_map = {
    'SP1': 1,
    'm7a': 2,
    'm9a': 3,
    'm6a': 4,
    'm8a': 5,
    'm3a': 6,
}


def validate_path(path_str, must_exist=True, is_dir_ok=True, is_file_ok=True, resolve_path=False) -> Path:
    path = Path(path_str)  # .expanduser().resolve()
    if resolve_path:
        path = path.resolve()

    if must_exist and not path.exists():
        raise FileNotFoundError(f'Path does not exist: {path}')

    if path.exists():
        if path.is_dir() and not is_dir_ok:
            raise ValueError(f'Expected a file but got a directory: {path}')
        if path.is_file() and not is_file_ok:
            raise ValueError(f'Expected a directory but got a file: {path}')

    if not path.exists() and not must_exist:
        if not is_file_ok:
            path.mkdir(parents=True, exist_ok=True)
        else:
            parent = path.parent
            parent.mkdir(parents=True, exist_ok=True)

    return path


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


def write_csv_gz(df, path, keep_uncompressed: bool = False):
    if isinstance(df, pl.LazyFrame):
        df.sink_csv(path)
    if isinstance(df, pl.DataFrame):
        df.write_csv(path)
    if isinstance(df, pd.DataFrame):
        df.to_csv(path, index=False)

    gzip_file(path, remove_original=not keep_uncompressed)


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


class TruncatingFormatter(logging.Formatter):
    def __init__(self, *args, name_max=20, **kwargs):
        super().__init__(*args, **kwargs)
        self.name_max = name_max

    def format(self, record):
        if len(record.name) > self.name_max:
            end = self.name_max - len('...' + 'XX')
            record.name = record.name[:end] + '...' + record.name[-2:]
        # else:
        #     fill = '.' * (self.name_max - len(record.name))
        #     record.name = record.name + fill

        return super().format(record)


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
        # format = '[%(asctime)s | %(name)s | %(levelname)s: %(message)s'
        format = '%(asctime)s | %(name)s | %(levelname)7s - %(message)s'

    formatter = logging.Formatter(format, datefmt='%H:%M:%S')
    # formatter = TruncatingFormatter(format, name_max=16, datefmt='%H:%M:%S')

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

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_path = out_dir / f'{logger_name}_{timestamp}.log'

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


def symlink_original_files(g4x_obj, out_dir: Path | str) -> None:
    ignore_file_list = ['clustering_umap.csv.gz', 'dgex.csv.gz', 'transcript_panel.csv']
    data_dir = Path(g4x_obj.data_dir).resolve()

    for root, dirs, files in os.walk(data_dir):
        rel_root = Path(root).relative_to(data_dir)
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


def parse_input_manifest(file_path: Path, verbose: bool = False) -> pl.DataFrame:
    """
    Parse strings in `df[col]` of the form '<prefix>-<15mer>-<primer>'
    and add columns: gene_name, sequence, primer, primer_code, read.
    Rows that don't match the pattern get nulls in the new columns.
    """
    # print(f'Parsing transcript manifest: {file_path.name}')
    df = pl.read_csv(file_path)

    col = 'probe_name'
    if col not in df.columns:
        raise ValueError(f"transcript manifest must contain column '{col}'")

    parsed = df.with_columns(
        [
            pl.col(col).str.extract(PROBE_PATTERN, 1).alias('gene_name'),
            pl.col(col).str.extract(PROBE_PATTERN, 2).alias('sequence'),
            pl.col(col).str.extract(PROBE_PATTERN, 3).alias('primer'),
        ]
    )

    null_count = parsed.null_count()['sequence'][0]
    if null_count > 0:
        if verbose:
            print(f'{null_count} probes with invalid sequence format will be ignored:')
            null_seqs = parsed.filter(pl.col('sequence').is_null())['probe_name'].to_list()
            for ns in null_seqs:
                print(f'- {ns}')

    if 'panel_type' in df.columns:
        parsed = parsed.drop('panel_type')

    if 'read' in df.columns:
        if verbose:
            print("Using 'read' column provided by input manifest.")
        parsed = parsed.drop('read')
    else:
        plist = parsed['primer'].unique().to_list()
        ign_primer = [p for p in plist if p not in primer_read_map]
        if len(ign_primer) > 0:
            if verbose:
                print('Warning: the following primer names are not known and will be ignored:')
                for ip in ign_primer:
                    print(f'- {ip}')
        parsed = parsed.filter(pl.col('primer').is_in(primer_read_map.keys())).with_columns(
            pl.col('primer').replace(primer_read_map).cast(pl.Int8).alias('read')
        )

    if 'gene_name' in df.columns:
        if verbose:
            print("Using 'gene_name' column provided by input manifest.")
        parsed = parsed.drop('gene_name')

    manifest = parsed.join(df, on=col, how='left')

    manifest = manifest.with_columns(
        probe_type=(
            pl.when(pl.col('gene_name').str.to_lowercase() == 'gdna')
            .then(pl.lit('GCP'))
            .when(pl.col('gene_name').str.starts_with('NCS-'))
            .then(pl.lit('NCS'))
            .when(pl.col('gene_name').str.starts_with('NCP-'))
            .then(pl.lit('NCP'))
            .otherwise(pl.lit('targeting'))
        )
    )

    return manifest
