from __future__ import annotations

import logging
import os
import textwrap
from datetime import datetime
from pathlib import Path

import polars as pl

PACKAGE_LOGGER_NAME = __package__ or __name__.split('.')[0]
INDENT = 2 * ' '
PGAP = INDENT + '> '

log = logging.getLogger(__name__)


def configure_logging(
    *,
    logger_name: str = PACKAGE_LOGGER_NAME,
    level: int = logging.INFO,
    stream_log: bool = True,
    file_log: bool = False,
    out_dir: str | None = './g4x_helpers/logs',
    append_time: bool = True,
    file_mode: str = 'a',
    clear_handlers: bool = True,
    # format: str | None = None,
) -> logging.Logger:
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    if clear_handlers:
        logger.handlers.clear()

    level = level.upper() if isinstance(level, str) else level

    # stream_format = '%(g4x_name)s - %(message)s'
    stream_format = '%(asctime)s %(levelshort)1s: %(message)s'
    file_format = '%(asctime)s %(levelname)7s | %(g4x_name)s - %(message)s'

    stream_formatter = G4XFormatter(stream_format, datefmt='%H:%M:%S')
    file_formatter = G4XFormatter(file_format, datefmt='%H:%M:%S')

    if stream_log:
        sh = logging.StreamHandler()
        sh.setLevel(level)
        sh.setFormatter(stream_formatter)
        logger.addHandler(sh)

    if file_log:
        if out_dir is None:
            raise ValueError('out_dir must be provided when file_log=True')
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_path = out_dir / f'g4x_{timestamp}.log' if append_time else out_dir / 'g4x.log'

        fh = logging.FileHandler(log_path, mode=file_mode, encoding='utf-8')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(file_formatter)
        logger.addHandler(fh)

    logger.propagate = False
    return logger


class G4XFormatter(logging.Formatter):
    LEVEL_MAP = {
        'DEBUG': 'D',
        'INFO': 'I',
        'WARNING': 'W',
        'ERROR': 'E',
        'CRITICAL': 'C',
    }

    def format(self, record):
        record.levelshort = f'[{self.LEVEL_MAP.get(record.levelname, "?")}]'
        # record.g4x_name = f'g4x.{record.name.split(".")[-1]}'
        record.g4x_name = f'{record.name}'  # .split(".")[-1]}'
        return super().format(record)


def log_msg_wrapped(header: str, msg: str, *, prefix: str = '    ', level: int = logging.INFO):
    if isinstance(level, str):
        level = getattr(logging, level.upper())

    formatted = textwrap.indent(str(msg), prefix=prefix)
    log.log(level, f'{header}\n%s', formatted)


def log_with_path(
    message: str,
    path: str | list,
    *,
    after_path: str = '',
    level: int = logging.DEBUG,
):
    if isinstance(level, str):
        level = getattr(logging, level.upper())

    paths = path if isinstance(path, (list, tuple)) else [path]

    msg = message
    for p in paths:
        msg += f'\n{PGAP}{p}'

    if after_path:
        msg += f'\n{INDENT}{after_path}'

    log.log(level, msg)


def verbose_to_level(verbose: int) -> int:
    """Convert a verbosity level to a logging level.

    -1: disable logging
     0: warnings and above
     1: info and above
     2+: debug and above
    """
    if verbose < 0:
        return logging.CRITICAL + 1  # effectively disables logging

    levels = [
        logging.WARNING,
        logging.INFO,
        logging.DEBUG,
    ]

    return levels[min(verbose, len(levels) - 1)]


def default_workers(max_workers: int = 16, reserve: int = 1) -> int:
    cpu = os.cpu_count() or 1
    return min(max_workers, max(1, cpu - reserve))


def write_table_to_csv(table: pl.DataFrame | pl.LazyFrame, out_path, source_path: str = None):
    if not isinstance(table, (pl.DataFrame, pl.LazyFrame)):
        raise ValueError(f'Expected a Polars DataFrame or LazyFrame, got {type(table)}')
    out_path = Path(out_path)

    is_lazy = isinstance(table, pl.LazyFrame)
    is_gz = out_path.suffix == '.gz'
    compression = 'gzip' if is_gz else 'uncompressed'

    if source_path is not None:
        if is_lazy and Path(source_path) == out_path:
            table = table.collect()
            is_lazy = False

    if is_lazy:
        table.sink_csv(out_path, compression=compression)
    else:
        table.write_csv(out_path, compression=compression)


def get_image_shape(img_path):
    import glymur
    import tifffile

    if img_path.suffix == '.tiff':
        with tifffile.TiffFile(img_path) as tif:
            series = tif.series[0]
            return series.shape
    elif img_path.suffix == '.jp2':
        return glymur.Jp2k(img_path).shape
    else:
        raise ValueError(f'Unsupported image format: {img_path.suffix}')


def kv_line_gap(key, value, separator=' - ', gap=2):
    value = '<undefined>' if not value else value
    line = f'{key:<{gap}}'
    line += separator
    line += f'{value}'

    return line


def pretty_dict_str(d, separator=' - '):
    max_len = max([len(k) for k in d.keys()])
    msg = ''
    for k, v in d.items():
        msg += kv_line_gap(k, v, separator=separator, gap=max_len) + '\n'
    return msg


def peek(smp):
    import spaceplot as sp

    axs = sp.montage_plot(3, panel_size=4.5, layout='compressed')

    downsample = 4

    img = smp.load_nuclear_image()[::downsample, ::downsample]
    mask = smp.load_segmentation()[::downsample, ::downsample]

    df = smp.src.TxTable.load().sample(100_000, with_replacement=True)

    axs[0].imshow(img, cmap='gray')
    axs[1].scatter(df['x_pixel_coordinate'] / downsample, df['y_pixel_coordinate'] / downsample, s=0.1, alpha=0.5)
    axs[1].invert_yaxis()
    axs[2].imshow(mask)

    axs.layout(ticks=False, margins=0)
    sp.show()
