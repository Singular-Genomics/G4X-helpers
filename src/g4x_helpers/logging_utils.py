import logging
import textwrap
from datetime import datetime
from pathlib import Path

PACKAGE_LOGGER_NAME = __package__ or __name__.split('.')[0]
LOGGER = logging.getLogger(__name__)
INDENT = 2 * ' '
PGAP = INDENT + '> '


def configure_g4x_logging(
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
    # if format is None:
    stream_format = '%(g4x_name)s - %(message)s'
    file_format = '%(asctime)s %(levelname)7s | %(g4x_name)s - %(message)s'

    stream_formatter = G4XFormatter(stream_format)
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
    def format(self, record):
        # record.g4x_name = f'g4x.{record.name.split(".")[-1]}'
        record.g4x_name = f'{record.name}'  # .split(".")[-1]}'
        return super().format(record)


def log_msg_wrapped(
    header: str, msg: str, logger: logging.Logger | None = None, *, prefix: str = '    ', level: int = logging.INFO
):
    log = logger or LOGGER

    if isinstance(level, str):
        level = getattr(logging, level.upper())

    formatted = textwrap.indent(str(msg), prefix=prefix)
    log.log(level, f'{header}\n%s', formatted)


def log_with_path(
    message: str,
    path: str | list,
    logger: logging.Logger | None = None,
    *,
    after_path: str = '',
    level: int = logging.DEBUG,
):
    log = logger or LOGGER

    if isinstance(level, str):
        level = getattr(logging, level.upper())

    paths = path if isinstance(path, (list, tuple)) else [path]

    msg = message
    for p in paths:
        msg += f'\n{PGAP}{p}'

    if after_path:
        msg += f'\n{INDENT}{after_path}'

    log.log(level, msg)
