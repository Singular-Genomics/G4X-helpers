import logging
from datetime import datetime
from pathlib import Path

PACKAGE_LOGGER_NAME = __package__ or __name__.split('.')[0]


def configure_g4x_logging(
    *,
    logger_name: str = PACKAGE_LOGGER_NAME,
    level: int = logging.INFO,
    stream_log: bool = True,
    file_log: bool = False,
    out_dir: str | Path | None = None,
    file_mode: str = 'a',
    clear_handlers: bool = True,
    format: str | None = None,
) -> logging.Logger:
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    if clear_handlers:
        logger.handlers.clear()

    if format is None:
        format = '%(asctime)s | %(g4x_name)s | %(levelname)7s - %(message)s'

    formatter = G4XFormatter(format, datefmt='%H:%M:%S')

    if stream_log:
        sh = logging.StreamHandler()
        sh.setLevel(level)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    if file_log:
        if out_dir is None:
            raise ValueError('out_dir must be provided when file_log=True')
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_path = out_dir / f'g4x_{timestamp}.log'

        fh = logging.FileHandler(log_path, mode=file_mode, encoding='utf-8')
        fh.setLevel(level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    logger.propagate = False
    return logger


class G4XFormatter(logging.Formatter):
    def format(self, record):
        record.g4x_name = f'g4x.{record.name.split(".")[-1]}'
        return super().format(record)
