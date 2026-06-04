import logging
from typing import TYPE_CHECKING

from . import io
from .schema import definition as sd

if TYPE_CHECKING:
    from .g4x_output import G4Xoutput
    from .schema.validator import BaseValidator

LOGGER = logging.getLogger(__name__)


def demux(
    smp: 'G4Xoutput',
    manifest: str | None = None,
    out_dir: str | None = None,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
    **kwargs,
) -> None:

    from .modules.demux import demux_raw_features

    log = logger or LOGGER

    out_dir = smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)

    manifest_path = smp.src.Manifest.p if manifest is None else manifest
    manifest_file = _collect_input(manifest_path, sd.Manifest)
    raw_features_file = smp.src.RawFeatures

    tx_table = demux_raw_features(
        raw_features=raw_features_file.load(lazy=True), manifest=manifest_file.parse(), logger=log, **kwargs
    )

    smp.reroute_source(sd.Manifest, out_dir, overwrite=overwrite)
    smp.reroute_source(sd.TxTable, out_dir, overwrite=overwrite)
    tx_table.write_csv(smp.src.TxTable.p, compression='gzip')
    manifest_file.load().write_csv(smp.src.Manifest.p)


def _collect_input(
    path: str,
    validator: 'BaseValidator',
    validate: bool = True,
):
    path_valid = io.pathval.validate_file_path(path)
    in_obj = validator(target_path=path_valid)

    if validate and not in_obj.is_valid:
        raise ValueError(f'Provided {validator.__name__} is not valid!\n{in_obj.report_validation()}')

    return in_obj
