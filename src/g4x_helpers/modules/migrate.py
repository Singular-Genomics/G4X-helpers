import logging
from typing import Literal

from .. import io, schema
from .. import logging_utils as logut
from .. import sample_ops as ops
from .. import utils as ut
from ..g4x_output import G4Xoutput
from ..roi import Roi
from ..schema.legacy import migrators as mig

log = logging.getLogger(__name__)


def migrate_sample(
    sample_dir: str,
    out_dir: str,
    roi_coords: tuple | None = None,
    protein_subset: list | None = None,
    downstream: bool = True,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
) -> None:

    sample_dir = io.pathval.validate_dir_path(sample_dir)
    out_dir = io.pathval.validate_dir_path(out_dir)

    if (out_dir / 'sample.g4x').exists():
        raise FileExistsError(
            'Output directory already contains a sample.g4x file. Aborting migration to prevent overwriting existing data.'
        )

    logut.log_with_path('Starting migration for:', sample_dir, level='INFO')

    roi = None
    if roi_coords is not None:
        roi = Roi(xlims=(roi_coords[0], roi_coords[2]), ylims=(roi_coords[1], roi_coords[3]))
        roi_sz_um = roi.width_um, roi.height_um
        x0, x1, y0, y1 = roi.extent
        log.info(f'Roi provided with xlims={x0, x1}, ylims={y0, y1}, size in um: {roi_sz_um}')

    basic_migrators, roi_migrators = _gather_migrators(sample_dir)

    if not all(m.is_migratable for m in basic_migrators + roi_migrators):
        log.error('Not all migrators are migratable. Aborting migration.')
        status(sample_dir)
        return

    for m in basic_migrators:
        m.migrate(out_dir)

    for m in roi_migrators:
        m.migrate(out_dir, roi=roi, protein_subset=protein_subset)

    log.info('All migrators completed migration. Starting post-processing...')

    smp = G4Xoutput(smp_dir=out_dir)

    if downstream:
        ops.aggregate(smp, out_dir=out_dir, overwrite=True, backend=backend)
        ops.sc_process(smp, out_dir=out_dir, overwrite=True, backend=backend)
        ops.viewer_zarr(smp, out_dir=out_dir, overwrite=True)

    logut.log_msg_wrapped(header='Migration completed. Migrated data is available at\n', msg=smp, level='INFO')


def status(
    sample_dir,
) -> None:

    basic_migrators, roi_migrators = _gather_migrators(sample_dir)
    migrators = basic_migrators + roi_migrators

    res = {}
    for m in migrators:
        icon = '✓' if m.is_migratable else '✗'
        status = 'is migratable' if m.is_migratable else 'can not be migrated'
        res[f'{icon} {m._name}'] = status

    print(ut.pretty_dict_str(res, separator=' '))


def _gather_migrators(sample_dir):
    sg4x = mig.SampleG4X_Migrator(root=sample_dir)
    basic_migrators = [sg4x]
    roi_migrators = []
    smp_meta = sg4x.build_sample_g4x()

    basic_migrators.extend(
        [
            mig.SampleSheet_Migrator(root=sample_dir),
            mig.Manifest_Migrator(root=sample_dir),
            mig.QCSummary_Migrator(root=sample_dir),
            mig.Metrics_Migrator(root=sample_dir),
        ]
    )
    roi_migrators = [
        mig.HnEDir_Migrator(root=sample_dir),
        mig.Segmentation_Migrator(root=sample_dir),
        mig.BeadMask_Migrator(root=sample_dir),
        mig.RawFeatures_Migrator(root=sample_dir),
        mig.TxTable_Migrator(root=sample_dir),
    ]

    assay_type = schema.ut.detect_assay_type(smp_meta)
    if assay_type != 'tx_only':
        roi_migrators.append(mig.Protein_Migrator(root=sample_dir))

    return basic_migrators, roi_migrators
