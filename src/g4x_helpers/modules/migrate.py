import logging

from .. import io, schema
from .. import utils as ut
from ..roi import Roi
from ..schema.legacy import migrators as mig

log = logging.getLogger(__name__)


def migrate_raw_data(
    legacy_dir: str,
    out_dir: str,
    roi_coords: tuple | None = None,
    protein_subset: list | None = None,
) -> None:

    legacy_dir = io.pathval.validate_dir_path(legacy_dir)
    out_dir = io.pathval.validate_dir_path(out_dir)

    if (out_dir / 'sample.g4x').exists():
        raise FileExistsError(
            'Output directory already contains a sample.g4x file. Aborting migration to prevent overwriting existing data.'
        )

    ut.log_with_path('Starting migration for:', legacy_dir.resolve(), level='INFO')

    roi = None
    if roi_coords is not None:
        roi = Roi(xlims=(roi_coords[0], roi_coords[2]), ylims=(roi_coords[1], roi_coords[3]))
        roi_sz_um = roi.width_um, roi.height_um
        x0, x1, y0, y1 = roi.extent
        log.info(f'Roi provided with xlims={x0, x1}, ylims={y0, y1}, size in um: {roi_sz_um}')

    basic_migrators, roi_migrators = _gather_migrators(legacy_dir)

    if not all(m.is_migratable for m in basic_migrators + roi_migrators):
        msg = status(legacy_dir)
        error = 'Not all source files are migratable. Aborting migration!'
        ut.log_msg_wrapped(f'{error}\n', msg, level='ERROR')
        raise mig.MigrationError(error)

    for m in basic_migrators:
        m.migrate(out_dir)

    for m in roi_migrators:
        m.migrate(out_dir, roi=roi, protein_subset=protein_subset)

    ut.log_with_path('Migration completed! Raw data is available at:', path=out_dir, level='INFO')


def status(
    legacy_dir,
) -> None:

    basic_migrators, roi_migrators = _gather_migrators(legacy_dir)
    migrators = basic_migrators + roi_migrators

    res = {}
    for m in migrators:
        icon = '✓' if m.is_migratable else '✗'
        status = 'is migratable' if m.is_migratable else 'can not be migrated'
        res[f'{icon} {m._name}'] = status

    return ut.pretty_dict_str(res, separator=' ')


def _gather_migrators(legacy_dir):
    sg4x = mig.SampleG4X_Migrator(root=legacy_dir)
    basic_migrators = [sg4x]
    roi_migrators = []
    smp_meta = sg4x.build_sample_g4x()

    basic_migrators.extend(
        [
            mig.SampleSheet_Migrator(root=legacy_dir),
            mig.Manifest_Migrator(root=legacy_dir),
            mig.QCSummary_Migrator(root=legacy_dir),
            mig.Metrics_Migrator(root=legacy_dir),
        ]
    )
    roi_migrators = [
        mig.HnEDir_Migrator(root=legacy_dir),
        mig.Segmentation_Migrator(root=legacy_dir),
        mig.BeadMask_Migrator(root=legacy_dir),
        mig.RawFeatures_Migrator(root=legacy_dir),
        mig.TxTable_Migrator(root=legacy_dir),
    ]

    assay_type = schema.ut.detect_assay_type(smp_meta)
    if assay_type != 'tx_only':
        roi_migrators.append(mig.Protein_Migrator(root=legacy_dir))

    return basic_migrators, roi_migrators
