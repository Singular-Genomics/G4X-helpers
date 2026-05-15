import json
import logging
import re
import shutil

import numpy as np
import polars as pl

from .. import c, io
from . import definition as sd
from .validator import BaseValidator, DirectoryValidator, FileValidator

LOGGER = logging.getLogger(__name__)


class DataMigrator(BaseValidator):
    IS_OPTIONAL = False
    VERSION_VALIDATORS = {}
    VERSION_CLASS_PATTERN = re.compile(r'.*_V\d+$')
    COPY_CURRENT = True

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        auto_validators = {}
        for name, obj in cls.__dict__.items():
            if isinstance(obj, type) and issubclass(obj, BaseValidator) and cls.VERSION_CLASS_PATTERN.match(name):
                key = name[0] + name[1:]

                auto_validators[key] = obj

        explicit_validators = getattr(cls, 'VERSION_VALIDATORS', {})

        cls.VERSION_VALIDATORS = {
            **auto_validators,
            **explicit_validators,
        }

    def _current_validator(self):
        for cls in type(self).__mro__:
            if cls is DataMigrator:
                continue

            if issubclass(cls, BaseValidator):
                return cls(root=self.root, format=self.format)

        raise TypeError(f'Could not infer "current" validator for {self._name}')

    def versions(self):
        versions = {'current': self._current_validator()}
        versions.update(
            {name: validator_cls(root=self.root) for name, validator_cls in type(self).VERSION_VALIDATORS.items()}
        )
        return versions

    @property
    def _name(self):
        return type(self).__name__.removesuffix('_Migrator')

    @property
    def target_validator(self):
        return type(self).__bases__[1]

    @property
    def valid_versions(self):
        return list({k for k, v in self.versions().items() if v.is_valid})

    @property
    def mig_version(self):
        if len(self.valid_versions) == 0:
            if self.IS_OPTIONAL:
                return 'missing_optional'
            else:
                raise ValueError(f'No valid versions found for {self._name} in {self.root}')
        if 'current' in self.valid_versions:
            mig_version = 'current'
        else:
            sorted_legacy_versions = sorted(self.valid_versions, key=lambda x: int(x.split('_V')[-1]))
            mig_version = sorted_legacy_versions[-1]
        return mig_version

    @property
    def migrator(self):
        return self.versions()[self.mig_version]

    @property
    def is_migratable(self):
        return self.migration_status[0]

    @property
    def migration_status(self):
        if self.IS_OPTIONAL:
            return True, f'{self._name} is optional, migration not required.'

        has_legacy = len(self.valid_versions) > 0
        if not has_legacy:
            msg = f'{self._name} has no valid versions detected.'
        else:
            msg = f'{self._name} has migratable versions: {self.valid_versions}.'

        return has_legacy, msg

    @property
    def is_file(self):
        return isinstance(self, FileValidator)

    @property
    def is_folder(self):
        return isinstance(self, DirectoryValidator)

    def copy_if_current(self, out_path):
        if self.is_file:
            file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)
            shutil.copy(self.migrator.p, file_out)
        elif self.is_folder:
            shutil.copytree(self.migrator.p, out_path / self.DEFAULT_TARGET_PATH)
        else:
            raise TypeError(f'Cannot copy {self._name} because it is neither a file nor a folder.')

    def migrate(self, out_path, logger: logging.Logger | None = None, *args, **kwargs):
        log = logger or LOGGER
        kwargs['logger'] = log
        log.info(f'Initializing migration of {self._name}')

        if not self.valid_versions and self.IS_OPTIONAL:
            log.warning(f'No valid versions of {self._name} found. Skipping migration of optional data.')
            return

        if not self.is_migratable:
            raise ValueError(f'{self._name} is not migratable. {self.migration_status[1]}')

        out_path = io.pathval.validate_dir_path(out_path)

        if self.COPY_CURRENT and 'current' in self.valid_versions:
            log.debug(f'{self._name} has correct schema, copying file without migration.')
            self.copy_if_current(out_path)
            return

        log.debug(f'Detected legacy version: "{self.mig_version}"')

        migrate_method = getattr(self, '_migrate_method', False)
        if not migrate_method:
            raise NotImplementedError(f'{self._name} does not have a _migrate_method defined.')
        else:
            try:
                result = migrate_method(out_path, *args, **kwargs)
                target = self.target_validator(root=out_path, format=self.format)
                if target.is_valid:
                    log.debug(f'✓ Migration successful for {self._name}!')
                else:
                    raise ImportError(f'Migrated {self._name} did not pass validation! {target.validation()}.')
                return result

            except Exception as e:
                raise ImportError(f'Error migrating {self._name} from {self.target_path}!\nreason: {e}') from None


class SampleG4X_Migrator(DataMigrator, sd.SampleG4X):
    class RunMeta_V1(sd.FileValidator):
        DEFAULT_TARGET_PATH = 'run_meta.json'

        KEYS = [
            'machine',
            'run_id',
            'platform',
            'fc',
            'lane',
            'time_of_creation',
            'transcript_panel',
            'protein_panel',
            'software',
            'software_version',
        ]

    @property
    def smp_sheet(self):
        return SampleSheet_Migrator(root=self.root)

    @property
    def legacy_smp_id(self):
        option_a = 'metrics/transcript_core_metrics.csv'
        option_b = 'metrics/core_metrics.csv'

        if (self.root / option_a).exists():
            option = option_a
        elif (self.root / option_b).exists():
            option = option_b
        else:
            raise FileNotFoundError(f'Neither {option_a} nor {option_b} found in {self.root}')

        smp_id = None
        try:
            smp_id = pl.read_csv(self.root / option)['sample_id'].item()
        except Exception as e:
            raise ValueError(f'Error reading sample_id from {option} in {self.root}. Reason: {e}') from None
        return smp_id

    @property  # this overrides the default behaviour
    def migration_status(self):
        if not self.smp_sheet.is_valid:
            return False, f'{self._name} SampleSheet is not valid.'

        if self.legacy_smp_id is None:
            return False, f'{self._name} Could not determine legacy sample_id.'

        if not len(self.valid_versions) > 0:
            return False, f'{self._name} No valid legacy versions are available.'

        return True, f'{self._name} is migratable.'

    def _migrate_method(self, out_path, **kwargs):

        with open(self.migrator.p, 'r') as f:
            run_meta = json.load(f)

        if run_meta['protein_panel'] == []:
            run_meta['protein_panel'] = None

        return io.create_sample_g4x(
            sample_id=self.legacy_smp_id,
            run_meta=run_meta,
            ssheet=self.smp_sheet.p,
            out_path=out_path / self.DEFAULT_TARGET_PATH,
        )


class QCSummary_Migrator(DataMigrator, sd.QCSummary):
    IS_OPTIONAL = True
    COPY_CURRENT = False

    def _migrate_method(self, out_path, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.p.name)
        shutil.copy(self.p, file_out)


class SampleSheet_Migrator(DataMigrator, sd.SampleSheet):
    def _migrate_method(self, out_path, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)
        shutil.copy(self.p, file_out)


class Segmentation_Migrator(DataMigrator, sd.Segmentation):
    COPY_CURRENT = False

    class Segmentation_V1(sd.Segmentation):
        DEFAULT_TARGET_PATH = 'segmentation/segmentation_mask.npz'

    def _migrate_method(self, out_path: str, roi=None, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if roi is not None:
            masks = crop_segmentations(self.migrator, roi)
            np.savez(file_out, **masks)
        else:
            shutil.copy(self.migrator.p, file_out)


class BeadMask_Migrator(DataMigrator, sd.BeadMask):
    IS_OPTIONAL = True
    COPY_CURRENT = False

    class BeadMask_V1(sd.BeadMask):
        DEFAULT_TARGET_PATH = 'protein/bead_mask.npz'

    def _migrate_method(self, out_path: str, roi=None, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if roi is not None:
            masks = crop_bead_mask(self.migrator, roi)
            np.savez(file_out, **masks)
        else:
            shutil.copy(self.migrator.p, file_out)


class Manifest_Migrator(DataMigrator, sd.Manifest):
    SCHEMA = {
        'probe': pl.String,
        'probe_id': pl.String,
        'gene_name': pl.String,
        'panel_type': pl.String,
        'probe_type': pl.String,
        'read_num': pl.Int32,
    }

    class Manifest_V1(sd.Manifest):
        DEFAULT_TARGET_PATH = 'transcript_panel.csv'

        SCHEMA = {
            'target_condensed': pl.String,
            'panel_type': pl.String,
        }

        schema_rename = {'target_condensed': 'gene_name'}

        def convert(self):
            df = self.load().rename(self.schema_rename)
            return df

    class Manifest_V2(sd.Manifest):
        DEFAULT_TARGET_PATH = 'transcript_panel.csv'

        SCHEMA = {
            'probe_name': pl.String,
            'gene_name': pl.String,
            'panel_type': pl.String,
            # 'probe_id': pl.String, # this is the only columns thats sometimes missing
        }

        schema_rename = {'probe_name': 'probe'}

        def convert(self):
            df = self.load().rename(self.schema_rename)
            return df

    def populate_missing_columns(self, df: pl.DataFrame) -> pl.DataFrame:
        has_probe_id = 'probe_id' in df.columns
        has_probe_type = 'probe_type' in df.columns

        dummy_value = '<unknown>'
        for col, dtype in self.SCHEMA.items():
            if col not in df.columns:
                dummy_value = dummy_value if dtype == pl.String else None
                df = df.with_columns(pl.lit(dummy_value).cast(dtype).alias(col))

        if not has_probe_id:
            df = (
                df.with_columns(gene_probe_idx=pl.int_range(0, pl.len()).over('gene_name'))
                .with_columns(pl.col('gene_probe_idx').cast(pl.Utf8).str.zfill(4) + '-' + pl.col('gene_name'))
                .drop('probe_id')
                .rename({'gene_probe_idx': 'probe_id'})
            )

        if not has_probe_type:
            df = df.with_columns(
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

        return self.order_output_columns(df)

    def order_output_columns(self, df):
        first_cols = [col for col in self.SCHEMA.keys() if col in df.columns]
        other_cols = [col for col in df.columns if col not in first_cols]
        return df.select(first_cols + other_cols)

    def _migrate_method(self, out_path: str, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)
        df = self.migrator.convert()
        df = self.populate_missing_columns(df)
        df.write_csv(file_out)


class RawFeatures_Migrator(DataMigrator, sd.RawFeatures):
    COPY_CURRENT = False

    class DummyFallback_V0(sd.RawFeatures):
        DEFAULT_TARGET_PATH = 'rna/transcript_table.csv.gz'

        SCHEMA = {
            'x_pixel_coordinate': pl.Float64,
            'y_pixel_coordinate': pl.Float64,
            'z_level': pl.Int64,
            'gene_name': pl.Boolean,
            'confidence_score': pl.String,
            'cell_id': pl.String,
        }

        def convert(self):
            return pl.LazyFrame(
                {
                    'TXUID': ['<unknown>'],
                    'sequence': ['<unknown>'],
                    'confidence_score': [None],
                    'y_pixel_coordinate': [None],
                    'x_pixel_coordinate': [None],
                    'z_level': [None],
                },
                schema=sd.RawFeatures.SCHEMA,
            )

    class RawFeatures_V1(sd.RawFeatures):
        DEFAULT_TARGET_PATH = 'diagnostics/transcript_table.parquet'

        SCHEMA = {
            'x_coord_shift': pl.Float64,
            'y_coord_shift': pl.Float64,
            'z': pl.Int64,
            'demuxed': pl.Boolean,
            'transcript_condensed': pl.String,
            'meanQS': pl.Float64,
            'cell_id': pl.UInt32,
            'sequence_to_demux': pl.String,
            'transcript': pl.String,
            'TXUID': pl.String,
        }

        rename_and_flip = {
            'sequence_to_demux': 'sequence',
            'meanQS': 'confidence_score',
            'x_coord_shift': 'y_pixel_coordinate',
            'y_coord_shift': 'x_pixel_coordinate',
            'z': 'z_level',
        }

        def convert(self):
            return self.load(lazy=True).rename(self.rename_and_flip)

    class RawFeatures_V2(sd.RawFeatures):
        DEFAULT_TARGET_PATH = 'rna/transcript_table.parquet'

        SCHEMA = {
            'y_pixel_coordinate': pl.Float64,
            'x_pixel_coordinate': pl.Float64,
            'z_level': pl.Int64,
            'demuxed': pl.Boolean,
            'probe_name': pl.String,
            'gene_name': pl.String,
            'confidence_score': pl.Float64,
            'in_nucleus': pl.UInt16,
            'cell_id': pl.UInt16,
            'sequence': pl.String,
            'TXUID': pl.String,
        }

        def convert(self):
            return self.load(lazy=True)

    def _migrate_method(self, out_path: str, roi=None, **kwargs):
        log = kwargs.get('logger', LOGGER)
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if 'current' in self.valid_versions:
            df = self.load(lazy=True)
        else:
            df = self.migrator.convert()

        if self.valid_versions == ['DummyFallback_V0']:
            log.warning(
                'Only fallback version available. Creating dummy output file to pass validators. Demuxing will be unavailable for this sample.'
            )

        elif roi is not None:
            df = crop_tx_features(df, roi)

        df.select(self.SCHEMA.keys()).sink_parquet(file_out)


class TxTable_Migrator(DataMigrator, sd.TxTable):
    COPY_CURRENT = False

    class TxTable_V1(sd.TxTable):
        SCHEMA = {
            'y_pixel_coordinate': pl.Float64,
            'x_pixel_coordinate': pl.Float64,
            'z_level': pl.Int64,
            'gene_name': pl.Boolean,
            'confidence_score': pl.String,
        }

        flipped_coord_order = ['x_pixel_coordinate', 'y_pixel_coordinate']

        col_rename = {
            'gene_name': c.GENE_ID_NAME,
        }

        flip_coords = {
            'x_pixel_coordinate': 'y_pixel_coordinate',
            'y_pixel_coordinate': 'x_pixel_coordinate',
        }

        drop_cols = ['cell_id', 'in_nucleus']

        def convert(self):
            lf = self.load(lazy=True)
            schema = lf.collect_schema().names()

            if 'TXUID' not in schema:
                lf = lf.with_row_index(name='TXUID')

            for col in self.drop_cols:
                if col in schema:
                    lf = lf.drop(col)

            coord_order = schema[0:2]
            if coord_order == self.flipped_coord_order:
                # print(f'Flipping xy-coordinates for {self._name}')
                self.col_rename.update(self.flip_coords)

            lf = lf.rename(self.col_rename)

            return lf

    def _migrate_method(self, out_path: str, roi=None, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if 'current' in self.valid_versions:
            df = self.load(lazy=True)
        else:
            df = self.migrator.convert()

        if roi is not None:
            df = crop_tx_features(df, roi)

        df.select(self.SCHEMA.keys()).sink_csv(file_out, compression='gzip')


class HnEDir_Migrator(DataMigrator, sd.HnEDir):
    COPY_CURRENT = False

    class HnEDir_V1(sd.HnEDir):
        DEFAULT_TARGET_PATH = 'h_and_e'
        EXPECTED_DIRS = {}

        FILE_MAP = {
            c.NUCLEAR_STAIN: ['nuclear.jp2'],
            c.CYTOPLASMIC_STAIN: ['cytoplasmic.jp2', 'eosin.jp2'],
            c.H_AND_E: ['h_and_e.jp2'],
        }

    def _migrate_method(self, out_path: str, roi=None, **kwargs):
        log = kwargs.get('logger', LOGGER)
        out_dir = io.pathval.ensure_dir(out_path / self.DEFAULT_TARGET_PATH)

        for img in self.migrator.mapped_files.keys():
            log.debug(f'Converting {img} to OME-TIFF format.')
            img_type = 'rgb' if img == c.H_AND_E else 'grey'
            migrate_image(self.migrator, img, out_dir, img_type=img_type, roi=roi)


class Protein_Migrator(DataMigrator, sd.ProteinDir):
    COPY_CURRENT = False

    class ProteinDir_V1(sd.ProteinDir):
        # this class handles the rare edge case where the protein panel is missing
        # but the individual protein images are present.
        # In this case, we can create a dummy panel and migrate the images as usual.

        EXPECTED_DIRS = {}

        @property
        def panel(self):
            return sd.ProteinPanel(root=self.root, validate_absence=True)

        @sd.validation_test
        def images_match_panel(self):
            return True

        def infer_panel(self) -> pl.DataFrame:
            proteins = sorted([img for img in self.mapped_files.keys()])
            df = pl.DataFrame({'target': proteins}).with_columns(pl.lit('standard').alias('panel_type'))
            return df

    class ProteinDir_V2(sd.ProteinDir):
        EXPECTED_DIRS = {}

    def _migrate_method(self, out_path: str, n_images: int = None, roi=None, **kwargs):
        log = kwargs.get('logger', LOGGER)
        out_dir = io.pathval.ensure_dir(out_path / self.DEFAULT_TARGET_PATH)

        migrator = self.migrator

        img_list = list(migrator.mapped_files)
        keep_list = img_list if n_images is None else img_list[0:n_images]
        keep_signals = [img for img in keep_list]
        total = len(keep_list)

        for protein_img in keep_list:
            log.debug(f'Converting {protein_img} ({keep_signals.index(protein_img) + 1}/{total})')
            migrate_image(migrator, protein_img, out_dir, img_type='grey', roi=roi)

        if self.valid_versions == ['ProteinDir_V1']:
            panel_df = migrator.infer_panel()
        else:
            panel_df = migrator.panel.load()

        panel_df = panel_df.filter(pl.col('target').is_in(keep_signals))
        panel_df.write_csv(out_path / self.migrator.panel.DEFAULT_TARGET_PATH)


class Metrics_Migrator(DataMigrator, DirectoryValidator):
    DEFAULT_TARGET_PATH = 'metrics'
    IS_OPTIONAL = True
    COPY_CURRENT = True


def migrate_image(migrator, img_name, out_dir, img_type='auto', roi=None):
    io.convert.jp2_to_ometiff(
        in_file=migrator.mapped_files[img_name],
        out_file=f'{out_dir}/{img_name}.ome.tiff',
        img_type=img_type,
        create_thumb=True,
        report_size=False,
        extent=roi.extent_array if roi else None,
    )


def crop_tx_features(df: pl.LazyFrame, roi):
    df = df.filter(
        pl.col('x_pixel_coordinate').is_between(*roi.xlims, closed='left'),
        pl.col('y_pixel_coordinate').is_between(*roi.ylims, closed='left'),
    ).with_columns(
        pl.col('x_pixel_coordinate') - roi.xlims[0],
        pl.col('y_pixel_coordinate') - roi.ylims[0],
    )
    return df


def crop_bead_mask(bead_mask, roi):

    cropped = {}
    arr = bead_mask.load()
    cropped[bead_mask.DEFAULT_KEY] = roi.crop_array(arr)
    return cropped


def crop_segmentations(segmentation, roi):
    def remove_boundary_labels(labels, inplace=False):
        arr = labels if inplace else labels.copy()

        boundary_labels = np.unique(
            np.concatenate(
                [  # Collect all labels touching the boundary
                    arr[0, :],  # top row
                    arr[-1, :],  # bottom row
                    arr[:, 0],  # left column
                    arr[:, -1],  # right column
                ]
            )
        )
        # Create a mask for boundary labels and set them to 0
        mask = np.isin(arr, boundary_labels)
        arr[mask] = 0

        return arr

    cleaned_masks = {}

    if segmentation.main_key == 'nuclei_exp':
        segmentation.main_key = 'nuclei'

    main_seg = segmentation.load()
    main_crop = roi.crop_array(main_seg)
    main_crop_cleaned = remove_boundary_labels(main_crop)

    keep_labels = np.unique(main_crop_cleaned)
    if len(keep_labels) == 1 and keep_labels[0] == 0:
        raise ValueError('No valid segmentation labels remain after cropping.')

    cleaned_masks[segmentation.main_key] = main_crop_cleaned

    for key in segmentation.available_keys:
        if key == segmentation.main_key:
            continue

        sub_seg = segmentation.load(key=key)
        sub_crop = roi.crop_array(sub_seg)
        mask = ~np.isin(sub_crop, keep_labels)
        sub_crop_cleaned = sub_crop.copy()
        sub_crop_cleaned[mask] = 0

        sub_labels = np.unique(sub_crop_cleaned)
        if not np.isin(keep_labels, sub_labels).all():
            raise ValueError(
                f'Sub-segmentation "{key}" contains labels not present in main segmentation after cropping.'
            )
        cleaned_masks[key] = sub_crop_cleaned

    return cleaned_masks
