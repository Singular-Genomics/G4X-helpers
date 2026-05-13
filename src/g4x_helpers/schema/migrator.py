import json
import re
import shutil

import numpy as np
import polars as pl

from .. import c, io
from . import definition as sd
from .validator import BaseValidator, FileValidator, FolderValidator


class DataMigrator(BaseValidator):
    IS_OPTIONAL = False
    VERSION_VALIDATORS = {}
    VERSION_CLASS_PATTERN = re.compile(r'.*_V\d+$')

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        auto_validators = {}
        for name, obj in cls.__dict__.items():
            if isinstance(obj, type) and issubclass(obj, BaseValidator) and cls.VERSION_CLASS_PATTERN.match(name):
                key = name[0].lower() + name[1:]

                auto_validators[key] = obj

        explicit_validators = getattr(cls, 'VERSION_VALIDATORS', {})

        cls.VERSION_VALIDATORS = {
            **auto_validators,
            **explicit_validators,
        }

    def versions(self):
        versions = {'current': self._current_validator()}
        versions.update(
            {name: validator_cls(root=self.root) for name, validator_cls in type(self).VERSION_VALIDATORS.items()}
        )
        return versions

    def _current_validator(self):
        for cls in type(self).__mro__:
            if cls is DataMigrator:
                continue

            if issubclass(cls, BaseValidator):
                return cls(root=self.root)

        raise TypeError(f'Could not infer "current" validator for {type(self).__name__}')

    @property
    def target_validator(self):
        return type(self).__bases__[1]

    @property
    def valid_versions(self):
        return list({k for k, v in self.versions().items() if v.is_valid})

    @property
    def mig_version(self):
        if len(self.valid_versions) == 0:
            raise ValueError(f'No valid versions found for {type(self).__name__} in {self.root}')
        if 'current' in self.valid_versions:
            mig_version = 'current'
        else:
            mig_version = self.valid_versions[-1]
        return self.versions()[mig_version]

    @property
    def is_migratable(self):
        if self.IS_OPTIONAL:
            return True, ''

        has_legacy = len(self.valid_versions) > 0
        err = 'No valid legacy versions are available.' if not has_legacy else ''

        return has_legacy, err

    @property
    def is_file(self):
        return isinstance(self, FileValidator)

    @property
    def is_folder(self):
        return isinstance(self, FolderValidator)

    def copy_if_current(self, out_path):
        if 'current' in self.valid_versions:
            print(f"Version 'current' is valid for {self.__class__.__name__}, copying file without migration.")
            if self.is_file:
                shutil.copy(self.p, out_path / self.DEFAULT_TARGET_PATH)
            elif self.is_folder:
                shutil.copytree(self.p, out_path / self.DEFAULT_TARGET_PATH)
            else:
                raise TypeError(f'Cannot copy {type(self).__name__} because it is neither a file nor a folder.')

    def validate_migration(self, out_path):
        target = self.target_validator(root=out_path)
        if target.is_valid:
            print(f'Migration successful for {type(self).__name__}!')
        else:
            raise ImportError(f'Migrated {type(self).__name__} did not pass validation! {target.validation()}.')

    def migrate(self, out_path, check_current: bool = True, *args, **kwargs):
        if not self.valid_versions and self.IS_OPTIONAL:
            print(f'No valid versions of {type(self).__name__} found. Skipping migration of optional data.')
            return

        test, msg = self.is_migratable
        if not test:
            raise ValueError(f'{type(self).__name__} is not migratable. {msg}')

        out_path = io.pathval.validate_dir_path(out_path)

        if check_current:
            self.copy_if_current(out_path)

        migrate_method = getattr(self, '_migrate_method', False)
        if not migrate_method:
            raise NotImplementedError(f'{type(self).__name__} does not have a _migrate_method defined.')
        else:
            try:
                result = migrate_method(out_path, *args, **kwargs)
                self.validate_migration(out_path)
                return result

            except Exception as e:
                raise ImportError(
                    f'Error migrating {type(self).__name__} from {self.target_rel}!\nreason: {e}'
                ) from None


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
    def is_migratable(self):

        if not self.smp_sheet.is_valid:
            return False, 'SampleSheet is not valid.'

        if self.legacy_smp_id is None:
            return False, 'Could not determine legacy sample_id.'

        if not len(self.valid_versions) > 0:
            return False, 'No valid legacy versions are available.'

        return True, ''

    def _migrate_method(self, out_path, **kwargs):

        with open(self.mig_version.p, 'r') as f:
            run_meta = json.load(f)

        if run_meta['protein_panel'] == []:
            run_meta['protein_panel'] = None

        _ = io.create_sample_g4x(
            sample_id=self.legacy_smp_id,
            run_meta=run_meta,
            ssheet=self.smp_sheet.p,
            out_path=out_path / self.DEFAULT_TARGET_PATH,
        )


class QCSummary_Migrator(DataMigrator, sd.QCSummary):
    IS_OPTIONAL = True

    def _migrate_method(self, out_path, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)
        shutil.copy(self.p, file_out)


class SampleSheet_Migrator(DataMigrator, sd.SampleSheet):
    def _migrate_method(self, out_path, **kwargs):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)
        shutil.copy(self.p, file_out)


class Segmentation_Migrator(DataMigrator, sd.Segmentation):
    class Segmentation_V1(sd.Segmentation):
        DEFAULT_TARGET_PATH = 'segmentation/segmentation_mask.npz'

    def crop_segmentations(self, roi):
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

        segmentation = self.mig_version
        main_seg = segmentation.load()
        main_crop = roi.crop_array(main_seg)
        main_crop_cleaned = remove_boundary_labels(main_crop)

        keep_labels = np.unique(main_crop_cleaned)

        cleaned_masks[segmentation.main_key] = main_crop_cleaned

        for key in segmentation.available_keys:
            if key == segmentation.main_key:
                continue

            sub_seg = segmentation.load(key=key)
            sub_crop = roi.crop_array(sub_seg)
            mask = ~np.isin(sub_crop, keep_labels)
            sub_crop_cleaned = sub_crop.copy()
            sub_crop_cleaned[mask] = 0
            cleaned_masks[key] = sub_crop_cleaned

        return cleaned_masks

    def _migrate_method(self, out_path: str, roi=None):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if roi is not None:
            masks = self.crop_segmentations(roi)
            np.savez(file_out, **masks)
        else:
            shutil.copy(self.mig_version.p, file_out)


class BeadMask_Migrator(DataMigrator, sd.BeadMask):
    IS_OPTIONAL = True

    class BeadMask_V1(sd.BeadMask):
        DEFAULT_TARGET_PATH = 'protein/bead_mask.npz'

    def crop_bead_mask(self, roi):
        bead_mask = self.mig_version

        cropped = {}
        arr = bead_mask.load()
        cropped[bead_mask.DEFAULT_KEY] = roi.crop_array(arr)
        return cropped

    def _migrate_method(self, out_path: str, roi=None):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if roi is not None:
            masks = self.crop_bead_mask(roi)
            np.savez(file_out, **masks)
        else:
            shutil.copy(self.mig_version.p, file_out)


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
        df = self.mig_version.convert()
        df = self.populate_missing_columns(df)
        df.write_csv(file_out)


class RawFeatures_Migrator(DataMigrator, sd.RawFeatures):
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

        parquet_rename = {
            'sequence_to_demux': 'sequence',
            'meanQS': 'confidence_score',
            'x_coord_shift': 'y_pixel_coordinate',
            'y_coord_shift': 'x_pixel_coordinate',
            'z': 'z_level',
        }

        def convert(self):
            return self.load(lazy=True).rename(self.parquet_rename)

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

    def _migrate_method(self, out_path: str, roi=None):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if self.valid_versions == ['dummyFallback_V0']:
            print(
                'Only fallback version available. Creating dummy output file to pass validators, but demuxing will be unavailable for this sample.'
            )
            data = {k: '<unknown>' for k in self.SCHEMA.keys()}
            df = pl.LazyFrame(data)

        else:
            if 'current' in self.valid_versions:
                df = self.load(lazy=True)
            else:
                df = self.mig_version.convert()

            if roi is not None:
                df = df.filter(
                    pl.col('x_pixel_coordinate').is_between(*roi.xlims, closed='left'),
                    pl.col('y_pixel_coordinate').is_between(*roi.ylims, closed='left'),
                ).with_columns(
                    pl.col('x_pixel_coordinate') - roi.xlims[0],
                    pl.col('y_pixel_coordinate') - roi.ylims[0],
                )

        df.sink_parquet(file_out)


class TxTable_Migrator(DataMigrator, sd.TxTable):
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
                print(f'Flipping xy-coordinates for {self.__class__.__name__}')
                self.col_rename.update(self.flip_coords)

            lf = lf.rename(self.col_rename)

            return lf

    def _migrate_method(self, out_path: str, roi=None):
        file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)

        if 'current' in self.valid_versions:
            df = self.load(lazy=True)
        else:
            df = self.mig_version.convert()

        if roi is not None:
            df = df.filter(
                pl.col('x_pixel_coordinate').is_between(*roi.xlims, closed='left'),
                pl.col('y_pixel_coordinate').is_between(*roi.ylims, closed='left'),
            ).with_columns(
                pl.col('x_pixel_coordinate') - roi.xlims[0],
                pl.col('y_pixel_coordinate') - roi.ylims[0],
            )

        df.sink_csv(file_out, compression='gzip')


class HnE_Migrator(DataMigrator, sd.HnEDir):
    class HnEDir_V1(sd.FileValidator):
        DEFAULT_TARGET_PATH = 'h_and_e'

        FILE_MAP = {
            c.NUCLEAR_STAIN: ['nuclear.jp2'],
            c.CYTOPLASMIC_STAIN: ['cytoplasmic.jp2', 'eosin.jp2'],
            c.H_AND_E: ['h_and_e.jp2'],
        }

        @property
        def existing_imgs(self):
            existing_imgs = {}
            for file_name, file_paths in self.FILE_MAP.items():
                for f in file_paths:
                    query = self.p / f
                    if query.exists():
                        existing_imgs[file_name] = query
                        break
            return existing_imgs

        @sd.validation_test
        def imgs_complete(self):
            missing_imgs = set(self.FILE_MAP) - set(self.existing_imgs)
            if missing_imgs:
                print(f'Missing images: {missing_imgs}')
                return False
            return True

    def _migrate_method(self, out_path: str, roi=None):
        out_dir = io.pathval.ensure_dir(out_path / self.DEFAULT_TARGET_PATH)

        migrator = self.mig_version

        is_all_jp2 = all([f.suffix == '.jp2' for f in migrator.existing_imgs.values()])
        if not is_all_jp2:
            print('Not all images are in JP2 format. Conversion not possible.')
            return

        for img, in_file in migrator.existing_imgs.items():
            img_type = 'rgb' if img == c.H_AND_E else 'grey'
            io.convert.jp2_to_ometiff(
                in_file=in_file,
                out_file=f'{out_dir}/{img}.ome.tiff',
                img_type=img_type,
                create_thumb=True,
                report_size=False,
                extent=roi.extent_array if roi else None,
            )


class Protein_Migrator(DataMigrator, sd.ProteinDir):
    class ProteinDir_V0(sd.ProteinDir):
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

        def create_protein_panel(self: str) -> pl.DataFrame:
            img_list = self.existing_images
            proteins = sorted([img.stem for img in img_list])
            df = pl.DataFrame({'target': proteins}).with_columns(pl.lit('standard').alias('panel_type'))
            return df

    class ProteinDir_V1(sd.ProteinDir):
        EXPECTED_DIRS = {}

    def _migrate_method(self, out_path: str, roi=None):
        out_dir = io.pathval.ensure_dir(out_path / self.DEFAULT_TARGET_PATH)

        migrator = self.mig_version

        is_all_jp2 = all([f.suffix == '.jp2' for f in migrator.existing_images])
        if not is_all_jp2:
            print('Not all images are in JP2 format. Conversion not possible.')
            return

        for protein_img in migrator.existing_images:
            io.convert.jp2_to_ometiff(
                in_file=protein_img,
                out_file=f'{out_dir}/{protein_img.stem}.ome.tiff',
                img_type='grey',
                create_thumb=True,
                report_size=False,
                extent=roi.extent_array if roi else None,
            )

        if self.valid_versions == ['proteinDir_V0']:
            panel_df = migrator.create_protein_panel()
            panel_df.write_csv(out_path / migrator.panel.DEFAULT_TARGET_PATH)
        else:
            shutil.copyfile(migrator.panel.p, out_path / migrator.panel.DEFAULT_TARGET_PATH)
