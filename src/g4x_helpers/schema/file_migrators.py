import shutil

import numpy as np
import polars as pl

from .. import c, io
from . import definition


class SampleG4X_Migrator(definition.SampleG4X):
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
        return definition.SampleSheet(root=self.root)

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

        smp_id = pl.read_csv(self.root / option)['sample_id'].item()
        return smp_id

    @property
    def is_migratable(self):
        if self.is_valid and self.smp_sheet.is_valid and self.legacy_smp_id is not None:
            return True
        return False

    def migrate(self, out_path: str):
        out_path = io.pathval.validate_dir_path(out_path)
        _ = io.create_sample_g4x(
            sample_id=self.legacy_smp_id,
            run_meta=self.load(),
            ssheet=self.smp_sheet.p,
            out_path=out_path / definition.SampleG4X.DEFAULT_TARGET_PATH,
        )


class HnE_Migrator(definition.FolderValidator):
    DEFAULT_TARGET_PATH = c.HE_DIR

    FILE_MAP = {
        c.NUCLEAR_STAIN: ['nuclear.jp2'],
        c.CYTOPLASMIC_STAIN: ['cytoplasmic.jp2', 'eosin.jp2'],
        c.H_AND_E: ['h_and_e.jp2'],
    }

    @definition.validation_test
    def imgs_complete(self):
        existing_imgs = self.existing_imgs
        missing_imgs = set(self.FILE_MAP) - set(existing_imgs)
        if missing_imgs:
            print(f'Missing images: {missing_imgs}')
            return False
        return True

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

    def migrate(self, out_path: str, roi=None):
        out_path = io.pathval.validate_dir_path(out_path)
        out_dir = io.pathval.ensure_dir(out_path / self.DEFAULT_TARGET_PATH)

        is_all_jp2 = all([f.suffix == '.jp2' for f in self.existing_imgs.values()])
        if not is_all_jp2:
            print('Not all images are in JP2 format. Conversion not possible.')
            return

        for img, in_file in self.existing_imgs.items():
            img_type = 'rgb' if img == c.H_AND_E else 'grey'
            io.convert.jp2_to_ometiff(
                in_file=in_file,
                out_file=f'{out_dir}/{img}.ome.tiff',
                img_type=img_type,
                create_thumb=True,
                report_size=False,
                extent=roi.extent_array if roi else None,
            )


class Protein_Migrator(definition.ProteinDir):
    EXPECTED_DIRS = {}

    def migrate(self, out_path: str, roi=None):
        out_path = io.pathval.validate_dir_path(out_path)

        out_dir = io.pathval.ensure_dir(out_path / self.DEFAULT_TARGET_PATH)

        is_all_jp2 = all([f.suffix == '.jp2' for f in self.existing_images])
        if not is_all_jp2:
            print('Not all images are in JP2 format. Conversion not possible.')
            return

        for protein_img in self.existing_images:
            io.convert.jp2_to_ometiff(
                in_file=protein_img,
                out_file=f'{out_dir}/{protein_img.stem}.ome.tiff',
                img_type='grey',
                create_thumb=True,
                report_size=False,
                extent=roi.extent_array if roi else None,
            )

        shutil.copyfile(self.panel.p, out_path / self.panel.DEFAULT_TARGET_PATH)


class RawFeatures_Migrator(definition.BaseValidator):
    class RawFeatures_v1(definition.TableValidator):
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
            return self.load(lazy=True).rename(self.parquet_rename).select(definition.RawFeatures.SCHEMA.keys())

    class RawFeatures_v2(definition.TableValidator):
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
            return self.load(lazy=True).select(definition.RawFeatures.SCHEMA.keys())

    def versions(self):
        return {'v1': self.RawFeatures_v1(root=self.root), 'v2': self.RawFeatures_v2(root=self.root)}

    @property
    def valid_versions(self):
        return list({k for k, v in self.versions().items() if v.is_valid})

    @property
    def is_migratable(self):
        return len(self.valid_versions) > 0

    def migrate(self, out_path: str, roi=None):
        out_path = io.pathval.validate_dir_path(out_path)

        legacy_format = self.versions()[self.valid_versions[-1]]
        df = legacy_format.convert()

        if roi is not None:
            df = df.filter(
                pl.col('x_pixel_coordinate').is_between(*roi.xlims, closed='left'),
                pl.col('y_pixel_coordinate').is_between(*roi.ylims, closed='left'),
            ).with_columns(
                pl.col('x_pixel_coordinate') - roi.xlims[0],
                pl.col('y_pixel_coordinate') - roi.ylims[0],
            )

        file_out = io.pathval.ensure_parent_dir(out_path / definition.RawFeatures.DEFAULT_TARGET_PATH)
        df.sink_parquet(file_out)


# TODO collect all forms of manifests
class Manifest_Migrator(definition.TableValidator):
    DEFAULT_TARGET_PATH = 'transcript_panel.csv'

    SCHEMA = {
        'probe_name': pl.String,
        'gene_name': pl.String,
        'probe_id': pl.String,
        'panel_type': pl.String,
    }

    schema_rename = {'probe_name': 'probe'}

    def convert(self):
        return self.load().rename(self.schema_rename)

    def migrate(self, out_path: str):
        out_path = io.pathval.validate_dir_path(out_path)
        file_out = io.pathval.ensure_parent_dir(out_path / definition.Manifest.DEFAULT_TARGET_PATH)
        self.convert().write_csv(file_out)


class Segmentation_Migrator(definition.Segmentation):
    DEFAULT_TARGET_PATH = 'segmentation/segmentation_mask.npz'

    def migrate(self, out_path: str, roi=None):
        out_path = io.pathval.validate_dir_path(out_path)
        file_out = io.pathval.ensure_parent_dir(out_path / definition.Segmentation.DEFAULT_TARGET_PATH)

        if roi is not None:
            masks = crop_segmentations(self, roi)
            np.savez(file_out, **masks)
        else:
            shutil.copy(self.p, file_out)


class SampleSheet_Migrator(definition.SampleSheet):
    def migrate(self, out_path: str):
        out_path = io.pathval.validate_dir_path(out_path)
        file_out = io.pathval.ensure_parent_dir(out_path / definition.SampleSheet.DEFAULT_TARGET_PATH)
        shutil.copy(self.p, file_out)


class BeadMask_Migrator(definition.BeadMask):
    DEFAULT_TARGET_PATH = 'protein/bead_mask.npz'

    def migrate(self, out_path: str, roi=None):
        out_path = io.pathval.validate_dir_path(out_path)
        file_out = io.pathval.ensure_parent_dir(out_path / definition.BeadMask.DEFAULT_TARGET_PATH)

        if roi is not None:
            masks = crop_bead_mask(self, roi)
            np.savez(file_out, **masks)
        else:
            shutil.copy(self.p, file_out)


def crop_bead_mask(bead_mask, roi):
    cropped = {}
    arr = bead_mask.load()
    cropped[bead_mask.DEFAULT_KEY] = roi.crop_array(arr)
    return cropped


def crop_segmentations(segmentation, roi):
    cleaned_masks = {}

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
