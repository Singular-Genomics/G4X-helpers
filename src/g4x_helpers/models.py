import json
import logging
import os
from dataclasses import InitVar, dataclass  # , asdict, field
from pathlib import Path
from typing import Optional, Union  # , Tuple, Any, Generator, Iterable, Iterator, List, Literal

import anndata
import glymur
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
import tifffile
from geopandas.geodataframe import GeoDataFrame
from packaging import version

import g4x_helpers.plotting as g4pl
import g4x_helpers.segmentation as reseg
import g4x_helpers.utils as utils

glymur.set_option('lib.num_threads', 16)


@dataclass()
class G4Xoutput:
    run_base: Union[Path, str]
    sample_id: str = None
    log_level: InitVar[int] = logging.INFO

    def __post_init__(self, log_level):
        self.run_base = Path(self.run_base)

        if self.sample_id is None:
            self.sample_id = self.run_base.name

        ## this only sets up a stream handler, no file handlers
        _ = self.setup_logger(stream_level=log_level, clear_handlers=True)

        with open(self.run_base / 'run_meta.json', 'r') as f:
            self.run_meta = json.load(f)

        self._populate_attrs()
        self._get_shape()
        self._get_cell_coords()
        self._get_tissue_crop()

        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / 'transcript_panel.csv', index_col=0, header=0)
            self.transcript_panel_dict = transcript_panel.to_dict()['panel_type']
            self.genes = list(self.transcript_panel_dict.keys())

        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / 'protein_panel.csv', index_col=0, header=0)
            self.protein_panel_dict = protein_panel.to_dict()['panel_type']
            self.proteins = list(self.protein_panel_dict.keys())

        roi_data = utils.Roi(extent=self.extent, name='data')
        roi_tissue = utils.Roi(extent=self.tissue_extent, name='tissue', sub_rois=2)
        self.rois = {
            'data': roi_data,
            'tissue': roi_tissue,
            'A': roi_tissue.sub_rois[0],
            'B': roi_tissue.sub_rois[1],
            'C': roi_tissue.scale(0.5),
            'D': roi_tissue.sub_rois[2],
            'E': roi_tissue.sub_rois[3],
        }

        self.cached_images = {}
        self.cached_thumbs = {}

    @property
    def logger(self) -> logging.Logger:
        return logging.getLogger(f'{self.sample_id}_G4XOutput')

    def _populate_attrs(self):
        for key in (
            'machine',
            'run_id',
            'fc',
            'lane',
            'platform',
            'software',
            'software_version',
            'transcript_panel',
            'protein_panel',
        ):
            setattr(self, key, self.run_meta.get(key))

        self.includes_transcript = True if self.transcript_panel is not None else False
        self.includes_protein = True if self.protein_panel is not None else False

    # region methods
    def setup_logger(
        self,
        *,
        stream_logger: Optional[bool] = True,
        stream_level: Optional[int] = logging.INFO,
        file_logger: Optional[bool] = False,
        file_level: Optional[int] = logging.INFO,
        clear_handlers: Optional[bool] = True,
    ) -> None:
        _ = utils.setup_logger(
            f'{self.sample_id}_G4XOutput',
            stream_logger=stream_logger,
            stream_level=stream_level,
            file_logger=file_logger,
            file_level=file_level,
            clear_handlers=clear_handlers,
        )

    def thumb(
        self, size=3, ax=None, crop='tissue', scale=1, scale_bar: bool = True, show_tissue_bounds=True, background=True
    ):
        if ax is None:
            fig, ax = plt.subplots(figsize=(size, size), dpi=300, layout='compressed')

        ax.set_title(self.sample_id)

        if not hasattr(self, 'he_thumb'):
            self.he_thumb = self._get_image(parent_directory='h_and_e', pattern='h_and_e_thumbnail.jpg')
        plot_image = self.he_thumb.copy()

        roi = self.rois['data'] if crop == 'data' else self.rois[crop]

        # if background:
        #     if isinstance(background, bool):
        #         background = '#F7F3F9'  # HE Image background color
        #     ax.set_facecolor(background)

        if scale != 1:
            roi = roi.scale(scale)

        display_size = np.max([roi.width / 5, roi.height / 5])
        plot_image = g4pl.downsample_mask(plot_image, ax=ax, mask_length=display_size, downsample='auto')

        ax.imshow(plot_image, extent=self.rois['data'].extent, origin='lower')

        if show_tissue_bounds:
            if crop == 'data':
                crop_roi = self.rois['tissue']
                crop_roi.add_to_plot(ax, color='orange', label=True, pad=-0.1, label_position='top_center')

        if scale_bar:
            um = roi.width * 0.3125
            rounded_scale = round((um / 10) / 50) * 100
            g4pl.scale_bar(
                ax, length=rounded_scale, unit='micron', theme='light', background=True, lw=1, background_alpha=0.75
            )

        ax.set_xlim(roi.xlims)
        ax.set_ylim(roi.ylims)
        ax.set_xticks([])
        ax.set_yticks([])

    def load_adata(
        self, *, remove_nontargeting: Optional[bool] = True, load_clustering: Optional[bool] = True
    ) -> anndata.AnnData:
        adata = sc.read_h5ad(self.run_base / 'single_cell_data' / 'feature_matrix.h5')

        # TODO this may conflict with the last line?!
        # adata.obs_names = adata.obs['cell_id']
        adata.var_names = adata.var['gene_id']

        adata.obs['sample_id'] = adata.uns['sample_id'] = self.sample_id
        adata.uns['software_version'] = self.software_version

        if remove_nontargeting:
            adata = adata[:, adata.var.query(" probe_type == 'targeting' ").index].copy()

        if load_clustering:
            try:
                df = pd.read_csv(self.run_base / 'single_cell_data' / 'clustering_umap.csv.gz', index_col=0, header=0)
                adata.obs = adata.obs.merge(df, how='left', left_index=True, right_index=True)
                umap_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
                adata.obsm['X_umap'] = adata.obs[[f'{umap_key}_1', f'{umap_key}_2']].to_numpy(dtype=float)
            except Exception as e:
                self.logger.exception('No single-cell clustering results found.')
                self.logger.exception(e)

        adata.obs_names = f'{self.sample_id}-' + adata.obs['cell_id'].str.split('-').str[1]
        return adata

    def load_protein_image(self, protein: str) -> np.ndarray:
        if not self.protein_panel:
            self.logger.error('No protein results.')
            return None
        assert protein in self.proteins, print(f'{protein} not in protein panel.')
        return self._get_image('protein', f'{protein}.*')

    def load_he_image(self) -> np.ndarray:
        return self._get_image('h_and_e', 'h_and_e.*')

    def load_nuclear_image(self) -> np.ndarray:
        return self._get_image('h_and_e', 'nuclear.*')

    def load_eosin_image(self) -> np.ndarray:
        return self._get_image('h_and_e', 'eosin.*')

    def load_segmentation(self, expanded: Optional[bool] = True) -> np.ndarray:
        arr = np.load(self.run_base / 'segmentation' / 'segmentation_mask.npz')
        if expanded:
            return arr['nuclei_exp']
        else:
            return arr['nuclei']

    def load_transcript_table(self, return_polars: Optional[bool] = False) -> Union[pd.DataFrame, pl.DataFrame]:
        reads = pl.read_csv(self.run_base / 'rna' / 'transcript_table.csv.gz')
        if not return_polars:
            reads = reads.to_pandas()
        return reads

    def list_contents(self, subdir=None):
        if subdir is None:
            subdir = ''

        list_path = self.run_base / subdir
        output = os.listdir(list_path)

        contents = {'dirs': [], 'files': []}
        for item in output:
            if os.path.isdir(list_path / item):
                contents['dirs'].append(item)
            if os.path.isfile(list_path / item):
                contents['files'].append(item)

        return contents

    def add_roi(self, roi_name: str, xlims: tuple, ylims: tuple):
        add_roi = utils.Roi(xlims=xlims, ylims=ylims, name=roi_name)
        self.rois[roi_name] = add_roi

    def run_g4x_segmentation(
        self, labels: GeoDataFrame | np.ndarray, out_dir=None, include_channels=None, exclude_channels=None
    ):
        mask = labels

        if isinstance(mask, GeoDataFrame):
            print('Rasterizing provided GeoDataFrame')
            mask = reseg.rasterize_polygons(gdf=mask, target_shape=self.shape)

        print('Extracting mask properties')
        segmentation_props = reseg.get_mask_properties(self, mask)

        print('Assigning transcripts to mask labels')
        reads_new_labels = reseg.assign_tx_to_mask_labels(self, mask)

        print('Extracting image signals')
        cell_by_protein = reseg.extract_image_signals(
            self, mask, include_channels=include_channels, exclude_channels=exclude_channels
        )

        print('Building output data structures')
        cell_metadata = reseg._make_cell_metadata(segmentation_props, cell_by_protein)
        cell_by_gene = reseg._make_cell_by_gene(segmentation_props, reads_new_labels)
        adata = reseg._make_adata(cell_by_gene, cell_metadata)

        if out_dir:
            print(f'Saving output files to {out_dir}')
        else:
            print(f'No output directory specified, saving to ["custom"] directories in {self.run_base}.')

        outfile = utils._create_custom_out(self, out_dir, 'segmentation', 'segmentation_mask.npz')
        np.savez(outfile, cell_labels=mask)

        outfile = utils._create_custom_out(self, out_dir, 'rna', 'transcript_table.csv.gz')
        reads_new_labels.write_csv(outfile)

        outfile = utils._create_custom_out(self, out_dir, 'single_cell_data', 'cell_by_transcript.csv.gz')
        cell_by_gene.write_csv(outfile)

        outfile = utils._create_custom_out(self, out_dir, 'single_cell_data', 'cell_by_protein.csv.gz')
        cell_by_protein.write_csv(outfile)

        outfile = utils._create_custom_out(self, out_dir, 'single_cell_data', 'feature_matrix.h5')
        adata.write_h5ad(outfile)

        outfile = utils._create_custom_out(self, out_dir, 'single_cell_data', 'cell_metadata.csv.gz')
        adata.obs.to_csv(outfile)

    def plot_signals(self, signals, size=3, crop='tissue', thumbnail=True, **kwargs):
        if not isinstance(signals, list):
            signals = [signals]

        axs = g4pl.montage_plot(len(signals), panel_size=size, n_cols=4, layout='compressed', hspace=0.1)

        if len(signals) == 1:
            axs = [axs]

        for ax, signal in zip(axs, signals):
            self.plot_signal(signal=signal, size=size, ax=ax, crop=crop, thumbnail=thumbnail, **kwargs)

    def plot_signal(
        self, signal=None, size=4.5, ax=None, crop='tissue', scale_bar=True, thumbnail=False, clip=True, scale=1
    ):
        if ax is None:
            fig, ax = plt.subplots(figsize=(size, size), layout='compressed')
        else:
            fig = ax.figure

        ax.set_title(signal)

        if not thumbnail and signal in self.cached_images:  # print(f'Using cached image for {signal}')
            plot_image, vmax, vmin = self.cached_images[signal].values()

        elif thumbnail and signal in self.cached_thumbs:  # print(f'Using cached thumb for {signal}')
            plot_image, vmax, vmin = self.cached_thumbs[signal].values()

        else:  # print(f'Loading image for {signal}')
            plot_image, vmax, vmin = self.load_image(signal, thumbnail=thumbnail).values()

        if not clip:
            vmax = vmin = None
        elif isinstance(clip, tuple):
            vmin, vmax = clip

        roi = self.rois['data'] if crop == 'data' else self.rois[crop]

        if scale != 1:
            roi = roi.scale(scale)

        fct = 5 if thumbnail else 1
        roi_size = np.max([roi.width / fct, roi.height / fct])

        plot_image = g4pl.downsample_mask(plot_image, mask_length=roi_size, ax=ax, downsample='auto')

        # print(f'Plotting {signal} with vmin={vmin}, vmax={vmax} and size {plot_image.shape}')
        im = ax.imshow(
            plot_image, vmin=vmin, vmax=vmax, cmap=g4pl.cm.lipari, extent=self.rois['data'].extent, origin='lower'
        )
        g4pl._add_colorbar(ax=ax, cax=im, loc='bottom')
        # fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.025, fraction=0.05, aspect=35)

        if scale_bar:
            um = roi.width * 0.3125
            rounded_scale = round((um / 10) / 50) * 100
            g4pl.scale_bar(
                ax, length=rounded_scale, unit='micron', theme='dark', background=True, lw=1, background_alpha=0.5
            )

        ax.set_xlim(roi.xlims)
        ax.set_ylim(roi.ylims)
        ax.set_xticks([])
        ax.set_yticks([])

    def load_image(self, signal=str, thumbnail=False):
        def _get_image(img_file):
            if img_file.suffix == '.png':
                return plt.imread(img_file)
            else:
                return glymur.Jp2k(img_file)[:]

        if signal in ('h_and_e', 'nuclear', 'eosin'):
            folder = 'h_and_e'
        else:
            folder = 'protein'

        img_folder = self.run_base / folder
        filename = f'{signal}_thumbnail.png' if thumbnail else f'{signal}.jp2'
        img_file = img_folder / filename

        if thumbnail and signal not in self.cached_thumbs:
            image = _get_image(img_file)

            vmax = np.quantile(image, 0.995)
            vmin = vmax * 0.05

            self.cached_thumbs[signal] = {'img': image, 'vmax': vmax, 'vmin': vmin}
            return self.cached_thumbs[signal]

        elif not thumbnail and signal not in self.cached_images:
            image = _get_image(img_file)

            vmax = np.quantile(image, 0.995)
            vmin = vmax * 0.05

            self.cached_images[signal] = {'img': image, 'vmax': vmax, 'vmin': vmin}
            return self.cached_images[signal]

        elif thumbnail and signal in self.cached_thumbs:
            return self.cached_thumbs[signal]

        elif not thumbnail and signal in self.cached_images:
            return self.cached_images[signal]

    # region internal
    def _get_image(self, parent_directory: str, pattern: str) -> np.ndarray:
        img_path = next((self.run_base / parent_directory).glob(pattern), None)
        if img_path is None:
            self.logger.error('Image file does not exist.')
            return None
        if img_path.suffix in ('.jp2', '.jpg'):
            img = glymur.Jp2k(img_path)[:]
        else:
            img = tifffile.imread(img_path)
        return img

    def _get_coord_order(self, verbose: bool = False) -> str:
        critical_version = version.parse('2.11.1')

        detected_caretta_version = version.parse(self.software_version)
        if verbose:
            print(f'Detected Caretta version: {detected_caretta_version}')

        if detected_caretta_version >= critical_version:
            return 'yx'
        else:
            return 'xy'

    def _get_shape(self):
        image_path = self.run_base / 'h_and_e' / 'nuclear.jp2'

        image = glymur.Jp2k(image_path)
        self.shape = image.shape
        self.extent = (0, self.shape[1], 0, self.shape[0])

    def _get_cell_coords(self):
        cell_meta_path = self.run_base / 'single_cell_data' / 'cell_metadata.csv.gz'

        if self._get_coord_order() == 'yx':
            coord_select = ['cell_y', 'cell_x']
        else:
            coord_select = ['cell_x', 'cell_y']

        self.cell_coords = pl.read_csv(cell_meta_path).select(coord_select).to_numpy()

    def _get_tissue_crop(self, ubins=100, threshold=0.1, expand=0.1, order='yx', **kwargs):
        def lims_to_extent(lims, lims_order='yx'):
            if lims_order == 'yx':
                y0, y1 = lims[0]
                x0, x1 = lims[1]
            elif lims_order == 'xy':
                x0, x1 = lims[0]
                y0, y1 = lims[1]

            return (x0, x1, y0, y1)

        lims = utils.find_tissue(
            self.cell_coords,
            ubins=ubins,
            threshold=threshold,
            expand=expand,
            order=order,
            smp_shape=self.shape,
            **kwargs,
        )
        self.tissue_lims = lims['lims']
        self.tissue_center = lims['centroid']
        self.tissue_extent = lims_to_extent(lims['lims'])

    def __repr__(self):
        machine_num = self.machine.removeprefix('g4-').lstrip('0')
        mac_run_id = f'G{machine_num.zfill(2)}-{self.run_id}'

        repr_string = f'G4X_output @ {self.run_base}\n'
        repr_string += f'Sample: \033[1m{self.sample_id}\033[0m of {mac_run_id}, {self.fc}\n\n'

        if self.includes_transcript:
            repr_string += f'Transcript panel with {len(self.genes)} genes\t[{", ".join(self.genes[0:5])} ... ]\n'

        if self.includes_protein:
            repr_string += f'Protein panel with {len(self.proteins)} proteins\t[{", ".join(self.proteins[0:5])} ... ]\n'

        return repr_string
