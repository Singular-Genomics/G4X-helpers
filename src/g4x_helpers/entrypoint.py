import argparse
import atexit
import json
import os
import shutil
import signal
import sys
import math
import tarfile
from pathlib import Path
from datetime import datetime
import geopandas
import numpy as np

from g4x_helpers.models import G4Xoutput
from g4x_helpers.utils import verbose_to_log_level, gzip_file, delete_existing

SUPPORTED_MASK_FILETYPES = ['.npz', '.npy', '.geojson']


def try_load_segmentation(
    segmentation_mask: Path, expected_shape: tuple[int], segmentation_mask_key: str | None = None
) -> np.ndarray | geopandas.GeoDataFrame:
    ## load new segmentation
    suffix = segmentation_mask.suffix.lower()
    if suffix == '.npz':
        with np.load(segmentation_mask) as labels:
            available_keys = list(labels.keys())

            if segmentation_mask_key:  # if a key is specified
                if segmentation_mask_key not in labels:
                    raise KeyError(f"Key '{segmentation_mask_key}' not found in .npz; available keys: {available_keys}")
                seg = labels[segmentation_mask_key]

            else:  # if no key specified,
                if len(labels) == 1:  # and only one key is available, use that key
                    seg = labels[available_keys[0]]
                else:  # and multiple keys are available, raise an error
                    if len(labels) > 1:
                        raise ValueError(
                            f'Multiple keys found in .npz: {available_keys}. Please specify a key using segmentation_mask_key.'
                        )
                seg = labels

    elif suffix == '.npy':
        # .npy: directly returns the array, no context manager available
        if segmentation_mask_key is not None:
            print('file is .npy, ignoring provided segmentation_mask_key.')
        seg = np.load(segmentation_mask)

    elif suffix == '.geojson':
        seg = geopandas.read_file(segmentation_mask)

        if segmentation_mask_key is not None:
            if segmentation_mask_key not in seg.columns:
                raise KeyError(
                    f"Column '{segmentation_mask_key}' not found in GeoJSON; available columns: {seg.columns.tolist()}"
                )

            # ensure that a coliumn named 'label' exists
            seg['label'] = seg[segmentation_mask_key]

        else:
            if 'label' not in seg.columns:
                raise ValueError(
                    "No column named 'label' found in GeoJSON. Please specify which column to use for labels via segmentation_mask_key."
                )

    # only validate shape for numpy arrays
    if isinstance(seg, np.ndarray):
        assert seg.shape == expected_shape, (
            f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}'
        )

    return seg


def launch_resegment():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument(
        '--segmentation_mask',
        help='Path to new segmentation mask. Supported files types are: .npy, .npz, .geojson.',
        action='store',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--sample_id', help='sample_id (Optional)', action='store', type=str, required=False, default=None
    )
    parser.add_argument(
        '--out_dir',
        help='Output directory where new files will be saved. Will be created if it does not exist. If not provided, the files in run_base will be updated in-place.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--segmentation_mask_key',
        help='Key in npz/geojson where mask/labels should be taken from (Optional)',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--threads',
        help='Number of threads to use for processing. [4]',
        action='store',
        type=int,
        required=False,
        default=4,
    )
    parser.add_argument(
        '--verbose',
        help='Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]',
        action='store',
        type=int,
        default=1,
    )

    args = parser.parse_args()

    ## preflight checks
    run_base = Path(args.run_base)
    assert run_base.exists(), f'{run_base} does not appear to exist.'
    segmentation_mask = Path(args.segmentation_mask)
    assert segmentation_mask.suffix in SUPPORTED_MASK_FILETYPES, (
        f'{segmentation_mask.suffix} not a supported file type.'
    )

    ## initialize G4X sample
    sample = G4Xoutput(
        run_base=run_base, sample_id=args.sample_id, out_dir=args.out_dir, log_level=verbose_to_log_level(args.verbose)
    )
    print(sample)

    ## load new segmentation
    labels = try_load_segmentation(segmentation_mask, sample.shape, args.segmentation_mask_key)

    ## run intersection with new segmentation
    _ = sample.intersect_segmentation(labels=labels, out_dir=args.out_dir, n_threads=args.threads)


def launch_update_bin():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument(
        '--bin_file', help='Path to G4X-Viewer segmentation bin file.', action='store', type=str, required=True
    )
    parser.add_argument('--out_path', help='Output file path', action='store', type=str, required=True)
    parser.add_argument(
        '--metadata',
        help='Path to metadata table where clustering and/or embedding information will be extracted. Table must contain a header.',
        action='store',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--cellid_key',
        help='Column name in metadata containing cell IDs that match with bin_file. If not provided, assumes that first column in metadata contains the cell IDs.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--cluster_key',
        help='Column name in metadata containing cluster IDs. Automatically assigns new colors if cluster_color_key is not provided. If not provided, skips updating cluster IDs.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--cluster_color_key',
        help='Column name in metadata containing cluster colors. Colors must be provided as hex codes. If provided, cluster_key must also be provided. If not provided, skips updating cluster colors.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--emb_key',
        help='Column name in metadata containing embedding. Parser will look for {emb_key}_1 and {emb_key}_2. If not provided, skips updating embedding.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--verbose',
        help='Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]',
        action='store',
        type=int,
        default=1,
    )

    args = parser.parse_args()

    ## preflight checks
    bin_file = Path(args.bin_file)
    if not bin_file.exists():
        raise FileNotFoundError(f'{bin_file} does not appear to exist.')
    metadata = Path(args.metadata)
    if not metadata.exists():
        raise FileNotFoundError(f'{metadata} does not appear to exist.')
    out_path = Path(args.out_path)
    out_dir = out_path.parent
    os.makedirs(out_dir, exist_ok=True)

    ## run converter
    from g4x_helpers.g4x_viewer.bin_generator import seg_updater

    _ = seg_updater(
        bin_file=bin_file,
        metadata_file=metadata,
        out_path=out_path,
        cellid_key=args.cellid_key,
        cluster_key=args.cluster_key,
        cluster_color_key=args.cluster_color_key,
        emb_key=args.emb_key,
        log_level=verbose_to_log_level(args.verbose),
    )


def launch_new_bin():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument(
        '--out_dir',
        help='Output directory where new files will be saved. Will be created if it does not exist. If not provided, the files in run_base will be updated in-place.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--threads',
        help='Number of threads to use for processing. [4]',
        action='store',
        type=int,
        required=False,
        default=4,
    )
    parser.add_argument(
        '--verbose',
        help='Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]',
        action='store',
        type=int,
        default=1,
    )

    args = parser.parse_args()

    ## preflight checks
    run_base = Path(args.run_base)
    assert run_base.exists(), f'{run_base} does not appear to exist.'

    ## initialize G4X sample
    sample = G4Xoutput(run_base=run_base, out_dir=args.out_dir, log_level=verbose_to_log_level(args.verbose))
    print(sample)

    ## set up the data
    try:
        adata = sample.load_adata()
        emb_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
        cluster_key = sorted([x for x in adata.obs.columns if 'leiden' in x])[0]
        sample.logger.info('Successfully loaded adata with clustering information.')
    except Exception:
        adata = sample.load_adata(load_clustering=False)
        emb_key = None
        cluster_key = None
        sample.logger.info('Clustering information was not found, cell coloring will be random.')

    mask = sample.load_segmentation()

    if args.out_dir is not None:
        os.makedirs(args.out_dir, exist_ok=True)
        outfile = Path(args.out_dir) / f'{sample.sample_id}.bin'
    else:
        outfile = run_base / 'g4x_viewer' / f'{sample.sample_id}.bin'

    ## make new bin file
    from g4x_helpers.g4x_viewer.bin_generator import seg_converter

    sample.logger.info('Making G4X-Viewer bin file.')
    _ = seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=outfile,
        cluster_key=cluster_key,
        emb_key=emb_key,
        protein_list=[f'{x}_intensity_mean' for x in sample.proteins],
        n_threads=args.threads,
        logger=sample.logger,
    )
    sample.logger.debug(f'G4X-Viewer bin --> {outfile}')


def launch_tar_viewer():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        '--viewer_dir', help='Path to G4X-viewer folder to tar', action='store', type=str, required=True
    )
    parser.add_argument('--out_path', help='Output file path', action='store', type=str, required=True)

    args = parser.parse_args()

    print('Checking files.')
    viewer_dir = Path(args.viewer_dir)
    assert viewer_dir.exists(), f'{viewer_dir} does not appear to exist.'

    bin_path = list(viewer_dir.glob('*.bin'))
    assert len(bin_path) == 1, 'Either no bin file was found in viewer_dir or multiple bin files were found.'
    bin_path = bin_path[0]

    sample_id = bin_path.stem

    ome_tiff_path = viewer_dir / f'{sample_id}.ome.tiff'
    assert ome_tiff_path.is_file(), 'ome.tiff file does not exist.'

    run_meta_path = viewer_dir / f'{sample_id}_run_metadata.json'
    assert run_meta_path.is_file(), 'run_metadata.json file does not exist.'

    tx_path = viewer_dir / f'{sample_id}.tar'
    assert tx_path.is_file(), 'transcript tar file does not exist.'

    # --- H&E paths
    h_and_e_path = viewer_dir / 'h_and_e'
    created_he_dir = False
    if not h_and_e_path.exists():
        h_and_e_path.mkdir(parents=True, exist_ok=True)
        created_he_dir = True

    orig_he = viewer_dir / f'{sample_id}_HE.ome.tiff'
    assert orig_he.is_file(), 'fH&E ome.tiff file does not exist.'
    moved_he = h_and_e_path / f'{sample_id}_HE.ome.tiff'

    moved = False
    restored = False

    def restore_he():
        nonlocal moved, restored, created_he_dir
        if restored:
            return
        try:
            if moved and moved_he.exists():
                if orig_he.exists():
                    orig_he.unlink()
                shutil.move(str(moved_he), str(orig_he))
                print('Restored H&E file to original location.')

            # NEW: remove the h_and_e folder if we created it and it’s empty
            try:
                if created_he_dir and h_and_e_path.exists() and not any(h_and_e_path.iterdir()):
                    h_and_e_path.rmdir()
                    print('Removed empty h_and_e directory.')
            except Exception as e:
                print(f'WARNING: failed to remove h_and_e directory: {e}', file=sys.stderr)

            restored = True
            moved = False
        except Exception as e:
            print(f'WARNING: failed to restore H&E file: {e}', file=sys.stderr)

    atexit.register(restore_he)

    def _handle_signal(sig):
        import signal as _signal

        sig_name = _signal.Signals(sig).name
        print(f'Received {sig_name}; cleaning up…', file=sys.stderr)
        restore_he()
        sys.exit(128 + sig)

    signal.signal(signal.SIGINT, _handle_signal)
    signal.signal(signal.SIGTERM, _handle_signal)

    # Move H&E into subfolder for packaging
    shutil.move(str(orig_he), str(moved_he))
    moved = True

    try:
        print('Making metadata.')
        metadata = {
            'protein_image_src': ome_tiff_path.name,
            'protein_image_data_src': run_meta_path.name,
            'he_images_src': h_and_e_path.name,
            'cell_segmentation_src': bin_path.name,
            'transcript_src': tx_path.name,
        }
        with open(viewer_dir / 'dataset.config.json', 'w') as f:
            json.dump(metadata, f)

        print('Tarring folder.')
        out_tar = Path(args.out_path)
        if not out_tar.exists():
            out_tar.mkdir(parents=True, exist_ok=True)
        out_tar = out_tar / f'{sample_id}_g4x_viewer.tar'
        with tarfile.open(out_tar, 'w') as tar:
            tar.add(viewer_dir, arcname=viewer_dir.name)

    finally:
        restore_he()
        try:
            atexit.unregister(restore_he)
        except Exception:
            pass


def launch_redemux():
    from g4x_helpers.demux import batched_dot_product_hamming_matrix, demux, load_manifest
    from g4x_helpers.g4x_viewer.tx_generator import tx_converter
    import polars as pl
    from tqdm import tqdm

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument(
        '--manifest',
        help='Path to manifest for demuxing. The manifest must be a 3-column csv with the following header: target,sequence,read.',
        action='store',
        type=str,
        required=True,
    )
    parser.add_argument(
        '--batch_size',
        help='Number of transcripts to process per batch.',
        action='store',
        type=int,
        required=False,
        default=1_000_000,
    )
    parser.add_argument(
        '--out_dir',
        help='Output directory where new files will be saved. Will be created if it does not exist. If not provided, the files in run_base will be updated in-place.',
        action='store',
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        '--threads',
        help='Number of threads to use for processing. [4]',
        action='store',
        type=int,
        required=False,
        default=4,
    )
    parser.add_argument(
        '--verbose',
        help='Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]',
        action='store',
        type=int,
        default=1,
    )

    args = parser.parse_args()

    ## preflight checks
    run_base = Path(args.run_base)
    assert run_base.exists(), f'{run_base} does not appear to exist.'
    manifest = Path(args.manifest)
    assert manifest.exists(), f'{manifest} does not appear to exist.'
    demux_batch_size = args.batch_size

    ## make output directory with symlinked files from original
    out_dir = Path(args.out_dir)
    batch_dir = out_dir / 'diagnostics' / 'batches'
    batch_dir.mkdir(parents=True, exist_ok=True)

    ignore_file_list = ['clustering_umap.csv.gz', 'dgex.csv.gz', 'transcript_panel.csv']
    for root, dirs, files in os.walk(run_base):
        rel_root = Path(root).relative_to(run_base)
        if str(rel_root) == 'metrics':
            continue
        dst_root = out_dir / rel_root
        dst_root.mkdir(parents=True, exist_ok=True)
        for f in files:
            if f in ignore_file_list:
                continue
            src_file = Path(root) / f
            dst_file = dst_root / f
            if dst_file.exists():
                dst_file.unlink()
            dst_file.symlink_to(src_file)

    ## initialize G4X sample
    sample = G4Xoutput(run_base=run_base, out_dir=out_dir, log_level=verbose_to_log_level(args.verbose))
    print(sample)

    ## update metadata and transcript panel file
    shutil.copy(manifest, out_dir / 'transcript_panel.csv')
    with open(out_dir / 'run_meta.json', 'r') as f:
        meta = json.load(f)
    meta['transcript_panel'] = manifest.name
    meta['redemuxed_timestamp'] = f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
    (out_dir / 'run_meta.json').unlink()
    with open(out_dir / 'run_meta.json', 'w') as f:
        _ = json.dump(meta, f)
    meta = {'run_metadata': meta}
    (out_dir / 'g4x_viewer' / f'{sample.sample_id}_run_metadata.json').unlink()
    with open(out_dir / 'g4x_viewer' / f'{sample.sample_id}_run_metadata.json', 'w') as f:
        _ = json.dump(meta, f)

    ## load the new manifest file that we will demux against
    sample.logger.info('Loading manifest file.')
    manifest, probe_dict = load_manifest(manifest)
    seq_reads = manifest['read'].unique().to_list()
    seq_reads = [int(x.split('_')[-1]) if isinstance(x, str) else x for x in seq_reads]

    ## do the re-demuxing
    sample.logger.info('Performing re-demuxing.')
    num_features = pl.scan_parquet(sample.feature_table_path).select(pl.len()).collect().item()
    num_expected_batches = math.ceil(num_features / args.batch_size)
    cols_to_select = [
        'x_coord_shift',
        'y_coord_shift',
        'z',
        'demuxed',
        'transcript_condensed',
        'meanQS',
        'sequence_to_demux',
        'transcript',
        'TXUID',
    ]
    for i, feature_batch in tqdm(
        enumerate(sample.stream_features(args.batch_size, cols_to_select)),
        total=num_expected_batches,
        desc='Demuxing transcripts',
        position=0,
    ):
        feature_batch = feature_batch.with_columns(pl.col('TXUID').str.split('_').list.last().cast(int).alias('read'))
        redemuxed_feature_batch = []
        for seq_read in seq_reads:
            feature_batch_read = feature_batch.filter(pl.col('read') == seq_read)
            manifest_read = manifest.filter(pl.col('read') == seq_read)
            if len(feature_batch_read) == 0 or len(manifest_read) == 0:
                continue
            seqs = feature_batch_read['sequence_to_demux'].to_list()
            codes = manifest_read['sequence'].to_list()
            codebook_target_ids = np.array(manifest_read['target'].to_list())
            hammings = batched_dot_product_hamming_matrix(seqs, codes, batch_size=demux_batch_size)
            feature_batch_read = demux(hammings, feature_batch_read, codebook_target_ids, probe_dict)
            redemuxed_feature_batch.append(feature_batch_read)
        pl.concat(redemuxed_feature_batch).write_parquet(batch_dir / f'batch_{i}.parquet')

    ## set run_base to the redemux output folder
    sample.run_base = out_dir

    ## concatenate results into final csv and parquet
    sample.logger.info('Writing updated transcript table.')
    final_tx_table_path = out_dir / 'rna' / 'transcript_table.csv'
    final_tx_pq_path = out_dir / 'diagnostics' / 'transcript_table.parquet'
    delete_existing(final_tx_table_path)
    delete_existing(final_tx_table_path.with_suffix('.csv.gz'))
    delete_existing(final_tx_pq_path)
    tx_table = pl.scan_parquet(list(batch_dir.glob('*.parquet')))
    tx_table.sink_parquet(final_tx_pq_path)
    _ = (
        tx_table.filter(pl.col('demuxed_new'))
        .drop(
            'transcript',
            'transcript_new',
            'transcript_condensed',
            'demuxed',
            'demuxed_new',
            'sequence_to_demux',
            'TXUID',
        )
        .rename(
            {
                'transcript_condensed_new': 'gene_name',
                'x_coord_shift': 'y_pixel_coordinate',
                'y_coord_shift': 'x_pixel_coordinate',
                'z': 'z_level',
                'meanQS': 'confidence_score',
            }
        )
        .sink_csv(final_tx_table_path)
    )
    _ = gzip_file(final_tx_table_path)
    shutil.rmtree(batch_dir)

    ## now regenerate the secondary files
    sample.logger.info('Regenerating downstream files.')
    labels = sample.load_segmentation()
    sample.logger.info('Intersecting with existing segmentation.')
    _ = sample.intersect_segmentation(labels=labels, out_dir=out_dir, n_threads=args.threads)
    sample.logger.info('Generating viewer transcript file.')
    _ = tx_converter(
        sample,
        out_path=out_dir / 'g4x_viewer' / f'{sample.sample_id}.tar',
        n_threads=args.threads,
        logger=sample.logger,
    )

    sample.logger.info('Completed redemux.')
