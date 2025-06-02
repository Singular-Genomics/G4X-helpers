import os
import argparse
import numpy as np
import pandas as pd
import geopandas
from pathlib import Path
from g4x_helpers.utils import verbose_to_log_level
from g4x_helpers.models import G4Xoutput
from g4x_helpers.g4x_viewer.bin_generator import seg_converter, seg_updater

SUPPORTED_MASK_FILETYPES = ['.npz', '.npy', '.geojson']


def try_load_segmentation(
    segmentation_mask: Path, 
    expected_shape: tuple[int], 
    segmentation_mask_key: str | None = None
) -> np.ndarray | geopandas.GeoDataFrame:
    
    ## load new segmentation
    if segmentation_mask.suffix == '.npz' or segmentation_mask.suffix == '.npy':
        with np.load(segmentation_mask) as labels:
            if segmentation_mask_key:
                if segmentation_mask.suffix != '.npz':
                    raise ValueError("segmentation_mask_key was provided, but .npz file was not provided.")
                if segmentation_mask_key not in labels:
                    raise KeyError(f"{segmentation_mask_key} does not appear to be in provided mask file, available keys = {list(labels.keys())}.")
                seg = labels[segmentation_mask_key]
            else:
                seg = labels
            assert seg.shape == expected_shape, f"provided mask shape {seg.shape} does not match with G4X sample shape {expected_shape}"
    else:
        seg = geopandas.read_file(segmentation_mask)

    return seg


def launch_resegment():

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument('--segmentation_mask', help='Path to new segmentation mask. Supported files types are: .npy, .npz, .geojson.', action='store', type=str, required=True)
    parser.add_argument('--sample_id', help='sample_id (Optional)', action='store', type=str, required=False, default=None)
    parser.add_argument('--out_dir', help='Output directory where new files will be saved. Will be created if it does not exist. If not provided, the files in run_base will be updated in-place.', action='store', type=str, required=False, default=None)
    parser.add_argument('--segmentation_mask_key', help='Key in npz where mask should be taken from (Optional)', action='store', type=str, required=False, default=None)
    parser.add_argument('--threads', help='Number of threads to use for processing. [4]', action='store', type=int, required=False, default= 4)
    parser.add_argument('--verbose', help= 'Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]', action= 'store', type= int, default= 1)

    args = parser.parse_args()

    ## preflight checks
    run_base = Path(args.run_base)
    assert run_base.exists(), f"{run_base} does not appear to exist."
    segmentation_mask = Path(args.segmentation_mask)
    assert segmentation_mask.suffix in SUPPORTED_MASK_FILETYPES, f"{segmentation_mask.suffix} not a supported file type."

    ## initialize G4X sample
    sample = G4Xoutput(run_base=run_base, sample_id=args.sample_id, out_dir=args.out_dir, log_level=verbose_to_log_level(args.verbose))
    print(sample)

    ## load new segmentation
    labels = try_load_segmentation(segmentation_mask, sample.shape, args.segmentation_mask_key)

    ## run intersection with new segmentation
    _= sample.intersect_segmentation(
        labels= labels,
        out_dir= args.out_dir,
        n_threads= args.threads
    )


def launch_update_bin():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--bin_file', help='Path to G4X-Viewer segmentation bin file.', action='store', type=str, required=True)
    parser.add_argument('--out_path', help='Output file path', action='store', type=str, required=True)
    parser.add_argument('--metadata', help='Path to metadata table where clustering and/or embedding information will be extracted. Table must contain a header.', action='store', type=str, required=True)
    parser.add_argument('--cellid_key', help='Column name in metadata containing cell IDs that match with bin_file. If not provided, assumes that first column in metadata contains the cell IDs.', action='store', type=str, required=False, default=None)
    parser.add_argument('--cluster_key', help='Column name in metadata containing cluster IDs. If not provided, skips updating cluster IDs.', action='store', type=str, required=False, default=None)
    parser.add_argument('--emb_key', help='Column name in metadata containing embedding. Parser will look for {emb_key}_0 and {emb_key}_1. If not provided, skips updating embedding.', action='store', type=str, required=False, default=None)
    parser.add_argument('--verbose', help= 'Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]', action= 'store', type= int, default= 1)
    
    args = parser.parse_args()

    ## preflight checks
    bin_file = Path(args.bin_file)
    if not bin_file.exists():
        raise FileNotFoundError(f"{bin_file} does not appear to exist.")
    metadata = Path(args.metadata)
    if not metadata.exists():
        raise FileNotFoundError(f"{metadata} does not appear to exist.")
    out_path = Path(args.out_path)
    out_dir = out_path.parent
    os.makedirs(out_dir, exist_ok=True)

    ## run converter
    _= seg_updater(
        bin_file=bin_file,
        metadata=metadata,
        out_path=out_path,
        cluster_key=args.cluster_key,
        emb_key=args.emb_key,
        log_level=verbose_to_log_level(args.verbose)
    )