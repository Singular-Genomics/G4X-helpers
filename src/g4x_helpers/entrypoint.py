import argparse
import numpy as np
import pandas as pd
import geopandas
from pathlib import Path
from g4x_helpers.models import G4Xoutput
from g4x_viewer.bin_generator import seg_converter

SUPPORTED_MASK_FILETYPES = ['.npz', '.npy', '.geojson']


def try_load_segmentation(
    segmentation_mask: Path, 
    expected_shape: tuple[int], 
    segmentation_mask_key: str | None = None
) -> np.ndarray | geopandas.GeoDataFrame:
    
    ## load new segmentation
    if segmentation_mask.suffix == '.npz' or segmentation_mask.suffix == '.npy':
        labels = np.load(segmentation_mask)
        if segmentation_mask_key:
            assert segmentation_mask_key in labels, print(f"{segmentation_mask_key} does not appear to be in provided mask file: {list(labels.keys())}.")
            labels = labels[segmentation_mask_key]
        assert labels.shape == expected_shape, print(f"provided mask shape ({labels.shape}) does not match with G4X sample shape ({expected_shape})")
    else:
        labels = geopandas.read_file(segmentation_mask)

    return labels


def launch_resegment():

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument('--segmentation_mask', help=f'Path to new segmentation mask. Supported files types are {SUPPORTED_MASK_FILETYPES}.', action='store', type=str, required=True)
    parser.add_argument('--sample_id', help='sample_id (Optional)', action='store', type=str, required=False, default=None)
    parser.add_argument('--out_dir', help='Output directory where new files will be saved. If not provided, the files in run_base will be updated in-place.', action='store', type=str, required=False, default=None)
    parser.add_argument('--segmentation_mask_key', help='Key in npz where mask should be taken from (Optional)', action='store', type=str, required=False, default=None)
    parser.add_argument('--threads', help='Number of threads to use for processing. [4]', action='store', type=int, required=False, default= 4)
    parser.add_argument('--verbose', help= 'Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]', action= 'store', type= int, default= 1)

    args = parser.parse_args()

    ## preflight checks
    run_base = Path(args.run_base)
    assert run_base.exists(), print(f"{run_base} does not appear to exist.")
    segmentation_mask = Path(args.segmentation_mask)
    assert segmentation_mask.suffix in SUPPORTED_MASK_FILETYPES, print(f"{segmentation_mask.suffix} not a supported file type.")

    ## initialize G4X sample
    sample = G4Xoutput(run_base=run_base, sample_id=args.sample_id, log_level=args.verbose)
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

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument('--out_path', help='Output file path', action='store', type=str, required=True)
    parser.add_argument('--segmentation_mask', help='Path to a custom segmentation mask. Supported files types are {SUPPORTED_MASK_FILETYPES}. If not provided, will use the default one in the G4X output folder.', action='store', type=str, required=True)
    parser.add_argument('--metadata', help='Path to metadata table where clustering and/or embedding information will extracted.', action='store', type=str, required=True)
    parser.add_argument('--cluster_key', help='Column name in metadata containing cluster IDs.', action='store', type=str, required=False, default=None)
    parser.add_argument('--emb_key', help='Column name in metadata containing embedding. Parser will look for emb_key_0 and emb_key_1', action='store', type=str, required=False, default=None)
    parser.add_argument('--threads', help='Number of threads to use for processing. [4]', action='store', type=int, required=False, default= 4)
    parser.add_argument('--verbose', help= 'Set logging level WARNING (0), INFO (1), or DEBUG (2). [1]', action= 'store', type= int, default= 1)
    
    args = parser.parse_args()

    ## preflight checks
    run_base = Path(args.run_base)
    assert run_base.exists(), print(f"{run_base} does not appear to exist.")
    if args.segmentation_mask:
        segmentation_mask = Path(args.segmentation_mask)
        assert segmentation_mask.suffix in SUPPORTED_MASK_FILETYPES, print(f"{segmentation_mask.suffix} not a supported file type.")

    ## initialize G4X sample
    sample = G4Xoutput(run_base=run_base, log_level=args.verbose)
    print(sample)

    ## load segmentation mask
    if args.segmentation_mask:
        labels = try_load_segmentation(segmentation_mask, sample.shape, args.segmentation_mask_key)
    else:
        labels = sample.load_segmentation(expanded= True)
    
    ## run converter
    _= seg_converter(
        adata= sample.load_adata(load_clustering= False),
        seg_mask= labels,
        outpath= args.out_path,
        metadata= args.metadata,
        cluster_key= args.cluster_key,
        emb_key= args.emb_key,
        protein_list= [f"{x}_intensity_mean" for x in sample.proteins],
        n_threads= args.threads
    )