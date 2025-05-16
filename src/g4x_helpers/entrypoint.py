import argparse
import numpy as np
import pandas as pd
import geopandas
from pathlib import Path
from g4x_helpers.models import G4Xoutput

SUPPORTED_MASK_FILETYPES = ['.npz', '.npy', '.csv']

def launch_resegmentation():

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('--run_base', help='Path to G4X sample output folder', action='store', type=str, required=True)
    parser.add_argument('--segmentation_mask', help='Path to new segmentation mask', action='store', type=str, required=True)
    parser.add_argument('--sample_id', help='sample_id (Optional)', action='store', type=str, required=False, default=None)
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
    if segmentation_mask.suffix == '.npz' or segmentation_mask.suffix == '.npy':
        labels = np.load(segmentation_mask)
        if args.segmentation_mask_key:
            assert args.segmentation_mask_key in labels, print(f"{args.segmentation_mask_key} does not appear to be in provided mask file: {list(labels.keys())}.")
            labels = labels[args.segmentation_mask_key]
        assert labels.shape == sample.shape, print(f"provided mask shape ({labels.shape}) does not match with G4X sample shape ({sample.shape})")
    else:
        labels = geopandas.read_file(segmentation_mask)

    ## run intersection with new segmentation
    _= sample.intersect_segmentation(
        labels= labels,
        n_threads= args.threads
    )