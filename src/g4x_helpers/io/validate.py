import json
from pathlib import Path

import polars as pl


def validate_raw_data(data_dir):
    data_dir = Path(data_dir)
    meta_path = data_dir / 'sample_meta.json'
    if not meta_path.exists():
        raise FileNotFoundError(f'G4X-helpers requires a sample_meta.json file: {meta_path}')

    with open(meta_path, 'r') as f:
        run_meta = json.load(f)

    # validate H&E images
    he_img_path = data_dir / 'h_and_e'
    if not he_img_path.exists():
        raise FileNotFoundError(f'H&E image directory not found: {he_img_path}')

    available_images = [f.stem for f in (he_img_path).glob('*.jp2')]
    listed_proteins = ['eosin', 'h_and_e', 'nuclear']
    missing = [p for p in listed_proteins if p not in available_images]

    if len(missing) > 0:
        raise FileNotFoundError(f'Missing h_and_e images for: {missing}')

    # TODO ensure this is in the meta
    # validate transcript panel
    if 'transcript_panel' in run_meta:
        tx_panel_path = data_dir / 'transcript_panel.csv'
        raw_feats = data_dir / 'rna' / 'raw_features.parquet'

        if not tx_panel_path.exists():
            raise FileNotFoundError(f'Transcript panel file not found: {tx_panel_path}')

        if not raw_feats.exists():
            raise FileNotFoundError(f'Raw features file not found: {raw_feats}')

    # TODO ensure this is in the meta
    # validate protein panel
    if 'protein_panel' in run_meta:
        pr_panel_path = data_dir / 'protein_panel.csv'
        pr_img_path = data_dir / 'protein'

        if not pr_panel_path.exists():
            raise FileNotFoundError(f'Protein panel file not found: {pr_panel_path}')

        if not pr_img_path.exists():
            raise FileNotFoundError(f'Protein image directory not found: {pr_img_path}')

        protein_panel = pl.read_csv(pr_panel_path).sort(
            by=['panel_type', pl.col('target').str.to_lowercase()], descending=[True, False]
        )

        available_images = [f.stem for f in pr_img_path.glob('*.jp2')]
        listed_proteins = protein_panel['target'].to_list()
        missing = [p for p in listed_proteins if p not in available_images]

        if len(missing) > 0:
            raise FileNotFoundError(f'Missing protein images for: {missing}')

    # validate segmentation
    seg_path = data_dir / 'segmentation'
    if not seg_path.exists():
        raise FileNotFoundError(f'Segmentation directory not found: {seg_path}')

    seg_file = seg_path / 'segmentation_mask.npz'
    if not seg_file.exists():
        raise FileNotFoundError(f'Segmentation file not found: {seg_file}')

    return run_meta
