import json
from pathlib import Path

import polars as pl

# required_files
SAMPLE_META = Path('sample_meta.json')
TX_PANEL = Path('transcript_panel.csv')
PR_PANEL = Path('protein_panel.csv')
RAW_FEATURES = Path('rna/raw_features.parquet')
SEG_MASK = Path('segmentation/segmentation_mask.npz')
NUC_IMAGE = Path('h_and_e/nuclear.jp2')
H_AND_E_IMAGE = Path('h_and_e/h_and_e.jp2')
CYT_IMAGE = Path('h_and_e/eosin.jp2')


class ValidationError(Exception):
    pass


def validate_raw_data(data_dir):
    data_dir = Path(data_dir)
    meta_path = data_dir / SAMPLE_META
    if not meta_path.exists():
        raise ValidationError(f'G4X-helpers requires a sample_meta.json file: {meta_path}')

    with open(meta_path, 'r') as f:
        run_meta = json.load(f)

    # validate H&E images
    he_img_path = data_dir / 'h_and_e'
    if not he_img_path.exists():
        raise ValidationError(f'H&E image directory not found: {he_img_path}')

    available_images = [f.stem for f in (he_img_path).glob('*.jp2')]
    listed_proteins = ['eosin', 'h_and_e', 'nuclear']
    missing = [p for p in listed_proteins if p not in available_images]

    if len(missing) > 0:
        raise ValidationError(f'Missing h_and_e images for: {missing}')

    # TODO ensure this is in the meta
    # validate transcript panel
    if 'transcript_panel' in run_meta:
        tx_panel_path = data_dir / TX_PANEL
        raw_feats = data_dir / RAW_FEATURES

        if not tx_panel_path.exists():
            raise ValidationError(f'Transcript panel file not found: {tx_panel_path}')

        if not raw_feats.exists():
            raise ValidationError(f'Raw features file not found: {raw_feats}')

    # TODO ensure this is in the meta
    # validate protein panel
    if 'protein_panel' in run_meta:
        pr_panel_path = data_dir / PR_PANEL
        pr_img_path = data_dir / 'protein'

        if not pr_panel_path.exists():
            raise ValidationError(f'Protein panel file not found: {pr_panel_path}')

        if not pr_img_path.exists():
            raise ValidationError(f'Protein image directory not found: {pr_img_path}')

        protein_panel = pl.read_csv(pr_panel_path).sort(
            by=['panel_type', pl.col('target').str.to_lowercase()], descending=[True, False]
        )

        available_images = [f.stem for f in pr_img_path.glob('*.jp2')]
        listed_proteins = protein_panel['target'].to_list()
        missing = [p for p in listed_proteins if p not in available_images]

        if len(missing) > 0:
            raise ValidationError(f'Missing protein images for: {missing}')

    # validate segmentation
    seg_file = data_dir / SEG_MASK
    if not seg_file.parent.exists():
        raise ValidationError(f'Segmentation directory not found: {seg_file.parent}')

    if not seg_file.exists():
        raise ValidationError(f'Segmentation file not found: {seg_file}')

    return run_meta
