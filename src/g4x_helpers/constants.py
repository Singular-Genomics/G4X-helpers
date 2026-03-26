from .utils import default_workers

# System
# TODO this isn't actually a constant, so it should go somewhere else
DEFAULT_THREADS = default_workers(max_workers=16, reserve=1)

# Files
SMP_META = 'sample.g4x'
TX_PANEL = 'transcript_panel.csv'
PR_PANEL = 'protein_panel.csv'
SUMMARY = 'summary_*.html'
SSHEET = 'samplesheet.csv'

NUCLEAR_STAIN = 'nuclear'
CYTOPLASMIC_STAIN = 'cytoplasmic'

HE_DIR = 'h_and_e'
NUC_IMG = f'{HE_DIR}/{NUCLEAR_STAIN}'
CYT_IMG = f'{HE_DIR}/{CYTOPLASMIC_STAIN}'
HNE_IMG = f'{HE_DIR}/h_and_e'

PR_DIR = 'protein'

RNA_DIR = 'rna'
RAW_FEATURES = f'{RNA_DIR}/raw_features.parquet'
FILE_TX_TABLE = f'{RNA_DIR}/transcript_table.csv.gz'

MASK_DIR = 'masks'
SEG_MASK = f'{MASK_DIR}/segmentation_mask.npz'
BEAD_MASK = f'{MASK_DIR}/bead_mask.npz'

FILE_VIEWER_ZARR = 'g4x-viewer.zarr'

SINGLE_CELL_DIR = 'single_cell_data'
FILE_CELL_METADATA = f'{SINGLE_CELL_DIR}/cell_metadata.csv.gz'
FILE_CELL_X_GENE = f'{SINGLE_CELL_DIR}/cell_by_gene.csv.gz'
FILE_CELL_X_PROTEIN = f'{SINGLE_CELL_DIR}/cell_by_protein.csv.gz'
FILE_CLUSTERING_UMAP = f'{SINGLE_CELL_DIR}/clustering_umap.csv.gz'
FILE_DGEX = f'{SINGLE_CELL_DIR}/dgex.csv.gz'
FILE_FEAT_MTX = f'{SINGLE_CELL_DIR}/sc_processed.h5ad'

PREFERRED_IMG_SUFFIX = '.ome.tiff'
ALT_IMG_SUFFIX = '.jp2'

# Physical
PIXEL_SIZE_MICRONS = 0.3125
PIXEL_PER_MICRON = 3.2

# Demux
PROBE_PATTERN = r'^(.*?)-([ACGT]{2,30})-([^-]+)$'
BASE_ORDER = 'CTGA'
DEFAULT_BATCH_SIZE = 1_000_000

primer_read_map = {
    'SP1': 1,
    'm7a': 2,
    'm9a': 3,
    'm6a': 4,
    'm8a': 5,
    'm3a': 6,
}

# Metadata column names
CELL_ID_NAME = 'seg_cell_id'
GENE_ID_NAME = 'gene_id'
CELL_COORD_X = 'cell_x'
CELL_COORD_Y = 'cell_y'

# Cell colors in G4X-viewer
UNASSIGNED_COLOR = '#BFBFBF'

SG_PALETTE = [
    '#72FFAB',
    '#A16CFD',
    '#FF7043',
    '#008FFF',
    '#D32F2F',
    '#7CB342',
    '#7F34BE',
    '#FFCA28',
    '#0C8668',
    '#FB4695',
    '#005EE1',
    '#28EDED',
    '#A17B64',
    '#FFFF58',
    '#BC29AE',
    '#006D8F',
    '#FFBAFF',
    '#FFD091',
    '#5C6BC0',
    '#F490B2',
]
