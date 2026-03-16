from .utils import default_workers

# System
# TODO this isn't actually a constant, so it should go somewhere else
DEFAULT_THREADS = default_workers(max_workers=16, reserve=1)

# Files
REQUIRED_SMP_META = 'sample.g4x'
REQUIRED_TX_PANEL = 'transcript_panel.csv'
REQUIRED_PR_PANEL = 'protein_panel.csv'
REQUIRED_RAW_FEATURES = 'rna/raw_features.parquet'
REQUIRED_SEG_MASK = 'masks/segmentation_mask.npz'
REQUIRED_BEAD_MASK = 'masks/bead_mask.npz'
REQUIRED_SUMMARY = 'summary_*.html'
REQUIRED_SSHEET = 'samplesheet.csv'

REQUIRED_NUC_IMG = 'h_and_e/nuclear.ome.tiff'
REQUIRED_CYT_IMG = 'h_and_e/cytoplasmic.ome.tiff'
REQUIRED_HnE_IMG = 'h_and_e/h_and_e.ome.tiff'

REQUIRED_HE_DIR = 'h_and_e'
REQUIRED_PR_DIR = 'protein'
REQUIRED_PR_SUFFIX = '.ome.tiff'
ALT_PR_SUFFIX = '.jp2'

DIRECTORY_SINGLE_CELL = 'single_cell_data'

FILE_TX_TABLE = 'rna/transcript_table.csv.gz'
FILE_VIEWER_ZARR = 'g4x-viewer.zarr'

FILE_CELL_METADATA = 'single_cell_data/cell_metadata.csv.gz'
FILE_CELL_X_GENE = 'single_cell_data/cell_by_gene.csv.gz'
FILE_CELL_X_PROTEIN = 'single_cell_data/cell_by_protein.csv.gz'
FILE_FEAT_MTX = 'single_cell_data/feature_matrix.h5'
FILE_CLUSTERING_UMAP = 'single_cell_data/clustering_umap.csv.gz'


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
