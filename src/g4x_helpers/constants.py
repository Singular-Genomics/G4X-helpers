import os

# System
DEFAULT_THREADS = max(1, (os.cpu_count() // 2 or 4))

# Physical
PIXEL_SIZE_MICRONS = 0.3125
PIXEL_PER_MICRON = 3.2

# Demux
PROBE_PATTERN = r'^(.*?)-([ACGT]{2,30})-([^-]+)$'

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
