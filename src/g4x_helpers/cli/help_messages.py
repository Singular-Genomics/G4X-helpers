from .. import constants as c

CLI_HELP = 'Post-processing tools and utilities for G4X-data\n\ndocs.singulargenomics.com'

RESEGMENT_HELP = (
    'Reprocess G4X-data with a new segmentation\n\n'
    'Takes new cell-labels from a custom segmentation output and re-assigns transcripts and protein signals to those cells. '
    'The operation recreates single-cell outputs and initializes a new G4X-viewer zarr store. '
    'Does not regenerate metrics.'
)

REDEMUX_HELP = (
    'Reprocess G4X-data with a new transcript manifest\n\n'
    'Generates a new "transcript_table.csv.gz" by demultiplexing the raw feature data against a provided list of probe sequences '
    'and mapping each feature to its corresponding target gene name. It then proceeds to regenerate single-cell outputs and initializes '
    'a new G4X-viewer zarr store. '
    'Does not regenerate metrics.'
)


MIGRATE_HELP = (
    'Migrate legacy data to the latest schemas for G4X-viewer & helpers\n\n'
    'Moves raw data files to a new location and applies standard post-processing steps to generate outputs with an updated format and a g4x-viewer.zarr compatible with G4X-viewer v4.'
)

VALIDATE_HELP = 'Validate G4X-data to ensure correct file and folder structure\n\n'

VIEWER_HELP = (
    'Modify metadata in a G4X-viewer zarr store\n\n'
    'Each data layer in the g4x-viewer.zarr has its own subcommand with options to import or export metadata'
)

shared_metadata_help = (
    '\n\nIt is recommended to inspect and modify an exported file.'
    '\n\n--import and --export options cannot be used together'
)

CELLS_META_HELP = (
    'Modify cell metadata in a G4X-viewer zarr store\n\n'
    f'A valid metadata file must include the following columns:\n\n'
    f'[{c.CELL_ID_NAME}, UMAP1, UMAP2] and at least one pair of [label, label_color]\n\n'
) + shared_metadata_help

IMAGES_META_HELP = (
    'Modify image metadata in a G4X-viewer zarr store\n\n'
    'A valid metadata file must include the following columns:\n\n'
    '[label, active, color, min, max, start, end] where min/max define the intensity range and start/end the pre-selected channel limits\n\n'
) + shared_metadata_help

TRANSCRIPTS_META_HELP = (
    'Modify transcript metadata in a G4X-viewer zarr store\n\n'
    'A valid metadata file must include the following columns: [gene_id, color]\n\n'
) + shared_metadata_help
