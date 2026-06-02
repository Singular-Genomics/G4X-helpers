CLI_HELP = 'Post-processing tools and utilities for G4X-data\n\ndocs.singulargenomics.com'

RESEGMENT_HELP = (
    'Reprocess G4X-data with a new segmentation\n\n'
    'Takes new cell-labels from a custom segmentation output and re-assigns transcripts and protein signals to those cells. '
    'The operation recreates most single-cell outputs and initializes a new G4X-viewer "segmentation.bin" file. '
    'Does not regenerate umaps/clustering or metrics.'
)

REDEMUX_HELP = (
    'Reprocess G4X-data with a new transcript manifest\n\n'
    'Generates a new "transcript_table.csv" by demultiplexing the raw feature data against a provided list of probe sequences '
    'and mapping each feature to its corresponding gene/target name. It then proceeds to regenerate single-cell outputs and initializes '
    'a new G4X-viewer "segmentation.bin" and "transcripts.tar" file.'
)


MIGRATE_HELP = (
    'Migrate legacy data to the latest schemas for G4X-viewer & helpers\n\n'
    'Moves raw data files to a new location and applies standard post-processing steps to generate single-cell outputs and a g4x-viewer.zarr compatible with G4X-viewer v4.'
)

VALIDATE_HELP = 'Validate G4X-data to ensure correct file and folder structure\n\n'

VIEWER_HELP = (
    'Modify metadata in a G4X-viewer zarr store\n\n'
    'Each data layer in the viewer has its own subcommand with options to import or export metadata'
)
