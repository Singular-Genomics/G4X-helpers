CLI_HELP = 'Helper models and post-processing tools for G4X-data\n\ndocs.singulargenomics.com'

RESEG_HELP = (
    'Reprocess G4X-data with a new segmentation\n\n'
    'Takes new cell-labels from a custom segmentation output and re-assigns transcripts and protein signals to those cells. '
    'The operation recreates most single-cell outputs and initializes a new G4X-viewer "segmentation.bin" file. '
    'Does not regenerate umaps/clustering or metrics.'
)

REDMX_HELP = (
    'Reprocess G4X-data with a new transcript manifest\n\n'
    'Generates a new "transcript_table.csv" by demultiplexing the raw feature data against a provided list of probe sequences '
    'and mapping each feature to its corresponding gene/target name. It then proceeds to regenerate single-cell outputs and initializes '
    'a new G4X-viewer "segmentation.bin" and "transcripts.tar" file.'
)


UDBIN_HELP = (
    'Update existing G4X-viewer segmentation file with new metadata\n\n'
    'Populates cluster or embedding information in an existing G4X-viewer segmentation.bin file. '
    'Only updates fields that are provided leaving other fields unchanged. Cell-IDs in the metadata table must match those in the bin file. '
    'Any cell-IDs in the metadata that are not found in the bin file will be assigned a default value.'
)

NWBIN_HELP = (
    'Initialize new segmentation file for G4X-viewer\n\n '
    'Converts cell-labels from "segmentation_mask.npz" into a "segmentation.bin" file compatible with the latest version of G4X-viewer. '
    'Cluster labels and umap embeddings are not populated in this operation. Please use "g4x-helpers update_bin" to add desired metadata.\n\n'
)

TARVW_HELP = (
    'Package G4X-viewer folder for distribution\n\n'
    'Creates a .tar archive of the "g4x_viewer" directory for easy upload and sharing.\n\n'
    'Archive file name: {SAMPLE_ID}_viewer.tar\n\n'
)

MIGRT_HELP = (
    'Migrate legacy data to comply with latest schemas for G4X-viewer & helpers\n\n'
    'Updates file names, headers and locations and creates new G4X-viewer files if necessary.\n\n'
    'A backup of the originals is created before migrating and operation will not proceed if an existing backup is present.\n\n'
)

VLDTE_HELP = 'Validate G4X-data to ensure correct file and folder structure\n\n'
