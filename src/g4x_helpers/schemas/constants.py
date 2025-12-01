files_to_migrate = {
    'bin': ('g4x_viewer/{sample_id}.bin', 'g4x_viewer/{sample_id}_segmentation.bin', True),
    'tar': ('g4x_viewer/{sample_id}.tar', 'g4x_viewer/{sample_id}_transcripts.tar', True),
    'ome.tiff': ('g4x_viewer/{sample_id}.ome.tiff', 'g4x_viewer/{sample_id}_multiplex.ome.tiff', False),
    'raw_tx_v2': ('diagnostics/transcript_table.parquet', 'rna/raw_features.parquet', True),
    'raw_tx_v3': ('rna/transcript_table.parquet', 'rna/raw_features.parquet', True),
    'tx_table': ('rna/transcript_table.csv.gz', 'rna/transcript_table.csv.gz', True),
    'tx_panel': ('transcript_panel.csv', 'transcript_panel.csv', True),
    'feat_mtx': ('single_cell_data/feature_matrix.h5', 'single_cell_data/feature_matrix.h5', True),
    'cell_meta': ('single_cell_data/cell_metadata.csv.gz', 'single_cell_data/cell_metadata.csv.gz', True),
}
