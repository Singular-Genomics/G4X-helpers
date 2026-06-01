from ...modules.viewer import cells


def cell_metadata(viewer_zarr, import_metadata, export_metadata, segmentation):
    if import_metadata is not None and export_metadata is not None:
        raise ValueError('--import-metadata and --export-metadata cannot be used together.')

    seg_group = cells.get_seg_group(viewer_zarr, segmentation)

    if export_metadata is not None:
        meta = cells.get_cell_metadata(seg_group)
        meta.write_csv(export_metadata)
    if import_metadata is not None:
        cells.apply_viewer_metadata(seg_group, import_metadata)
