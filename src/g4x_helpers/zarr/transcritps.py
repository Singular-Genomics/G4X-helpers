import numpy as np
import polars as pl
from numcodecs import Blosc


def write_transcripts(smp, root_group):
    tx_group = setup_tx_group(root_group)

    aggregation_level = 'gene_name'
    keep_cols = ['x_pixel_coordinate', 'y_pixel_coordinate', 'cell_id', aggregation_level]
    lf = smp.load_transcript_table(lazy=True, columns=keep_cols)
    df = lf.collect()

    img_resolution = smp.shape

    tile_specs = choose_square_tiling(
        image_resolution_hw=img_resolution,
        total_points=df.height,
        target_points_per_tile=5000,
        min_tile_size=256,
    )
    print(tile_specs)

    pyramid = build_tx_pyramid(tile_specs, image_resolution=img_resolution)

    for level, specs in pyramid.items():
        print(f'Level {level}: tile_size: {specs["tile_size"]} - scale: {specs["scale_fct"]}')

    tx_group.attrs['layer_config'] = {
        'layers': len(pyramid) - 1,
        'tile_size': pyramid[len(pyramid) - 1]['tile_size'],
        'layer_height': smp.shape[0],
        'layer_width': smp.shape[1],
    }

    gene_list = lf.unique('gene_name').sort('gene_name').collect()['gene_name'].to_list()

    N = len(gene_list)  # number of colors
    colors = np.random.randint(0, 256, size=(N, 3), dtype=np.uint8)
    gene_colors = {g: c.tolist() for g, c in zip(gene_list, colors)}
    tx_group.attrs['gene_colors'] = gene_colors

    pyramid = construct_tile_dfs(df, pyramid)
    write_tx_zarr(tx_group, pyramid)


def setup_tx_group(root_group):
    tx_group = root_group.create_group('transcripts', overwrite=True)

    return tx_group


def choose_square_tiling(image_resolution_hw, total_points, target_points_per_tile, min_tile_size=64):
    tiles_needed = np.ceil(total_points / target_points_per_tile).astype(int)

    H, W = image_resolution_hw  # (H, W)
    s = max(1, int(np.sqrt((W * H) / tiles_needed)))  # start guess

    def _num_tiles(s):
        nx = np.ceil(W / s)
        ny = np.ceil(H / s)
        return nx, ny

    while True:
        nx, ny = _num_tiles(s)
        if nx * ny >= tiles_needed:
            break
        s -= 1

    if s < min_tile_size:
        s = min_tile_size
        nx, ny = _num_tiles(s)

    res = {'tile_size': int(s), 'nx': int(nx), 'ny': int(ny), 'grid_total': int(nx * ny)}
    return res


def build_tx_pyramid(tile_specs, image_resolution):
    pyramid = {}
    level = 0
    tile_size = tile_specs['tile_size']

    # build pyramid until last tile size exceeds image resolution
    while tile_size <= max(image_resolution):
        scale = 2**level
        tile_size = tile_specs['tile_size'] * scale
        sampling_fct = 1 / (scale**2)
        pyramid[level] = {'tile_size': tile_size, 'scale_fct': scale, 'sampling_fct': sampling_fct}
        level += 1

    return pyramid


def construct_tile_dfs(df, pyramid):
    for level, specs in pyramid.items():
        tile_size, _, sampling_fct = specs.values()
        df_smp = (
            df.sample(fraction=sampling_fct, with_replacement=False)
            .with_columns(
                (pl.col('x_pixel_coordinate') / tile_size).cast(pl.Int32).alias('tile_y'),
                (pl.col('y_pixel_coordinate') / tile_size).cast(pl.Int32).alias('tile_x'),
            )
            .sort('tile_y', 'tile_x')
        )
        pyramid[level].update({'tile_df': df_smp})
    return pyramid


def write_tx_zarr(tx_group, pyramid):
    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    for level in pyramid:
        tile_df = pyramid[level]['tile_df']
        all_coords = tile_df.select(['x_pixel_coordinate', 'y_pixel_coordinate']).to_numpy()
        all_cell_ids = tile_df.select('cell_id').to_numpy()
        all_gene_names = tile_df.select('gene_name').to_numpy()
        all_tile_ids = tile_df.select(['tile_y', 'tile_x']).to_numpy()

        tile_ids = tile_df.unique(['tile_y', 'tile_x'], maintain_order=True).select(['tile_y', 'tile_x']).to_numpy()

        for tile in tile_ids:
            y = str(tile[0]).zfill(2)
            x = str(tile[1]).zfill(2)
            tile_path = f'p{level}/y{y}/x{x}'
            tile_group = tx_group.require_group(tile_path)

            idx = np.where((all_tile_ids == tile).all(axis=1))[0]

            coords = all_coords[idx].astype(np.int32)
            gene_names = all_gene_names[idx].astype('U10')
            cell_ids = all_cell_ids[idx].astype(np.int32)

            tile_group.create_array('position', data=coords, compressor=compressor)
            tile_group.create_array('gene_name', data=gene_names, compressor=compressor)
            tile_group.create_array('cell_id', data=cell_ids, compressor=compressor)
