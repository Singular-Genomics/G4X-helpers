import logging
from importlib.metadata import version as ver
from pathlib import Path

import rich_click as click

from .models import G4Xoutput
from .utils import setup_logger, verbose_to_log_level

# rich_click.rich_click.THEME = 'modern'
click.rich_click.MAX_WIDTH = 100
click.rich_click.COMMANDS_PANEL_TITLE = 'commands'
click.rich_click.OPTIONS_PANEL_TITLE = 'options'
click.rich_click.STYLE_OPTION = 'bold blue'
click.rich_click.STYLE_ARGUMENT = 'bold blue'
click.rich_click.STYLE_COMMAND = 'bold blue'
click.rich_click.STYLE_SWITCH = 'bold red'
click.rich_click.STYLE_METAVAR = 'bold red'
click.rich_click.STYLE_METAVAR_SEPARATOR = 'dim'
click.rich_click.STYLE_USAGE = 'bold yellow'
click.rich_click.STYLE_USAGE_COMMAND = 'bold'
click.rich_click.STYLE_HELPTEXT_FIRST_LINE = ''
click.rich_click.STYLE_HELPTEXT = 'dim'
click.rich_click.STYLE_OPTION_DEFAULT = 'dim'
click.rich_click.STYLE_REQUIRED_SHORT = 'bold yellow'
click.rich_click.STYLE_REQUIRED_LONG = 'bold yellow'
click.rich_click.STYLE_OPTIONS_PANEL_BORDER = 'dim'
click.rich_click.STYLE_COMMANDS_PANEL_BORDER = 'dim'
# click.rich_click.COMMANDS_BEFORE_OPTIONS = True

click.rich_click.ARGUMENTS_PANEL_TITLE = 'input/output'

click.rich_click.COMMAND_GROUPS = {
    'g4x-helpers': [
        {'name': 'commands', 'commands': ['redemux', 'resegment', 'update-bin', 'new-bin', 'tar-viewer']},
        # {"name": "utilities", "commands": ["log"]},
    ],
}

CLI_HELP = 'Helper models and post-processing tools for G4X-data\n\ndocs.singulargenomics.com'
RESEG_HELP = 'Reprocess G4X-output with a new segmentation'
REDMX_HELP = 'Reprocess G4X-output with a new transcript manifest'
UDBIN_HELP = 'Update G4X-viewer .bin file with new metadata'
NWBIN_HELP = 'Generate G4X-viewer .bin files from sample output'
TARVW_HELP = 'Package G4X-viewer folder for distribution'


# region cli
@click.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    invoke_without_command=True,
    add_help_option=True,
    help=CLI_HELP,
)
@click.argument(
    'sample-dir',
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help='Path to G4X sample output folder passed to commands.',
    # panel='input/output',
)
@click.argument(
    'out-dir',
    required=False,
    type=click.Path(exists=False, file_okay=False, dir_okay=True, writable=True, path_type=Path),
    default='./g4x_helpers',
    help='Output directory used by downstream commands.',
    # panel='input/output',
)
@click.option(
    '--sample-id',  #
    required=False,
    type=str,
    help='Sample ID (optional).',
)
@click.option(
    '--threads',
    required=False,
    type=int,
    default=4,
    show_default=True,
    help='Number of threads to use for processing.',
)
@click.option(
    '--verbose',
    required=False,
    type=int,
    default=1,
    show_default=True,
    help='Set logging level WARNING (0), INFO (1), or DEBUG (2).',
)
@click.option(
    '--version',
    is_flag=True,
    default=False,
    help='Display g4x-helpers version.',
)
@click.pass_context
def cli(ctx, sample_dir, out_dir, sample_id, threads, verbose, version):
    if version:
        v = ver('g4x_helpers')
        click.echo(f'g4x-helpers version: {v}')

    else:
        ctx.ensure_object(dict)

        out_dir.mkdir(parents=True, exist_ok=True)

        ctx.obj['sample_dir'] = sample_dir
        ctx.obj['out_dir'] = out_dir
        ctx.obj['threads'] = threads
        ctx.obj['verbose'] = verbose

        try:
            sample = G4Xoutput(run_base=ctx.obj['sample_dir'], sample_id=sample_id)
        except Exception as e:
            click.echo('')
            click.secho('Failed to load G4Xoutput:', fg='red', err=True, bold=True)
            raise click.ClickException(f'{type(e).__name__}: {e}')

        click.echo(sample)

        # No subcommand given → show help
        if ctx.invoked_subcommand is None:
            click.echo(ctx.get_help())
            # print('')
            ctx.exit()


############################################################
# region resegment
@cli.command(name='resegment', help=RESEG_HELP)
@click.argument(
    'segmentation-mask',
    panel='commands',
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to new segmentation mask. Supported file types: .npy, .npz, .geojson.',
)
@click.option(
    '--sample-id',  #
    required=False,
    type=str,
    default=None,
    help='Sample ID (optional).',
)
@click.option(
    '--segmentation-mask-key',
    required=False,
    type=str,
    default=None,
    help='Key in npz/geojson where mask/labels should be taken from (optional).',
)
@click.pass_context
def resegment(ctx, segmentation_mask, sample_id, segmentation_mask_key):
    from . import segmentation as seg

    click.echo('')
    click.echo(f'Resegmentation (threads={ctx.obj["threads"]})')

    out_dir = ctx.obj['out_dir']

    logger = setup_logger(
        logger_name='g4x_helpers',
        stream_logger=True,
        stream_level=verbose_to_log_level(ctx.obj['verbose']),
        file_logger=True,
        file_level=logging.DEBUG,
        out_dir=out_dir,
        clear_handlers=True,
    )

    # initialize G4X sample
    logger.info(f'Initializing G4X sample with run_base: {ctx.obj["sample_dir"]}, sample_id: {sample_id}')
    sample = G4Xoutput(run_base=ctx.obj['sample_dir'], sample_id=sample_id)
    click.echo(sample)

    # load new segmentation
    labels = seg.try_load_segmentation(segmentation_mask, sample.shape, segmentation_mask_key)
    logger.info(f'Loaded new segmentation with shape: {labels.shape}, n_labels: {len(set(labels.flatten()))}')

    # run intersection with new segmentation
    seg.intersect_segmentation(
        g4x_out=sample, labels=labels, out_dir=out_dir, n_threads=ctx.obj['threads'], logger=logger
    )
    click.echo('Resegmentation complete.')


############################################################
# region redemux
@cli.command(name='redemux', help=REDMX_HELP)
@click.argument(
    'manifest',
    panel='commands',
    required=False,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to manifest for demuxing. The manifest must be a 3-column CSV with the following header: target,sequence,read.',
)
@click.option(
    '--batch-size',
    default=1_000_000,
    show_default=True,
    type=int,
    help='Number of transcripts to process per batch.',
)
@click.pass_context
def redemux(ctx, manifest, batch_size):
    from . import demux as dmx

    click.echo('')
    click.echo(f'Redemuxing (threads={ctx.obj["threads"]})')

    click.echo(f'g4x-helpers sample_dir: {ctx.obj["sample_dir"]}')
    click.echo(f'g4x-helpers out_dir: {ctx.obj["out_dir"]}')

    # out_dir = ctx.obj['out_dir']

    # logger = setup_logger(
    #     logger_name='redemux',
    #     stream_logger=True,
    #     stream_level=verbose_to_log_level(ctx.obj['verbose']),
    #     file_logger=True,
    #     file_level=logging.DEBUG,
    #     out_dir=out_dir,
    #     clear_handlers=True,
    # )

    # # initialize G4X sample
    # logger.info(f'Initializing G4X sample with run_base: {ctx.obj["sample_dir"]}')
    # sample = G4Xoutput(run_base=ctx.obj['sample_dir'])
    # click.echo(sample)

    # # load new segmentation
    # dmx.redemux(
    #     g4x_out=sample,
    #     manifest=manifest,
    #     batch_size=batch_size,
    #     out_dir=out_dir,
    #     threads=ctx.obj['threads'],
    #     logger=logger,
    # )

    click.echo('Redemuxing complete.')


############################################################
# region update_bin
@cli.command(name='update-bin', help=UDBIN_HELP)
@click.option(
    '--bin-file', required=True, type=click.Path(exists=True), help='Path to G4X-Viewer segmentation bin file.'
)
@click.option(
    '--metadata',
    required=True,
    type=click.Path(exists=True),
    help='Path to metadata table where clustering and/or embedding information will be extracted. Table must contain a header.',
)
@click.option(
    '--cellid-key',
    default=None,
    type=str,
    help='Column name in metadata containing cell IDs that match with bin_file. If not provided, assumes that first column in metadata contains the cell IDs.',
)
@click.option(
    '--cluster-key',
    default=None,
    type=str,
    help='Column name in metadata containing cluster IDs. Automatically assigns new colors if cluster_color_key is not provided. If not provided, skips updating cluster IDs.',
)
@click.option(
    '--cluster-color-key',
    default=None,
    type=str,
    help='Column name in metadata containing cluster colors. Colors must be provided as hex codes. If provided, cluster_key must also be provided. If not provided, skips updating cluster colors.',
)
@click.option(
    '--emb-key',
    default=None,
    type=str,
    help='Column name in metadata containing embedding. Parser will look for {emb_key}_1 and {emb_key}_2. If not provided, skips updating embedding.',
)
@click.pass_context
def update_bin(ctx, bin_file, metadata, cellid_key, cluster_key, cluster_color_key, emb_key):
    from .g4x_viewer.bin_generator import seg_updater

    out_dir = ctx.obj['out_dir']
    out_dir = out_dir.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    _ = seg_updater(
        bin_file=bin_file,
        metadata_file=metadata,
        out_path=out_dir,
        cellid_key=cellid_key,
        cluster_key=cluster_key,
        cluster_color_key=cluster_color_key,
        emb_key=emb_key,
        log_level=verbose_to_log_level(ctx.obj['verbose']),
    )


############################################################
# region new_bin
@cli.command(name='new-bin', help=NWBIN_HELP)
@click.pass_context
def new_bin(ctx):
    from g4x_helpers.g4x_viewer.bin_generator import seg_converter

    run_base = ctx.obj['sample_dir']
    out_dir = ctx.obj['out_dir']

    ## initialize G4X sample
    sample = G4Xoutput(run_base=run_base)
    print(sample)

    ## set up the data
    try:
        adata = sample.load_adata()
        emb_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
        cluster_key = sorted([x for x in adata.obs.columns if 'leiden' in x])[0]
        sample.logger.info('Successfully loaded adata with clustering information.')
    except Exception:
        adata = sample.load_adata(load_clustering=False)
        emb_key = None
        cluster_key = None
        sample.logger.info('Clustering information was not found, cell coloring will be random.')

    mask = sample.load_segmentation()

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        outfile = Path(out_dir) / f'{sample.sample_id}.bin'
    else:
        outfile = run_base / 'g4x_viewer' / f'{sample.sample_id}.bin'

    sample.logger.info('Making G4X-Viewer bin file.')
    _ = seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=outfile,
        cluster_key=cluster_key,
        emb_key=emb_key,
        protein_list=[f'{x}_intensity_mean' for x in sample.proteins],
        n_threads=ctx.obj['n_threads'],
        logger=sample.logger,
    )
    sample.logger.debug(f'G4X-Viewer bin --> {outfile}')


############################################################
# region tar_viewer
@cli.command(name='tar-viewer', help=TARVW_HELP)
@click.option(
    '--viewer-dir',
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help='Path to G4X-viewer folder to tar.',
)
@click.pass_context
def tar_viewer(ctx, viewer_dir):
    from . import tar_viewer as tv

    out_dir = ctx.obj['out_dir']

    tv.tar_viewer(out_path=out_dir, viewer_dir=viewer_dir)


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
