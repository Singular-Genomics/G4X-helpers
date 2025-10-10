import rich_click as click

from . import __version__, utils
from . import main_features as main

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
    required=False,
    type=click.Path(exists=True, file_okay=False),
    # default='.',
    help='Path to G4X sample output folder passed to commands.',
    # panel='input/output',
)
@click.argument(
    'out-dir',
    required=False,
    type=click.Path(exists=False, file_okay=False, dir_okay=True, writable=True),
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
    '-t',
    '--threads',
    required=False,
    type=int,
    default=utils.DEFAULT_THREADS,
    show_default=True,
    help='Number of threads to use for processing.',
)
@click.option(
    '-v',
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
# @click.option('-v', '--verbose', count=True, help="Increase verbosity")
def cli(ctx, sample_dir, out_dir, sample_id, threads, verbose, version):
    if version:
        click.echo(f'g4x-helpers version: {__version__}')

    else:
        ctx.ensure_object(dict)

        # out_dir.mkdir(parents=True, exist_ok=True)

        ctx.obj['sample_dir'] = sample_dir
        ctx.obj['out_dir'] = out_dir
        ctx.obj['threads'] = threads
        ctx.obj['verbose'] = verbose
        ctx.obj['version'] = __version__

        if sample_dir:
            sample = utils.initialize_sample(sample_dir=sample_dir, sample_id=sample_id, n_threads=threads)
            ctx.obj['sample'] = sample
            click.echo(sample)

        # No subcommand given → show help
        if ctx.invoked_subcommand is None:
            click.echo(ctx.get_help())
            ctx.exit()


############################################################
# region resegment
@cli.command(name='resegment', help=RESEG_HELP)
@click.argument(
    'segmentation-mask',
    panel='commands',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to new segmentation mask. Supported file types: .npy, .npz, .geojson.',
)
@click.option(
    '--segmentation-mask-key',
    required=False,
    type=str,
    default=None,
    help='Key in npz/geojson where mask/labels should be taken from (optional).',
)
@click.pass_context
def resegment(ctx, segmentation_mask, segmentation_mask_key):
    main.resegment(
        g4x_out=ctx.obj['sample'],
        segmentation_mask=segmentation_mask,
        out_dir=ctx.obj['out_dir'],
        segmentation_mask_key=segmentation_mask_key,
        n_threads=ctx.obj['threads'],
        verbose=ctx.obj['verbose'],
    )


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
    main.redemux(
        g4x_out=ctx.obj['sample'],
        manifest=manifest,
        out_dir=ctx.obj['out_dir'],
        batch_size=batch_size,
        n_threads=ctx.obj['threads'],
        verbose=ctx.obj['verbose'],
    )


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
    main.update_bin(
        bin_file=bin_file,
        out_dir=ctx.obj['out_dir'],
        metadata=metadata,
        cellid_key=cellid_key,
        cluster_key=cluster_key,
        cluster_color_key=cluster_color_key,
        emb_key=emb_key,
        verbose=ctx.obj['verbose'],
    )


############################################################
# region new_bin
@cli.command(name='new-bin', help=NWBIN_HELP)
@click.pass_context
def new_bin(ctx):
    main.new_bin(
        g4x_out=ctx.obj['sample'],  #
        out_dir=ctx.obj['out_dir'],
        n_threads=ctx.obj['threads'],
        verbose=ctx.obj['verbose'],
    )


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
    main.tar_viewer(
        g4x_out=ctx.obj['sample'],
        viewer_dir=viewer_dir,
        out_dir=ctx.obj['out_dir'],
        verbose=ctx.obj['verbose'],
    )


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
