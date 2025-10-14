# import time as _time
# _CLI_T0 = _time.perf_counter()

import inspect
import textwrap

import rich_click as click

from . import __version__, utils
from . import main_features as main

# click.rich_click.THEME = 'modern'
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

click.rich_click.ARGUMENTS_PANEL_TITLE = 'in/out'

click.rich_click.COMMAND_GROUPS = {
    'g4x-helpers': [
        {'name': 'commands', 'commands': ['redemux', 'resegment', 'update_bin', 'new_bin', 'tar_viewer']},
        # {"name": "utilities", "commands": ["log"]},
    ],
}

# click.rich_click.OPTION_GROUPS = {
#     'g4x-helpers': [
#         {
#             'name': 'foo',  #
#             'options': ['--input', '--out-dir'],
#         }
#     ]
# }

CLI_HELP = 'Helper models and post-processing tools for G4X-data\n\ndocs.singulargenomics.com'

RESEG_HELP = (
    'Reprocess G4X-output with a new segmentation\n\n'
    'Takes custom cell segmentation labels as input and re-assigns transcripts and protein signals to those cells.'
)

REDMX_HELP = 'Reprocess G4X-output with a new transcript manifest'

UDBIN_HELP = 'Update existing G4X-viewer .bin file with new metadata'

NWBIN_HELP = 'Create new G4X-viewer .bin file from sample directory'

TARVW_HELP = (
    'Package G4X-viewer folder for distribution\n\n'
    'Creates a .tar archive of the G4X_DATA/g4x_viewer folder for easy upload and sharing.\n\n'
    'Archive file name: {SAMPLE_ID}_g4x_viewer.tar\n\n'
    'If no OUT_DIR is specified via g4x-helpers top-level command, the archive will be created in the G4X_DATA directory.'
)


# region cli
@click.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    invoke_without_command=True,
    add_help_option=True,
    help=CLI_HELP,
)
@click.argument(
    'g4x-data',  #
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help='Directory containing G4X-data for a single sample',
    # panel='input',
)
@click.argument(
    'out-dir',
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True),
    help='Output directory used by subcommands',
    # panel='options [input/output]',
)
@click.option(
    '--sample-id',  #
    required=False,
    type=str,
    help='Sample ID (used for naming outputs)',
)
@click.option(
    '-t',
    '--threads',
    required=False,
    type=int,
    default=utils.DEFAULT_THREADS,
    show_default=True,
    help='Number of threads to use for processing',
)
@click.option(
    '-v',
    '--verbose',
    required=False,
    type=int,
    default=1,
    show_default=True,
    help='Console logging level (0, 1, 2)',
)
@click.option(
    '--version',
    is_flag=True,
    default=False,
    help='Display g4x-helpers version',
)
@click.pass_context
# @click.option('-v', '--verbose', count=True, help="Increase verbosity")
def cli(ctx, g4x_data, out_dir, sample_id, threads, verbose, version):
    # startup_ms = (_time.perf_counter() - _CLI_T0) * 1000
    # click.secho(f'[startup] CLI initialized in {startup_ms:.1f} ms', dim=True)

    if version:
        click.echo(f'g4x-helpers version: {__version__}')

    else:
        ctx.ensure_object(dict)

        ctx.obj['g4x_data'] = g4x_data
        ctx.obj['out_dir'] = out_dir
        ctx.obj['threads'] = threads
        ctx.obj['verbose'] = verbose
        ctx.obj['version'] = __version__

        if g4x_data:
            sample = utils.initialize_sample(g4x_dir=g4x_data, sample_id=sample_id, n_threads=threads)
            ctx.obj['sample'] = sample

        # No subcommand given but G4X_DATA is set → show sample info
        if ctx.invoked_subcommand is None and g4x_data:
            click.secho('\nG4X-helpers successfully initialized with:\n', bold=True, dim=True)
            click.echo(textwrap.indent(repr(sample), prefix='    '))
            click.secho('please provide a subcommand, use -h for help\n', fg='blue')

        # No subcommand and no options given → show help
        elif ctx.invoked_subcommand is None and not g4x_data:
            click.echo(ctx.get_help())
            ctx.exit()

        elif ctx.invoked_subcommand and g4x_data:
            if not out_dir:
                click.secho('!! No output directory specified. Do you want to edit in-place?', fg='blue', bold=True)
                confirm = click.confirm('Confirm in-place editing', default=False)
                if confirm:
                    click.secho('In-place editing confirmed!')
                else:
                    click.secho('Operation cancelled by user.', fg='red')
                    ctx.exit()


############################################################
# region resegment
@cli.command(name='resegment', help=RESEG_HELP)
@click.argument(
    'cell-labels',
    panel='commands',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='File containing cell segmentation labels.\n\nsupported file types: [.npy, .npz, .geojson]',
)
@click.option(
    '--labels-key',
    required=False,
    type=str,
    default=None,
    help='Key/column in npz/geojson where labels should be taken from (optional)',
)
@click.pass_context
def resegment(ctx, cell_labels, labels_key):
    try:
        main.resegment(
            g4x_out=ctx.obj['sample'],
            cell_labels=cell_labels,
            out_dir=ctx.obj['out_dir'],
            labels_key=labels_key,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        func_name = inspect.currentframe().f_code.co_name
        utils._fail_message(func_name, e, trace_back=True)


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
    try:
        main.redemux(
            g4x_out=ctx.obj['sample'],
            manifest=manifest,
            out_dir=ctx.obj['out_dir'],
            batch_size=batch_size,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        func_name = inspect.currentframe().f_code.co_name
        utils._fail_message(func_name, e)


############################################################
# region update_bin
@cli.command(name='update_bin', help=UDBIN_HELP)
@click.option(
    '--metadata',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to metadata table where clustering and/or embedding information will be extracted. Table must contain a header.',
)
@click.option(
    '--bin-file',
    required=False,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to G4X-Viewer segmentation bin file.',
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
    try:
        main.update_bin(
            g4x_out=ctx.obj['sample'],
            bin_file=bin_file,
            out_dir=ctx.obj['out_dir'],
            metadata=metadata,
            cellid_key=cellid_key,
            cluster_key=cluster_key,
            cluster_color_key=cluster_color_key,
            emb_key=emb_key,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        func_name = inspect.currentframe().f_code.co_name
        utils._fail_message(func_name, e)


############################################################
# region new_bin
@cli.command(name='new_bin', help=NWBIN_HELP)
@click.pass_context
def new_bin(ctx):
    try:
        main.new_bin(
            g4x_out=ctx.obj['sample'],  #
            out_dir=ctx.obj['out_dir'],
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        func_name = inspect.currentframe().f_code.co_name
        utils._fail_message(func_name, e)


############################################################
# region tar_viewer
@cli.command(name='tar_viewer', help=TARVW_HELP)
@click.option(
    '--viewer-dir',
    required=False,
    type=click.Path(exists=True, file_okay=False),
    help='(optional) Path to G4X-viewer folder. If set, will tar specified folder instead of the one supplied by G4X_DIR.',
)
@click.pass_context
def tar_viewer(ctx, viewer_dir):
    try:
        main.tar_viewer(
            g4x_out=ctx.obj['sample'],  #
            viewer_dir=viewer_dir,
            out_dir=ctx.obj['out_dir'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        func_name = inspect.currentframe().f_code.co_name
        utils._fail_message(func_name, e)


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
