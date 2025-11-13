import inspect

import rich_click as click

from . import __version__, utils

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
click.rich_click.COMMANDS_BEFORE_OPTIONS = True

# click.rich_click.ARGUMENTS_PANEL_TITLE = 'in/out'

click.rich_click.COMMAND_GROUPS = {
    'g4x-helpers': [
        {'name': 'commands', 'commands': ['redemux', 'resegment', 'update_bin', 'new_bin', 'tar_viewer']},
        # {"name": "utilities", "commands": ["log"]},
    ],
}


# click.rich_click.OPTION_GROUPS = {
#     'g4x-helpers': [
#         {
#             'name': 'in/out',  #
#             'options': ['--input', '--output'],
#         },
#         {
#             'name': 'options',  #
#             'options': [
#                 '--sample-id',
#                 '--threads',
#                 '--verbose',
#                 '--version',
#                 '--help',
#             ],
#         },
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
    'Creates a .tar archive of the "g4x_viewer" directory for easy upload and sharing.\n\n'
    'Archive file name: {SAMPLE_ID}_g4x_viewer.tar\n\n'
)

VALIDATE_HELP = 'Validate G4X-data directory structure and files'


def io_opts(*, takes_input=False, writes_output=True):
    options = []

    def decorator(f):
        if takes_input:
            options.append(
                click.option(
                    '-i',
                    '--g4x-data',
                    required=True,
                    type=click.Path(exists=True, file_okay=False),
                    help='Directory containing G4X-data for a single sample',
                    panel='data i/o',
                )
            )
            options.append(
                click.option(
                    '-s',
                    '--sample-id',  #
                    required=False,
                    type=str,
                    help='Sample ID (used for naming outputs)',
                    panel='data i/o',
                )
            )
        if writes_output:
            options.append(
                click.option(
                    '-o',
                    '--output',
                    required=False,
                    type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True),
                    help='Output directory used. Will edit G4X-data in-place if not provided.',
                    panel='data i/o',
                )
            )

        for opt in reversed(options):
            f = opt(f)
        return f

    return decorator


# region cli
@click.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    invoke_without_command=True,
    add_help_option=True,
    help=CLI_HELP,
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
def cli(ctx, threads, verbose, version):
    if version:
        click.echo(f'g4x-helpers: {__version__}')
        ctx.exit()

    # No subcommand and no input given → show help
    if not ctx.invoked_subcommand:
        click.echo(ctx.get_help())
        ctx.exit()

    if ctx.invoked_subcommand:
        ctx.ensure_object(dict)
        ctx.obj['threads'] = threads
        ctx.obj['verbose'] = verbose
        ctx.obj['version'] = __version__


############################################################
# region resegment
@cli.command(name='resegment', help=RESEG_HELP)
@io_opts(takes_input=True, writes_output=True)
@click.option(
    '--cell-labels',
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
def cli_resegment(ctx, g4x_data, output, sample_id, cell_labels, labels_key):
    func_name = inspect.currentframe().f_code.co_name.removeprefix('cli_')

    g4x_obj, output = utils.initialize_sample(
        data_dir=g4x_data, output=output, sample_id=sample_id, n_threads=ctx.obj['threads']
    )

    try:
        with utils._spinner(f'Initializing {func_name} process...'):
            from .main_features import resegment

        resegment(
            g4x_obj=g4x_obj,
            out_dir=output,
            cell_labels=cell_labels,
            labels_key=labels_key,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        utils._fail_message(func_name, e, trace_back=False)


############################################################
# region redemux
@cli.command(name='redemux', help=REDMX_HELP)
@io_opts(takes_input=True, writes_output=True)
@click.option(
    '--manifest',
    required=True,
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
def cli_redemux(ctx, g4x_data, output, sample_id, manifest, batch_size):
    func_name = inspect.currentframe().f_code.co_name.removeprefix('cli_')
    g4x_obj, output = utils.initialize_sample(
        data_dir=g4x_data, output=output, sample_id=sample_id, n_threads=ctx.obj['threads']
    )
    try:
        with utils._spinner(f'Initializing {func_name} process...'):
            from .main_features import redemux

        redemux(
            g4x_obj=g4x_obj,
            out_dir=output,
            manifest=manifest,
            batch_size=batch_size,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        utils._fail_message(func_name, e)


############################################################
# region update_bin
@cli.command(name='update_bin', help=UDBIN_HELP)
@io_opts(takes_input=True, writes_output=True)
@click.option(
    '--metadata',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to metadata table where clustering and/or embedding information will be extracted. Table must contain a header.',
)
@click.option(
    '--bin-file',
    type=click.Path(exists=True, dir_okay=False),
    help='Path to a specific .bin file to update.',
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
def cli_update_bin(
    ctx, g4x_data, output, sample_id, metadata, bin_file, cellid_key, cluster_key, cluster_color_key, emb_key
):
    func_name = inspect.currentframe().f_code.co_name.removeprefix('cli_')
    g4x_obj, output = utils.initialize_sample(
        data_dir=g4x_data, output=output, sample_id=sample_id, n_threads=ctx.obj['threads']
    )
    try:
        with utils._spinner(f'Initializing {func_name} process...'):
            from .main_features import update_bin

        update_bin(
            g4x_obj=g4x_obj,
            out_dir=output,
            bin_file=bin_file,
            metadata=metadata,
            cellid_key=cellid_key,
            cluster_key=cluster_key,
            cluster_color_key=cluster_color_key,
            emb_key=emb_key,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        utils._fail_message(func_name, e)


############################################################
# region new_bin
@cli.command(name='new_bin', help=NWBIN_HELP)
@io_opts(takes_input=True, writes_output=True)
@click.pass_context
def cli_new_bin(ctx, g4x_data, output, sample_id):
    func_name = inspect.currentframe().f_code.co_name.removeprefix('cli_')
    g4x_obj, output = utils.initialize_sample(
        data_dir=g4x_data, output=output, sample_id=sample_id, n_threads=ctx.obj['threads']
    )
    try:
        with utils._spinner(f'Initializing {func_name} process...'):
            from .main_features import new_bin

        new_bin(
            g4x_obj=g4x_obj,  #
            out_dir=output,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        utils._fail_message(func_name, e, trace_back=True)


############################################################
# region tar_viewer
@cli.command(name='tar_viewer', help=TARVW_HELP)
@io_opts(takes_input=True, writes_output=True)
@click.option(
    '--viewer-dir',
    required=False,
    type=click.Path(exists=True, file_okay=False),
    help='(optional) Path to G4X-viewer folder. If set, will archive the specified folder instead of the one supplied by G4X-data.',
)
@click.pass_context
def cli_tar_viewer(ctx, g4x_data, output, sample_id, viewer_dir):
    func_name = inspect.currentframe().f_code.co_name.removeprefix('cli_')
    g4x_obj, output = utils.initialize_sample(
        data_dir=g4x_data, output=output, sample_id=sample_id, n_threads=ctx.obj['threads']
    )
    try:
        with utils._spinner(f'Initializing {func_name} process...'):
            from .main_features import tar_viewer

        tar_viewer(
            g4x_obj=g4x_obj,
            out_dir=output,
            viewer_dir=viewer_dir,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        utils._fail_message(func_name, e)


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
