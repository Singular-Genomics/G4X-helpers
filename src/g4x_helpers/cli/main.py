import inspect
from pathlib import Path

from .. import __version__
from .. import constants as c
from . import help_messages as hm
from . import setup
from .features import general_group as gfeats
from .features import viewer_group as vfeats
from .setup import click


# region cli
@click.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    invoke_without_command=True,
    add_help_option=True,
    help=hm.CLI_HELP,
)
# @click.option(
#     '--backend',
#     type=click.Choice(['auto', 'cpu', 'gpu'], case_sensitive=False),
#     default='auto',
#     show_default=True,
#     help='Execution backend for GPU-accelerated operations',
# )
@click.option(
    '-v',
    '--verbose',
    type=click.Choice([0, 1, 2], case_sensitive=False),
    default=1,
    required=False,
    show_default=True,
    help='Console logging level',
)
@click.option(
    '--version',
    is_flag=True,
    default=False,
    help='Display g4x-helpers version',
)
@click.pass_context
def cli(ctx, verbose, version):
    if version:
        click.echo(f'g4x-helpers: {__version__}')
        ctx.exit()

    # No subcommand and no input given → show help
    if not ctx.invoked_subcommand:
        click.echo(ctx.get_help())
        ctx.exit()

    if ctx.invoked_subcommand:
        ctx.ensure_object(dict)

        ctx.obj['backend'] = 'auto'  # backend
        ctx.obj['verbose'] = verbose
        ctx.obj['version'] = __version__


############################################################
# region redemux
name = 'redemux'


@cli.command(name=name, help=hm.REDEMUX_HELP)
@setup.g4x_data_opt()
@click.option(
    '--manifest',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to manifest for demuxing.\n\n Must contain a "probe" column with the format "geneid-sequence-primer"',
)
@click.option(
    '--batch-size',
    default=c.DEFAULT_BATCH_SIZE,
    show_default=True,
    type=int,
    help='Number of transcripts to process per batch.',
)
@setup.branch_opt(name)
@setup.no_downstream_opt(name)
@click.pass_context
def redemux(ctx, g4x_data, manifest, batch_size, branch, no_downstream):
    func_name = inspect.currentframe().f_code.co_name

    try:
        # with setup._spinner(f'Running {func_name} process...'):
        out_dir = setup.out_dir_from_branch(g4x_data, branch)

        gfeats.redemux(
            smp_dir=g4x_data,
            out_dir=out_dir,
            manifest=manifest,
            batch_size=batch_size,
            downstream=not no_downstream,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        setup._fail_message(func_name, e)


############################################################
# region resegment
name = 'resegment'


@cli.command(name=name, help=hm.RESEGMENT_HELP)
@setup.g4x_data_opt()
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
    help='Key/column in npz/geojson where labels should be taken from (optional, but required for .npz with multiple arrays)',
)
@setup.branch_opt(name)
@setup.no_downstream_opt(name)
@click.pass_context
def resegment(ctx, g4x_data, cell_labels, labels_key, branch, no_downstream):
    func_name = inspect.currentframe().f_code.co_name

    try:
        # with setup._spinner(f'Initializing {func_name} process...'):
        out_dir = setup.out_dir_from_branch(g4x_data, branch)

        gfeats.resegment(
            smp_dir=g4x_data,
            out_dir=out_dir,
            segmentation_mask=cell_labels,
            mask_key=labels_key,
            downstream=not no_downstream,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        setup._fail_message(func_name, e)


############################################################
# region migrate
name = 'migrate'


@cli.command(name=name, help=hm.MIGRATE_HELP)
@setup.g4x_data_opt()
@click.option(
    '-o',
    '--out-dir',
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default=None,
    help='Output directory for migration results',
)
@click.option(
    '-c',
    '--check',
    required=False,
    is_flag=True,
    help='Check if the sample is migratable; without moving data',
)
@click.option(
    '--roi',
    required=False,
    nargs=4,
    type=int,
    default=None,
    help='Region of interest for migration (x0, y0, x1, y1)',
)
@setup.no_downstream_opt(name)
@click.pass_context
def migrate(ctx, g4x_data, out_dir, check, roi, no_downstream):
    func_name = inspect.currentframe().f_code.co_name

    try:
        # with setup._spinner(f'Initializing {func_name} process...'):
        #     from .features.general_group import migrate

        if check:
            gfeats.migrate_check(smp_dir=g4x_data)
            return

        gfeats.migrate(
            smp_dir=g4x_data,
            out_dir=out_dir,
            roi_coords=roi,
            downstream=not no_downstream,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        setup._fail_message(func_name, e)


############################################################
# region validate
@cli.command(name='validate', help=hm.VALIDATE_HELP)
@setup.g4x_data_opt()
@click.pass_context
def validate(ctx, g4x_data):
    func_name = inspect.currentframe().f_code.co_name

    try:
        # with setup._spinner(f'Initializing {func_name} process...'):
        #     from .features.general_group import validate

        gfeats.validate(
            smp_dir=g4x_data,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        setup._fail_message(func_name, e)


# region viewer
@cli.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    add_help_option=True,
    help=hm.VIEWER_HELP,
)
@click.argument(
    'viewer-zarr',
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help='Path to a g4x-viewer.zarr',
)
@click.pass_context
def viewer(ctx, viewer_zarr):
    ctx.obj = {'viewer_zarr': viewer_zarr}


@viewer.command(name='images', help='Modify image metadata in a G4X-viewer zarr store')
@click.option(
    '--export-metadata',
    type=click.Choice(['auto', 'cpu', 'gpu'], case_sensitive=False),
    default='auto',
    show_default=True,
    help='Execution backend.',
)
@click.pass_context
def images(ctx):
    pass


@viewer.command(
    name='cells',
    help='Modify cell metadata in a G4X-viewer zarr store\n\n--import and --export options cannot be used together',
)
@click.option(
    '--import-metadata',
    type=click.Path(exists=True, dir_okay=False),
    help='CSV file containing cell metadata to import.',
)
@click.option(
    '--export-metadata',
    type=click.Path(exists=False, writable=True, dir_okay=False),
    show_default=True,
    help='Output CSV file for exported cell metadata.',
)
@click.option(
    '--segmentation',
    type=str,
    default='g4x_default_segmentation',
    required=False,
    show_default=False,
    help='Only required if multiple segmentations are available.',
)
@click.pass_context
def cells(ctx, import_metadata, export_metadata, segmentation):
    func_name = 'viewer/' + inspect.currentframe().f_code.co_name

    if not import_metadata and not export_metadata:
        click.echo('Please provide one of --import-metadata or --export-metadata options')
        ctx.exit(0)
    
    try:
        vfeats.cell_metadata(ctx.obj['viewer_zarr'], import_metadata, export_metadata, segmentation)

    except Exception as e:
        setup._fail_message(func_name, e)


@viewer.command(name='transcripts', help='Modify transcript metadata in a G4X-viewer zarr store')
@click.option(
    '--export-metadata',
    type=click.Choice(['auto', 'cpu', 'gpu'], case_sensitive=False),
    default='auto',
    show_default=True,
    help='Execution backend.',
)
@click.pass_context
def transcripts(ctx):
    pass


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
