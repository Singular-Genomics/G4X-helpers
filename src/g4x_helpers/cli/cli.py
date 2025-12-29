import inspect

from .. import __version__, utils
from . import cli_setup
from . import help_messages as hm

click = cli_setup.click


# region cli
@click.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    invoke_without_command=True,
    add_help_option=True,
    help=hm.CLI_HELP,
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
# @click.option('-v', '--verbose', default=2, count=True, help='Console logging level (0, 1, 2)')
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


def g4x_data_opt():
    return click.argument(
        'g4x-data',
        type=click.Path(exists=True, file_okay=False),
        help='Directory containing G4X-data for a single sample',
        # panel='data i/o',
    )


def in_place_opt(cmd_name: str = ''):
    return click.option(
        '-ip',
        '--in-place',
        is_flag=True,
        help=f'Edit G4X-data in-place if this flag is set.\n\nOtherwise creates a "g4x_helpers/{cmd_name}" folder.',
    )


############################################################
# region resegment
name = 'resegment'


@cli.command(name=name, help=hm.RESEG_HELP)
@g4x_data_opt()
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
@in_place_opt(name)
@click.pass_context
def resegment(ctx, g4x_data, cell_labels, labels_key, in_place):
    func_name = inspect.currentframe().f_code.co_name

    g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import resegment as main_resegment

        main_resegment(
            g4x_obj=g4x_obj,
            out_dir=out_dir,
            cell_labels=cell_labels,
            labels_key=labels_key,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region redemux
name = 'redemux'


@cli.command(name=name, help=hm.REDMX_HELP)
@g4x_data_opt()
@click.option(
    '--manifest',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to manifest for demuxing.\n\n Must contain a "probe_name" column with the format "geneid-sequence-primer"',
)
@click.option(
    '--batch-size',
    default=1_000_000,
    show_default=True,
    type=int,
    help='Number of transcripts to process per batch.',
)
@in_place_opt(name)
@click.pass_context
def redemux(ctx, g4x_data, manifest, batch_size, in_place):
    func_name = inspect.currentframe().f_code.co_name
    g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import redemux as main_redemux

        main_redemux(
            g4x_obj=g4x_obj,
            out_dir=out_dir,
            manifest=manifest,
            batch_size=batch_size,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region update_bin
name = 'update_bin'


@cli.command(name=name, help=hm.UDBIN_HELP)
@g4x_data_opt()
@click.option(
    '--metadata',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Path to metadata table with clustering and/or embedding information. Must contain cell-IDs that match those in the bin file.',
)
@click.option(
    '--cellid-key',
    default='cell_id',
    type=str,
    help='Column name in metadata containing cell-IDs.\n\n If not provided, looks for column named "cell_id"',
)
@click.option(
    '--cluster-key',
    default=None,
    type=str,
    help='Column name in metadata containing cluster IDs.\n\n If not provided, skips updating cluster IDs.',
)
@click.option(
    '--cluster-color-key',
    default=None,
    type=str,
    help='Column name in metadata containing cluster colors.\n\n (format: hex) Only active if cluster_key is updated.\n\n If not provided, assigns colors automatically.',
)
@click.option(
    '--emb-key',
    default=None,
    type=str,
    help='Column name in metadata containing 2D-embedding coordinates.\n\n Parser will look for {emb_key}_1 and {emb_key}_2.\n\n If not provided, skips updating embedding.',
)
@in_place_opt(name)
@click.pass_context
def update_bin(ctx, g4x_data, metadata, cellid_key, cluster_key, cluster_color_key, emb_key, in_place):
    func_name = inspect.currentframe().f_code.co_name
    g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import update_bin as main_update_bin

        main_update_bin(
            g4x_obj=g4x_obj,
            bin_file=g4x_obj.data_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_segmentation.bin',
            bin_out=out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_segmentation.bin',
            out_dir=out_dir,
            metadata=metadata,
            cellid_key=cellid_key,
            cluster_key=cluster_key,
            cluster_color_key=cluster_color_key,
            emb_key=emb_key,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region new_bin
name = 'new_bin'


@cli.command(name=name, help=hm.NWBIN_HELP)
@g4x_data_opt()
@in_place_opt(name)
@click.pass_context
def new_bin(ctx, g4x_data, in_place):
    func_name = inspect.currentframe().f_code.co_name
    g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import new_bin as main_new_bin

        main_new_bin(
            g4x_obj=g4x_obj,
            out_dir=out_dir,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region tar_viewer
name = 'tar_viewer'


@cli.command(name=name, help=hm.TARVW_HELP)
@g4x_data_opt()
@in_place_opt(name)
@click.pass_context
def tar_viewer(ctx, g4x_data, in_place):
    func_name = inspect.currentframe().f_code.co_name

    g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import tar_viewer as main_tar_viewer

        main_tar_viewer(
            g4x_obj=g4x_obj,
            out_dir=out_dir,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region migrate
@cli.command(name='migrate', help=hm.MIGRT_HELP)
@g4x_data_opt()
@click.option(
    '--restore',
    is_flag=True,
    help='Restores an existing migration-backup.',
)
@click.pass_context
def migrate(ctx, g4x_data, restore):
    func_name = inspect.currentframe().f_code.co_name

    g4x_obj, _ = cli_setup.initialize_sample(data_dir=g4x_data, in_place=False, n_threads=ctx.obj['threads'])

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import migrate as main_migrate

        main_migrate(
            g4x_obj=g4x_obj,
            restore=restore,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region validate
@cli.command(name='validate', help=hm.VLDTE_HELP)
@g4x_data_opt()
@click.pass_context
def validate(ctx, g4x_data):
    func_name = inspect.currentframe().f_code.co_name

    g4x_obj, _ = cli_setup.initialize_sample(data_dir=g4x_data, in_place=False, n_threads=ctx.obj['threads'])

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import validate as main_validate

        main_validate(
            g4x_obj=g4x_obj,
            n_threads=ctx.obj['threads'],
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
