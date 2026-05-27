import inspect

from .. import __version__
from .. import constants as c
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
    '--backend',
    type=click.Choice(['auto', 'cpu', 'gpu'], case_sensitive=False),
    default='auto',
    show_default=True,
    help='Execution backend.',
)
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
# @click.option('-v', '--verbose', type=int, default=2, count=True, help='Console logging level (0, 1, 2)')
def cli(ctx, backend, verbose, version):
    if version:
        click.echo(f'g4x-helpers: {__version__}')
        ctx.exit()

    # No subcommand and no input given → show help
    if not ctx.invoked_subcommand:
        click.echo(ctx.get_help())
        ctx.exit()

    if ctx.invoked_subcommand:
        ctx.ensure_object(dict)

        ctx.obj['backend'] = backend
        ctx.obj['verbose'] = verbose
        ctx.obj['version'] = __version__


############################################################
# region redemux
name = 'redemux'


@cli.command(name=name, help=hm.REDMX_HELP)
@cli_setup.g4x_data_opt()
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
@cli_setup.branch_opt(name)
@cli_setup.no_downstream_opt(name)
@click.pass_context
def redemux(ctx, g4x_data, manifest, batch_size, branch, no_downstream):
    func_name = inspect.currentframe().f_code.co_name

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import _create_branch
            from ..main_features import redemux as main_redemux

        if branch is not None:
            if branch == 'main':
                out_dir = g4x_data
            else:
                out_dir = _create_branch(g4x_data, branch)
        else:
            out_dir = None

        main_redemux(
            smp_dir=g4x_data,
            out_dir=out_dir,
            manifest=manifest,
            batch_size=batch_size,
            downstream=not no_downstream,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region resegment
name = 'resegment'


@cli.command(name=name, help=hm.RESEG_HELP)
@cli_setup.g4x_data_opt()
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
@cli_setup.branch_opt(name)
@cli_setup.no_downstream_opt(name)
@click.pass_context
def resegment(ctx, g4x_data, cell_labels, labels_key, branch, no_downstream):
    func_name = inspect.currentframe().f_code.co_name

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import _create_branch
            from ..main_features import resegment as main_resegment

        if branch is not None:
            if branch == 'main':
                out_dir = g4x_data
            else:
                out_dir = _create_branch(g4x_data, branch)
        else:
            out_dir = None

        main_resegment(
            smp_dir=g4x_data,
            out_dir=out_dir,
            segmentation_mask=cell_labels,
            mask_key=labels_key,
            downstream=not no_downstream,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


# region viewer
@cli.group(
    context_settings=dict(help_option_names=['-h', '--help']),
    invoke_without_command=True,
    add_help_option=True,
    help='viewer',
)
@click.pass_context
def viewer(ctx):
    pass


############################################################
# region migrate
name = 'migrate'


@cli.command(name=name, help=hm.MIGRT_HELP)
@cli_setup.g4x_data_opt()
@click.option(
    '-o',
    '--out-dir',
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default=None,
    help='Output directory for migration results',
)
@click.option(
    '--roi',
    required=False,
    nargs=4,
    type=int,
    default=None,
    help='Region of interest for migration (x0, y0, x1, y1)',
)
@cli_setup.no_downstream_opt(name)
@click.pass_context
def migrate(ctx, g4x_data, out_dir, roi, no_downstream):
    func_name = inspect.currentframe().f_code.co_name

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import migrate as main_migrate

        main_migrate(
            smp_dir=g4x_data, out_dir=out_dir, roi_coords=roi, downstream=not no_downstream, verbose=ctx.obj['verbose']
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


############################################################
# region validate
@cli.command(name='validate', help=hm.VLDTE_HELP)
@cli_setup.g4x_data_opt()
@click.pass_context
def validate(ctx, g4x_data):
    func_name = inspect.currentframe().f_code.co_name

    try:
        with cli_setup._spinner(f'Initializing {func_name} process...'):
            from ..main_features import validate as main_validate

        main_validate(
            smp_dir=g4x_data,
            verbose=ctx.obj['verbose'],
        )
    except Exception as e:
        cli_setup._fail_message(func_name, e)


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')

############################################################
# region create_zarr
# name = 'create_zarr'


# @cli.command(name=name, help=hm.UDBIN_HELP)
# cli_setup. @g4x_data_opt()
# @click.option(
#     '--metadata',
#     required=True,
#     type=click.Path(exists=True, dir_okay=False),
#     help='Path to metadata table with clustering and/or embedding information. Must contain cell-IDs that match those in the bin file.',
# )
# @click.option(
#     '--cellid-key',
#     default='cell_id',
#     type=str,
#     help='Column name in metadata containing cell-IDs.\n\n If not provided, looks for column named "cell_id"',
# )
# @click.option(
#     '--cluster-key',
#     default=None,
#     type=str,
#     help='Column name in metadata containing cluster IDs.\n\n If not provided, skips updating cluster IDs.',
# )
# @click.option(
#     '--cluster-color-key',
#     default=None,
#     type=str,
#     help='Column name in metadata containing cluster colors.\n\n (format: hex) Only active if cluster_key is updated.\n\n If not provided, assigns colors automatically.',
# )
# @click.option(
#     '--emb-key',
#     default=None,
#     type=str,
#     help='Column name in metadata containing 2D-embedding coordinates.\n\n Parser will look for {emb_key}_1 and {emb_key}_2.\n\n If not provided, skips updating embedding.',
# )
# cli_setup. @in_place_opt(name)
# @click.pass_context
# def update_bin(ctx, g4x_data, metadata, cellid_key, cluster_key, cluster_color_key, emb_key, in_place):
#     func_name = inspect.currentframe().f_code.co_name
#     g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
#     try:
#         with cli_setup._spinner(f'Initializing {func_name} process...'):
#             from ..main_features import update_bin as main_update_bin

#         main_update_bin(
#             g4x_obj=g4x_obj,
#             bin_file=g4x_obj.data_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_segmentation.bin',
#             bin_out=out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_segmentation.bin',
#             out_dir=out_dir,
#             metadata=metadata,
#             cellid_key=cellid_key,
#             cluster_key=cluster_key,
#             cluster_color_key=cluster_color_key,
#             emb_key=emb_key,
#             verbose=ctx.obj['verbose'],
#         )
#     except Exception as e:
#         cli_setup._fail_message(func_name, e)


############################################################
# region new_bin
# name = 'new_bin'


# @cli.command(name=name, help=hm.NWBIN_HELP)
# cli_setup. @g4x_data_opt()
# cli_setup. @in_place_opt(name)
# @click.pass_context
# def new_bin(ctx, g4x_data, in_place):
#     func_name = inspect.currentframe().f_code.co_name
#     g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
#     try:
#         with cli_setup._spinner(f'Initializing {func_name} process...'):
#             from ..main_features import new_bin as main_new_bin

#         main_new_bin(
#             g4x_obj=g4x_obj,
#             out_dir=out_dir,
#             n_threads=ctx.obj['threads'],
#             verbose=ctx.obj['verbose'],
#         )
#     except Exception as e:
#         cli_setup._fail_message(func_name, e)


# ############################################################
# # region tar_viewer
# name = 'tar_viewer'


# @cli.command(name=name, help=hm.TARVW_HELP)
# cli_setup. @g4x_data_opt()
# cli_setup. @in_place_opt(name)
# @click.pass_context
# def tar_viewer(ctx, g4x_data, in_place):
#     func_name = inspect.currentframe().f_code.co_name

#     g4x_obj, out_dir = cli_setup.initialize_sample(data_dir=g4x_data, in_place=in_place, n_threads=ctx.obj['threads'])
#     try:
#         with cli_setup._spinner(f'Initializing {func_name} process...'):
#             from ..main_features import tar_viewer as main_tar_viewer

#         main_tar_viewer(
#             g4x_obj=g4x_obj,
#             out_dir=out_dir,
#             verbose=ctx.obj['verbose'],
#         )
#     except Exception as e:
#         cli_setup._fail_message(func_name, e)
