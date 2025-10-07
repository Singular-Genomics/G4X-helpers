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
click.rich_click.COMMANDS_BEFORE_OPTIONS = True
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


click.rich_click.COMMAND_GROUPS = {
    'g4x-helpers': [
        {'name': 'main', 'commands': ['resegment', 'redemux']},
        # {"name": "utilities", "commands": ["log"]},
    ],
}


# region cli
@click.group(
    context_settings=dict(help_option_names=['-h', '--help']), invoke_without_command=True, add_help_option=True
)
@click.option(
    '--run_base',
    required=False,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help='Path to G4X sample output folder passed to commands.',
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
def cli(ctx, run_base, threads, verbose, version):
    if version:
        v = ver('g4x_helpers')
        click.echo(f'g4x-helpers version: {v}')

    else:
        ctx.ensure_object(dict)

        ctx.obj['run_base'] = run_base
        ctx.obj['threads'] = threads
        ctx.obj['verbose'] = verbose
        # No subcommand given → show help
        if ctx.invoked_subcommand is None:
            click.echo(ctx.get_help())
            # print('')
            ctx.exit()


############################################################
# region resegment
@cli.command(name='resegment')
@click.argument(
    'segmentation_mask',
    panel='commands',
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to new segmentation mask. Supported file types: .npy, .npz, .geojson.',
)
@click.option(
    '--sample_id',  #
    required=False,
    type=str,
    default=None,
    help='Sample ID (optional).',
)
@click.option(
    '--out_dir',
    required=False,
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help='Output directory for new files. Created if it does not exist. If not provided, modifies files in run_base.',
)
@click.option(
    '--segmentation_mask_key',
    required=False,
    type=str,
    default=None,
    help='Key in npz/geojson where mask/labels should be taken from (optional).',
)
@click.pass_context
def resegment(ctx, segmentation_mask, sample_id, out_dir, segmentation_mask_key):
    """Run segmentation re-processing for a G4X sample."""
    from . import segmentation as seg

    click.echo('')
    click.echo(f'Resegmentation (threads={ctx.obj["threads"]})')

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
    logger.info(f'Initializing G4X sample with run_base: {ctx.obj["run_base"]}, sample_id: {sample_id}')
    sample = G4Xoutput(run_base=ctx.obj['run_base'], sample_id=sample_id)
    click.echo(sample)

    # load new segmentation
    # labels = seg.try_load_segmentation(segmentation_mask, sample.shape, segmentation_mask_key)
    # logger.info(f'Loaded new segmentation with shape: {labels.shape}, n_labels: {len(set(labels.flatten()))}')

    # # run intersection with new segmentation
    # seg.intersect_segmentation(
    #     g4x_out=sample, labels=labels, out_dir=out_dir, n_threads=ctx.obj['threads'], logger=logger
    # )
    click.echo('Resegmentation complete.')


############################################################
# region redemux
@cli.command(name='redemux')
@click.argument(
    'manifest',
    panel='commands',
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
@click.option(
    '--out-dir',
    default=None,
    type=click.Path(file_okay=False),
    help='Output directory where new files will be saved. Will be created if it does not exist. '
    'If not provided, the files in run_base will be updated in-place.',
)
@click.pass_context
def redemux(ctx, manifest, batch_size, out_dir):
    """Run demux re-processing for a G4X sample."""
    from . import demux as dmx

    click.echo('')
    click.echo(f'Redemuxing (threads={ctx.obj["threads"]})')

    logger = setup_logger(
        logger_name='redemux',
        stream_logger=True,
        stream_level=verbose_to_log_level(ctx.obj['verbose']),
        file_logger=True,
        file_level=logging.DEBUG,
        out_dir=out_dir,
        clear_handlers=True,
    )

    # initialize G4X sample
    logger.info(f'Initializing G4X sample with run_base: {ctx.obj["run_base"]}')
    sample = G4Xoutput(run_base=ctx.obj['run_base'])
    click.echo(sample)

    # load new segmentation
    dmx.redemux(
        g4x_out=sample,
        manifest=manifest,
        batch_size=batch_size,
        out_dir=out_dir,
        threads=ctx.obj['threads'],
        logger=logger,
    )

    click.echo('Redemuxing complete.')


if __name__ == '__main__':
    cli(prog_name='g4x-helpers')
