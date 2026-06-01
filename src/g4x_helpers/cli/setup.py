import traceback
from contextlib import contextmanager

import rich_click as click
from rich.console import Console

from .. import constants

console = Console()


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

click.rich_click.ARGUMENTS_PANEL_TITLE = 'input'

click.rich_click.COMMAND_GROUPS = {
    'g4x-helpers': [
        {'name': 'commands', 'commands': ['viewer', 'redemux', 'resegment', 'migrate', 'validate']},
    ],
    '* viewer': [
        {'name': 'commands', 'commands': ['images', 'cells', 'transcripts']},
    ],
}

click.rich_click.OPTION_GROUPS = {
    'g4x-helpers viewer': [
        {'name': 'input', 'options': ['viewer-zarr']},
    ],
}


@contextmanager
def _spinner(message: str):
    with console.status(message, spinner='dots', spinner_style='red'):
        yield


def _fail_message(func_name, e, trace_back=False):
    click.echo('')
    click.secho(f'Failed {func_name}:', fg='red', err=True, bold=True)
    if trace_back:
        traceback.print_exc()
    raise click.ClickException(f'{type(e).__name__}: {e}')


def initialize_sample(
    data_dir: str, sample_id: str | None = None, in_place: bool = False, n_threads: int = constants.DEFAULT_THREADS
) -> None:
    msg = f'loading G4X-data from [blue]{data_dir}[/blue]'
    with _spinner(msg):
        import glymur

        from ..g4x_output import G4Xoutput

        glymur.set_option('lib.num_threads', n_threads)
        try:
            sample = G4Xoutput(data_dir=data_dir, sample_id=sample_id)
        except Exception as e:
            click.echo('\n')
            click.secho('Failed to load G4X-data:', fg='red', err=True, bold=True)
            raise click.ClickException(f'{e}')

    if in_place:
        out_dir = sample.data_dir
        click.secho('Editing in-place!', fg='blue', bold=True)
    else:
        out_dir = sample.data_dir / 'g4x_helpers'

    return sample, out_dir


def print_k_v(item, value, gap=2):
    value = '<undefined>' if not value else value
    click.secho(f'{item:<{gap}}', dim=True, nl=False)
    click.secho('- ', dim=True, nl=False)
    click.secho(f'{value}', fg='blue', bold=True)


def g4x_data_opt():
    return click.argument(
        'g4x-data',
        type=click.Path(exists=True, file_okay=False),
        help='Directory containing G4X-data for a single sample',
        # panel='data i/o',
    )


help_map = {
    'redemux': 'After demuxing completes, do not create single-cell outputs or initialize viewer files',
    'resegment': 'After aggregation, do not post-process single-cell outputs or initialize viewer files',
    'migrate': 'Only migrate raw data files and metadata, but do not create single-cell output or viewer files',
}


def no_downstream_opt(cmd_name: str = ''):
    return click.option(
        '--no-downstream',
        is_flag=True,
        help=f'{help_map.get(cmd_name, "")}',
    )


def in_place_opt(cmd_name: str = ''):
    return click.option(
        '-ip',
        '--in-place',
        is_flag=True,
        help=f'Edit G4X-data in-place if this flag is set.\n\nOtherwise creates a "g4x_helpers/{cmd_name}" folder.',
    )


def branch_opt(cmd_name: str = ''):
    return click.option(
        '-b',
        '--branch',
        is_flag=False,
        type=str,
        default=None,
        help=(
            f'Branch of processed data to use. If not specified, a branch named '
            f'"g4x-helpers/{cmd_name}" will be created or reused automatically. '
            f'Set to "main" to use the main branch and edit the original data in-place.'
        ),
    )


def out_dir_from_branch(g4x_data, branch):
    from .features.general_group import _create_branch

    if branch is not None:
        if branch == 'main':
            out_dir = g4x_data
        else:
            out_dir = _create_branch(g4x_data, branch)
    else:
        out_dir = None
    return out_dir
