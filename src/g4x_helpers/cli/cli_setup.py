import traceback
from contextlib import contextmanager

import rich_click as click
from rich.console import Console

from .. import utils

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
    data_dir: str, sample_id: str | None = None, in_place: bool = False, n_threads: int = utils.DEFAULT_THREADS
) -> None:
    msg = f'loading G4X-data from [blue]{data_dir}[/blue]'
    with _spinner(msg):
        import glymur

        from ..models import G4Xoutput

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
