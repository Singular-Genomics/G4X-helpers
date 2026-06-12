from pathlib import Path

import pytest
import yaml

from g4x_helpers.cli.main import cli

CLI_COMMANDS_PATH = Path(__file__).parent / 'cli_commands.yml'
CASE_FIELDS = ('args', 'expect_files', 'expect_dirs', 'assert_output_contains')


def load_cli_command_names(cfg_path: Path) -> list[str]:
    cfg = yaml.safe_load(cfg_path.read_text())
    return list(cfg['commands'])


def flatten_cli_parts(parts: list) -> list[str]:
    flat: list[str] = []
    for part in parts:
        if isinstance(part, list):
            flat.extend(flatten_cli_parts(part))
        else:
            flat.append(part)
    return flat


def format_cli_parts(parts: list | None, *, data_dir: Path) -> list[str]:
    if parts is None:
        return []

    return [str(part).format(data_dir=str(data_dir)) for part in flatten_cli_parts(parts)]


def load_cli_commands(cfg_path: Path, *, data_dir: Path) -> dict[str, dict[str, list[str]]]:
    cfg = yaml.safe_load(cfg_path.read_text())
    out: dict[str, dict[str, list[str]]] = {}
    for name, case in cfg['commands'].items():
        if not isinstance(case, dict):
            raise TypeError(f'CLI command case "{name}" must be a mapping')
        if 'args' not in case:
            raise ValueError(f'CLI command case "{name}" must define args')

        out[name] = {field: format_cli_parts(case.get(field), data_dir=data_dir) for field in CASE_FIELDS}
    return out


@pytest.fixture(scope='function')
def cli_commands(workdir):
    """
    Build the CLI invocations with absolute paths so they can be executed
    directly via CliRunner without changing cwd.
    """
    return load_cli_commands(CLI_COMMANDS_PATH, data_dir=workdir)


@pytest.mark.parametrize(
    'command',
    load_cli_command_names(CLI_COMMANDS_PATH),
)
def test_cli_commands_with_runner(command, cli_commands, runner):
    case = cli_commands[command]
    result = runner.invoke(cli, case['args'], catch_exceptions=False)

    assert result.exit_code == 0, result.output

    for expected in case['assert_output_contains']:
        assert expected in result.output

    for expected_file in case['expect_files']:
        assert Path(expected_file).is_file()

    for expected_dir in case['expect_dirs']:
        assert Path(expected_dir).is_dir()
