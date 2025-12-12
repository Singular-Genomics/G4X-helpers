from pathlib import Path

import pytest
import yaml

from g4x_helpers.cli.cli import cli


@pytest.fixture(scope='module')
def workdir(ensure_test_data):
    return ensure_test_data


def load_cli_commands(cfg_path: Path, *, data_dir: Path) -> dict[str, list[str]]:
    cfg = yaml.safe_load(cfg_path.read_text())
    out: dict[str, list[str]] = {}
    for name, parts in cfg['commands'].items():
        # Flatten because YAML anchors can introduce nested lists
        flat: list[str] = []
        for p in parts:
            if isinstance(p, list):
                flat.extend(p)
            else:
                flat.append(p)

        out[name] = [s.format(data_dir=str(data_dir)) for s in flat]
    return out


@pytest.fixture(scope='module')
def cli_commands(workdir):
    """
    Build the CLI invocations with absolute paths so they can be executed
    directly via CliRunner without changing cwd.
    """
    data_dir = Path(workdir)
    return load_cli_commands(Path(__file__).parent / 'cli_commands.yml', data_dir=data_dir)


@pytest.mark.filterwarnings('ignore:Type google\\._upb\\._message\\..*:DeprecationWarning')
@pytest.mark.parametrize(
    'command',
    ['migrate', 'tar_viewer', 'new_bin', 'update_bin', 'resegment', 'redemux'],
)
def test_cli_commands_with_runner(command, cli_commands, runner):
    result = runner.invoke(cli, cli_commands[command], catch_exceptions=False)

    assert result.exit_code == 0, result.output
