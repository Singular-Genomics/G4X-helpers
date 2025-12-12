import pytest

from g4x_helpers.cli.cli import cli

features = ['resegment', 'redemux', 'new_bin', 'update_bin', 'migrate', 'tar_viewer']


def test_cli_help_option(runner):
    """Exercise the top-level CLI to ensure the package is wired up."""
    result = runner.invoke(cli, ['--help'])

    assert result.exit_code == 0, result.output
    assert 'docs.singulargenomics.com' in result.output
    for feature in features:
        assert feature in result.output


def test_cli_version_option(runner):
    result = runner.invoke(cli, ['--version'])

    assert result.exit_code == 0, result.output
    assert result.output.strip().startswith('g4x-helpers:')


def test_import_main_features():
    try:
        from g4x_helpers import main_features as mf
    except ImportError as e:
        raise AssertionError('Failed to import main_features') from e

    for name in features:
        try:
            getattr(mf, name)
            print(f"Imported feature '{name}' successfully.")
        except AttributeError as e:
            raise AssertionError(f'Failed to import feature: {name}') from e


@pytest.mark.filterwarnings('ignore:Type google\\._upb\\._message\\..*:DeprecationWarning')
def test_import_G4Xoutput():
    try:
        from g4x_helpers import G4Xoutput as G4Xoutput
    except ImportError as e:
        raise AssertionError('Failed to import G4Xoutput class') from e
