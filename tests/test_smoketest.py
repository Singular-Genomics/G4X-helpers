import pytest


def test_cli_help_option(runner):
    """Exercise the top-level CLI to ensure the package is wired up."""
    from g4x_helpers.cli.main import cli

    result = runner.invoke(cli, ['--help'])

    assert result.exit_code == 0, result.output
    assert 'docs.singulargenomics.com' in result.output

    result = runner.invoke(cli, ['--version'])
    assert result.exit_code == 0, result.output
    assert result.output.strip().startswith('g4x-helpers:')


def test_imports():
    try:
        import g4x_helpers as g4x
    except ImportError as e:
        raise AssertionError('Failed to import g4x_helpers') from e

    for name in g4x.__all__:
        try:
            getattr(g4x, name)
            print(f"Imported module '{name}' successfully.")
        except AttributeError as e:
            raise AssertionError(f'Failed to import module: {name}') from e


def test_cli_imports():
    try:
        from g4x_helpers import cli
    except ImportError as e:
        raise AssertionError('Failed to import g4x_helpers') from e

    for name in cli.__all__:
        try:
            getattr(cli, name)
            print(f"Imported module '{name}' successfully.")
        except AttributeError as e:
            raise AssertionError(f'Failed to import module: {name}') from e
