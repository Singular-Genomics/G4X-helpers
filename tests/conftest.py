import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner


@pytest.fixture(scope='session', autouse=True)
def disable_io_cache():
    import g4x_helpers as g4x

    previous = g4x.io.cache()
    g4x.io.cache(use=False)
    yield
    g4x.io.cache(use=previous, clear=True)


def reset_test_data():
    tests_dir = Path('./tests').resolve()
    test_data_dir = tests_dir / 'datasets' / 'test_data'

    if test_data_dir.exists():
        shutil.rmtree(test_data_dir)

    subprocess.run(
        ['bash', str(tests_dir / 'scripts/untar_test_data.sh')],
        check=True,
    )

    return test_data_dir


@pytest.fixture(scope='session')
def ensure_test_data_archive():
    tests_dir = Path('./tests').resolve()
    test_tar = tests_dir / 'datasets' / 'test_data.tar.gz'

    if not test_tar.exists():
        subprocess.run(
            ['bash', str(tests_dir / 'scripts/get_test_data.sh')],
            check=True,
        )


@pytest.fixture(scope='function')
def workdir(ensure_test_data_archive):
    return reset_test_data()


@pytest.fixture(scope='session')
def runner():
    """Shared click CliRunner for CLI tests."""
    return CliRunner()
