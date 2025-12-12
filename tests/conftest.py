import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner


@pytest.fixture(scope='session', autouse=False)
def ensure_test_data():
    """
    Make sure the test data archive is present and extracted before CLI tests run.
    Downloads once (if missing) via get_test_data.sh, then extracts with untar_test_data.sh.
    """
    tests_dir = Path('./tests').resolve()
    test_tar = tests_dir / 'datasets' / 'test_data.tar'

    if not test_tar.exists():
        subprocess.run(['bash', str(tests_dir / 'scripts/get_test_data.sh')], check=True)

    subprocess.run(['bash', str(tests_dir / 'scripts/untar_test_data.sh')], check=True)
    return tests_dir / 'datasets' / 'test_data'


@pytest.fixture(scope='session')
def runner():
    """Shared click CliRunner for CLI tests."""
    return CliRunner()
