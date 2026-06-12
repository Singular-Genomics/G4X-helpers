import json
import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from g4x_helpers import G4Xoutput

TESTS_DIR = Path('./tests').resolve()
TEST_DATA_DIR = TESTS_DIR / 'datasets' / 'test_data'
TEST_DATA_ARCHIVE = TESTS_DIR / 'datasets' / 'test_data.tar.gz'


@pytest.fixture(scope='session', autouse=True)
def disable_io_cache():
    import g4x_helpers as g4x

    previous = g4x.io.cache()
    g4x.io.cache(use=False)
    yield
    g4x.io.cache(use=previous, clear=True)


def remove_test_data():
    if TEST_DATA_DIR.exists():
        shutil.rmtree(TEST_DATA_DIR)


def remove_test_data_source():
    if TEST_DATA_ARCHIVE.exists():
        shutil.rmtree(TEST_DATA_ARCHIVE.parent)


def reset_test_data():
    remove_test_data()

    subprocess.run(
        ['bash', str(TESTS_DIR / 'scripts/untar_test_data.sh')],
        check=True,
    )

    return TEST_DATA_DIR


@pytest.fixture(scope='session')
def ensure_test_data_archive():
    if not TEST_DATA_ARCHIVE.exists():
        subprocess.run(
            ['bash', str(TESTS_DIR / 'scripts/get_test_data.sh')],
            check=True,
        )
    yield
    remove_test_data_source()


@pytest.fixture(scope='function')
def workdir(ensure_test_data_archive):
    yield reset_test_data()
    remove_test_data()


@pytest.fixture(scope='function')
def pr_sample(workdir):
    return G4Xoutput(workdir)


@pytest.fixture(scope='function')
def tx_sample(pr_sample):
    smp = pr_sample

    new_meta = smp.smp_meta.copy()
    new_meta['protein_panel'] = None
    with open(smp.src.SampleG4X.p, 'w') as f:
        json.dump(new_meta, f)

    shutil.rmtree(smp.smp_dir / 'protein')
    smp.src.ProteinPanel.p.unlink()

    return G4Xoutput(smp.smp_dir)


@pytest.fixture(scope='session')
def runner():
    """Shared click CliRunner for CLI tests."""
    return CliRunner()
