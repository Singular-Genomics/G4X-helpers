import json
import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from g4x_helpers import G4Xoutput


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
