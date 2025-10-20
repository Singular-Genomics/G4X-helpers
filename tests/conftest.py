from __future__ import annotations

import os
from pathlib import Path

import pytest
from click.testing import CliRunner

_SAMPLE_DIR_ENV = 'G4X_HELPERS_TEST_SAMPLE_DIR'

_DEFAULT_SAMPLE_DIRS = [
    Path(os.getcwd()).parent / 'test_sample' / 'A01_crop',
]


def _resolve_sample_dir() -> Path | None:
    candidates: list[Path] = []

    env_value = os.environ.get(_SAMPLE_DIR_ENV)
    if env_value:
        candidates.append(Path(env_value).expanduser())

    candidates.extend(_DEFAULT_SAMPLE_DIRS)

    for path in candidates:
        if path and path.exists():
            return path

    return None


@pytest.fixture(scope='session')
def sample_dir() -> Path:
    path = _resolve_sample_dir()
    if path is None:
        pytest.skip(
            f'No test sample directory available. '
            f'Set ${_SAMPLE_DIR_ENV} to point at a valid G4X sample snapshot.'
        )
    return path


@pytest.fixture(scope='session')
def sample_id(sample_dir: Path) -> str:
    return sample_dir.name


@pytest.fixture(scope='session')
def cli_runner() -> CliRunner:
    return CliRunner()


@pytest.fixture()
def cli_out_dir(sample_dir: Path, request) -> Path:
    root = sample_dir / 'pytest_out'
    root.mkdir(exist_ok=True)

    # if hasattr(request.node, 'callspec'):
    #     subdir = request.node.callspec.id
    # else:
    #     subdir = request.node.name

    # safe_subdir = subdir.replace(os.sep, '_').replace('[', '_').replace(']', '')
    # out_dir = root #/ safe_subdir

    # if out_dir.exists():
    #     shutil.rmtree(out_dir)

    # out_dir.mkdir()
    return root
