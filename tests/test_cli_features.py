import subprocess

import pytest

COMMANDS = [
    ['uv', 'run', 'g4x-helpers', '-v', '2', 'migrate', '.'],
    ['uv', 'run', 'g4x-helpers', '-v', '2', 'tar_viewer', '.'],
    ['uv', 'run', 'g4x-helpers', '-v', '2', 'new_bin', '.'],
    [
        'uv',
        'run',
        'g4x-helpers',
        '-v',
        '2',
        'update_bin',
        '.',
        '--metadata',
        'single_cell_data/clustering_umap.csv.gz',
        '--cluster-key',
        'leiden_0.400',
        '--cellid-key',
        'label',
    ],
    # ["uv", "run", "g4x-helpers", "-v", "2", "resegment", ".", "--cell-labels", "segmentation/segmentation_mask.npz", "--labels-key", "nuclei_exp"],
    # ["uv", "run", "g4x-helpers", "-v", "2", "redemux", ".", "--manifest", "transcript_panel.csv"],
    # ["uv", "run", "g4x-helpers", "-v", "2", "tar_viewer", ".", "--in-place"],
    # ["uv", "run", "g4x-helpers", "-v", "2", "new_bin", ".", "--in-place"],
    # ["uv", "run", "g4x-helpers", "-v", "2", "update_bin", ".", "--metadata", "single_cell_data/clustering_umap.csv.gz", "--cluster-key", "leiden_0.400", "--cellid-key", "label", "--in-place"],
    # ["uv", "run", "g4x-helpers", "-v", "2", "resegment", ".", "--cell-labels", "segmentation/segmentation_mask.npz", "--labels-key", "nuclei_exp", "--in-place"],
    # ["uv", "run", "g4x-helpers", "-v", "2", "redemux", ".", "--manifest", "transcript_panel.csv", "--in-place"],
]


@pytest.fixture(scope='module')
def workdir(ensure_test_data, tmp_path_factory):
    """
    Provide a single mutable copy of the test data for all CLI commands
    to run against, mirroring the original shell script behaviour.
    """
    return ensure_test_data
    # src = ensure_test_data
    # dest = tmp_path_factory.mktemp("test_data_run")
    # shutil.copytree(src, dest / "test_data")
    # return dest / "test_data"


@pytest.mark.parametrize('cmd', COMMANDS)
def test_cli_commands(cmd, workdir):
    subprocess.run(cmd, cwd=workdir, check=True)
