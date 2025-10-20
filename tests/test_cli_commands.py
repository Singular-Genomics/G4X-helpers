from __future__ import annotations

from pathlib import Path

import pytest

from g4x_helpers.cli import cli

PROTEIN_MANIFEST = Path(__file__).parent / 'protein_primer_seqs_manifest.csv'

COMMAND_CASES = [
    pytest.param(
        'tar_viewer',
        (),
        ('g4x_viewer/{sample_id}.tar',),
        id='tar_viewer',
    ),
    pytest.param(
        'resegment',
        (
            '{sample_dir}/segmentation/segmentation_mask.npz',
            '--labels-key',
            'nuclei_exp',
        ),
        ('g4x_viewer/{sample_id}.bin',),
        id='resegment',
    ),
    pytest.param(
        'new_bin',
        (),
        ('g4x_viewer/{sample_id}.bin',),
        id='new_bin',
    ),
    pytest.param(
        'redemux',
        (str(PROTEIN_MANIFEST),),
        (
            'g4x_viewer/{sample_id}.bin',
            'g4x_viewer/{sample_id}.tar',
        ),
        id='redemux',
        marks=pytest.mark.xfail(
            reason=(
                'Cropped sample dataset triggers IndexError during redemux '
                'segmentation intersection.'
            ),
            strict=False,
        ),
    ),
    pytest.param(
        'update_bin',
        (
            '--metadata',
            '{sample_dir}/single_cell_data/clustering_umap.csv.gz',
            '--emb-key',
            'X_umap_0.300_0.800',
        ),
        ('g4x_viewer/{sample_id}.bin',),
        id='update_bin',
    ),
]


@pytest.mark.parametrize('command,extra_args,expected_relpaths', COMMAND_CASES)
def test_cli_commands(cli_runner, sample_dir, sample_id, cli_out_dir, command, extra_args, expected_relpaths):
    formatted_args = [
        '-t',
        '1',
        '--sample-id',
        sample_id,
        str(sample_dir),
        str(cli_out_dir),
        command,
    ]

    formatted_args.extend(
        [
            str(arg).format(sample_dir=sample_dir, sample_id=sample_id)
            for arg in extra_args
        ]
    )

    result = cli_runner.invoke(cli, formatted_args, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    command_out_dir = Path(cli_out_dir) / command
    assert command_out_dir.exists(), f'Expected command output directory: {command_out_dir}'

    for rel_path in expected_relpaths:
        expected_path = command_out_dir / rel_path.format(sample_id=sample_id)
        assert expected_path.exists(), f'Missing expected output: {expected_path}'
