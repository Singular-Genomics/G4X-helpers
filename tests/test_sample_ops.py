import shutil
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import g4x_helpers as g4x


@pytest.fixture(params=['pr_sample', 'tx_sample'], ids=['pr_sample', 'tx_sample'])
def sample(request):
    return request.getfixturevalue(request.param)


def create_test_manifest(smp: g4x.G4Xsample, out_dir: Path) -> Path:
    txtable = smp.src.TxTable.load()
    manifest_ori = smp.src.Manifest.load()

    top_gene = txtable.group_by('gene_id').agg(pl.len()).sort('len', descending=True)[0]['gene_id'].item()
    replaced = manifest_ori.with_columns(pl.col('gene_name').replace({top_gene: '000'}))

    test_manifest = out_dir / 'test_manifest.csv'
    replaced.write_csv(test_manifest)
    return test_manifest


def create_test_mask(smp: g4x.G4Xsample, out_dir: Path, n_drop: int = 100) -> Path:

    flipped = np.fliplr(smp.src.Segmentation.load())

    remove_ids = np.unique(flipped)[1 : n_drop + 1]
    flipped[np.isin(flipped, remove_ids)] = 0

    data_dict = {'test': flipped}

    test_mask = out_dir / 'test_cell_mask.npz'
    np.savez(test_mask, **data_dict)
    return test_mask


def test_demux_with_provided_manifest(sample):
    smp = sample
    test_manifest = create_test_manifest(smp, smp.smp_dir)

    g4x.ops.demux(smp, manifest=test_manifest, out_dir=None)

    txtable = smp.src.TxTable.load()
    assert not txtable.filter(pl.col('gene_id') == '000').is_empty()


def test_aggregate_with_provided_mask(sample):
    smp = sample
    original_metadata = smp.src.CellMetadata.load()
    original_cellxgene = smp.src.CellxGene.load()
    original_cellxprot = smp.src.CellxProt.load() if smp.src.pr_detected else None

    n_drop = 100
    test_mask = create_test_mask(smp, smp.smp_dir, n_drop=n_drop)

    g4x.ops.aggregate(smp, cell_mask=test_mask, out_dir=None)

    cm = original_metadata.height - smp.src.CellMetadata.load().height
    cg = original_cellxgene.height - smp.src.CellxGene.load().height

    assert cm == cg == n_drop

    if smp.src.pr_detected:
        cp = original_cellxprot.height - smp.src.CellxProt.load().height
        assert cp == n_drop


def test_sc_process(sample):
    smp = sample
    shutil.rmtree(smp.smp_dir / 'single_cell_data', ignore_errors=True)

    g4x.ops.aggregate(smp)
    g4x.ops.sc_process(smp)

    assert smp.src.ClusteringUmap.is_valid == smp.src.Dgex.is_valid == smp.src.AdataH5.is_valid is True


def test_viewer_zarr(sample):
    smp = sample
    shutil.rmtree(sample.src.ViewerZarr.p, ignore_errors=True)

    g4x.ops.viewer_zarr(smp)
    assert smp.src.ViewerZarr.is_valid
