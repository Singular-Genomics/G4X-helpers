import json
import shutil

import zarr

from .. import constants as c


def setup_viewer(smp, override: bool = True):
    store_path = smp.data_dir / c.FILE_VIEWER_ZARR

    if store_path.exists() and override:
        shutil.rmtree(store_path)

    root_group = zarr.open_group(store_path, mode='w', zarr_version=2)

    meta_path = smp.data_dir / c.FILE_SMP_META
    with open(meta_path, 'r') as f:
        run_metadata = json.load(f)

    (store_path / 'misc').mkdir(parents=True, exist_ok=True)
    shutil.copy(smp.data_dir / c.FILE_SUMMARY, store_path / 'misc' / c.FILE_SUMMARY)

    # TODO re-build full structure
    root_group.attrs['run_metadata'] = {'Sample Information': run_metadata}
    return root_group
