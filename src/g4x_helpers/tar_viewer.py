import atexit
import json
import shutil
import signal
import sys
import tarfile
from pathlib import Path


def tar_viewer(viewer_dir, out_path):
    print('Checking files.')
    viewer_dir = Path(viewer_dir)
    assert viewer_dir.exists(), f'{viewer_dir} does not appear to exist.'

    bin_path = list(viewer_dir.glob('*.bin'))
    assert len(bin_path) == 1, 'Either no bin file was found in viewer_dir or multiple bin files were found.'
    bin_path = bin_path[0]

    sample_id = bin_path.stem

    ome_tiff_path = viewer_dir / f'{sample_id}.ome.tiff'
    assert ome_tiff_path.is_file(), 'ome.tiff file does not exist.'

    run_meta_path = viewer_dir / f'{sample_id}_run_metadata.json'
    assert run_meta_path.is_file(), 'run_metadata.json file does not exist.'

    tx_path = viewer_dir / f'{sample_id}.tar'
    assert tx_path.is_file(), 'transcript tar file does not exist.'

    # --- H&E paths
    h_and_e_path = viewer_dir / 'h_and_e'
    created_he_dir = False
    if not h_and_e_path.exists():
        h_and_e_path.mkdir(parents=True, exist_ok=True)
        created_he_dir = True

    orig_he = viewer_dir / f'{sample_id}_HE.ome.tiff'
    assert orig_he.is_file(), 'fH&E ome.tiff file does not exist.'
    moved_he = h_and_e_path / f'{sample_id}_HE.ome.tiff'

    moved = False
    restored = False

    def restore_he():
        nonlocal moved, restored, created_he_dir
        if restored:
            return
        try:
            if moved and moved_he.exists():
                if orig_he.exists():
                    orig_he.unlink()
                shutil.move(str(moved_he), str(orig_he))
                print('Restored H&E file to original location.')

            # NEW: remove the h_and_e folder if we created it and it’s empty
            try:
                if created_he_dir and h_and_e_path.exists() and not any(h_and_e_path.iterdir()):
                    h_and_e_path.rmdir()
                    print('Removed empty h_and_e directory.')
            except Exception as e:
                print(f'WARNING: failed to remove h_and_e directory: {e}', file=sys.stderr)

            restored = True
            moved = False
        except Exception as e:
            print(f'WARNING: failed to restore H&E file: {e}', file=sys.stderr)

    atexit.register(restore_he)

    def _handle_signal(sig):
        import signal as _signal

        sig_name = _signal.Signals(sig).name
        print(f'Received {sig_name}; cleaning up…', file=sys.stderr)
        restore_he()
        sys.exit(128 + sig)

    signal.signal(signal.SIGINT, _handle_signal)
    signal.signal(signal.SIGTERM, _handle_signal)

    # Move H&E into subfolder for packaging
    shutil.move(str(orig_he), str(moved_he))
    moved = True

    try:
        print('Making metadata.')
        metadata = {
            'protein_image_src': ome_tiff_path.name,
            'protein_image_data_src': run_meta_path.name,
            'he_images_src': h_and_e_path.name,
            'cell_segmentation_src': bin_path.name,
            'transcript_src': tx_path.name,
        }
        with open(viewer_dir / 'dataset.config.json', 'w') as f:
            json.dump(metadata, f)

        print('Tarring folder.')
        out_tar = Path(out_path)
        if not out_tar.exists():
            out_tar.mkdir(parents=True, exist_ok=True)
        out_tar = out_tar / f'{sample_id}_g4x_viewer.tar'
        with tarfile.open(out_tar, 'w') as tar:
            tar.add(viewer_dir, arcname=viewer_dir.name)

    finally:
        restore_he()
        try:
            atexit.unregister(restore_he)
        except Exception:
            pass
