from __future__ import annotations

import inspect


def _supports_param(func, name):
    try:
        return name in inspect.signature(func).parameters
    except (TypeError, ValueError):
        return False


def _normalize_chunks(chunks):
    if chunks == 'auto':
        return None
    return chunks


def _coerce_compressor(compressor):
    if not _is_zarr_v3() or compressor is None:
        return compressor

    try:
        from numcodecs import Blosc as NCBlosc
    except Exception:
        NCBlosc = None

    try:
        import zarr.codecs as zcodecs
    except Exception:
        return compressor

    if NCBlosc is not None and isinstance(compressor, NCBlosc):
        blosc_codec = getattr(zcodecs, 'BloscCodec', None)
        if blosc_codec is not None:
            shuffle_map = {0: 'noshuffle', 1: 'shuffle', 2: 'bitshuffle'}
            shuffle = shuffle_map.get(compressor.shuffle, compressor.shuffle)
            return blosc_codec(cname=compressor.cname, clevel=compressor.clevel, shuffle=shuffle)

    return compressor


def create_array(group, name, data, compressor=None, chunks=None):
    kwargs = {}
    chunks = _normalize_chunks(chunks)
    compressor = _coerce_compressor(compressor)

    if chunks is not None:
        if _supports_param(group.create_array, 'chunks'):
            kwargs['chunks'] = chunks
        elif _supports_param(group.create_array, 'chunk_shape'):
            kwargs['chunk_shape'] = chunks

    if compressor is not None:
        if _supports_param(group.create_array, 'compressors'):
            kwargs['compressors'] = [compressor]
        elif _supports_param(group.create_array, 'compressor'):
            kwargs['compressor'] = compressor

    return group.create_array(name, data=data, **kwargs)


def _is_zarr_v3():
    try:
        import zarr
    except Exception:
        return False
    return str(getattr(zarr, '__version__', '')).startswith('3')


def get_storage_options(chunks=None, compressor=None):
    options = {}
    chunks = _normalize_chunks(chunks)
    compressor = _coerce_compressor(compressor)

    if _is_zarr_v3():
        if chunks is not None:
            options['chunk_shape'] = chunks
        if compressor is not None:
            options['compressors'] = [compressor]
    else:
        if chunks is not None:
            options['chunks'] = chunks
        if compressor is not None:
            options['compressor'] = compressor

    return options
