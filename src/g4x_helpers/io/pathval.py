from pathlib import Path


def _ingest_path(path_str, *, must_exist: bool = True, resolve: bool = False) -> Path:
    path = Path(path_str).expanduser()
    if resolve:
        path = path.resolve()

    if must_exist and not path.exists():
        raise FileNotFoundError(f'Path does not exist: {path}')

    return path


def validate_file_path(path, *, resolve: bool = False) -> Path:
    """
    Validate that a path is a file.
    """
    path = _ingest_path(path, must_exist=True, resolve=resolve)

    if not path.is_file():
        raise ValueError(f'Expected file, got directory: {path}')

    return path


def validate_dir_path(path, *, resolve: bool = False) -> Path:
    """
    Validate that a path is a directory.
    """
    path = _ingest_path(path, must_exist=True, resolve=resolve)

    if not path.is_dir():
        raise ValueError(f'Expected directory, got file: {path}')

    return path


def validate_file_parent(path, *, resolve: bool = False) -> Path:
    """
    Validate that the parent directory of a file path exists.
    """
    path = _ingest_path(path, must_exist=False, resolve=resolve)
    _ = validate_dir_path(path.parent, resolve=resolve)
    return path


def validate_dir_parent(path, *, resolve: bool = False) -> Path:
    """
    Validate that the parent directory of a directory path exists.
    """
    path = _ingest_path(path, must_exist=False, resolve=resolve)
    _ = validate_dir_path(path.parent, resolve=resolve)
    return path


def ensure_dir(path, *, resolve: bool = False) -> Path:
    """
    Ensure a directory exists.
    """
    path = _ingest_path(path, must_exist=False, resolve=resolve)

    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)

    return validate_dir_path(path, resolve=resolve)


def ensure_parent_dir(path, *, resolve: bool = False) -> Path:
    """
    Ensure the parent directory of a file path exists.
    """
    path = _ingest_path(path, must_exist=False, resolve=resolve)

    if not path.parent.exists():
        ensure_dir(path.parent, resolve=resolve)

    return path
