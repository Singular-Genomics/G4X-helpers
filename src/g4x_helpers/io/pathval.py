from pathlib import Path


def _normalize_path(path_str, *, resolve: bool = False) -> Path:
    """Internal helper: expand ~ and optionally resolve."""
    path = Path(path_str).expanduser()
    if resolve:
        path = path.resolve()
    return path


def validate_file_path(path, *, must_exist: bool = True, resolve: bool = False) -> Path:
    """
    Validate that a path is a file.
    """
    path = _normalize_path(path, resolve=resolve)

    if must_exist and not path.exists():
        raise FileNotFoundError(f'File does not exist: {path}')

    if path.exists() and not path.is_file():
        raise ValueError(f'Expected file, got directory: {path}')

    return path


def validate_dir_path(path, *, must_exist: bool = True, resolve: bool = False) -> Path:
    """
    Validate that a path is a directory.
    """
    path = _normalize_path(path, resolve=resolve)

    if must_exist and not path.exists():
        raise FileNotFoundError(f'Directory does not exist: {path}')

    if path.exists() and not path.is_dir():
        raise ValueError(f'Expected directory, got file: {path}')

    return path


def ensure_dir(path, *, resolve: bool = False) -> Path:
    """
    Ensure a directory exists.
    """
    path = _normalize_path(path, resolve=resolve)

    if path.exists():
        if not path.is_dir():
            raise ValueError(f'Expected directory, got file: {path}')
        return path

    path.mkdir(parents=True, exist_ok=True)
    return path


def ensure_parent_dir(path, *, resolve: bool = False) -> Path:
    """
    Ensure the parent directory of a file path exists.
    """
    path = _normalize_path(path, resolve=resolve)

    parent = path.parent

    if parent.exists():
        if not parent.is_dir():
            raise ValueError(f'Parent is not a directory: {parent}')
        return path

    parent.mkdir(parents=True, exist_ok=True)
    return path
