from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from importlib import metadata
from typing import Any, Literal


class RapidsStatus(str, Enum):
    NOT_CHECKED = 'not_checked'
    NOT_INSTALLED = 'not_installed'
    IMPORT_FAILED = 'import_failed'
    NO_GPU = 'no_gpu'
    OK = 'ok'


@dataclass(frozen=True)
class RapidsCheck:
    status: RapidsStatus
    rsc_version: str | None = None
    detail: str | None = None

    @property
    def available(self) -> bool:
        return self.status is RapidsStatus.OK

    def __repr__(self) -> str:
        gap = 11
        lines = [f'{"status":<{gap}} : {self.status.value}']
        if self.rsc_version is not None:
            lines.append(f'{"rsc_version":<{gap}} : {self.rsc_version}')
        lines.append(f'{"detail":<{gap}} : {self.detail}')
        return '\n'.join(lines)


@dataclass(frozen=True)
class ComputeBackend:
    kind: Literal['cpu', 'gpu']
    rapids: RapidsCheck
    preference: Literal['cpu', 'gpu', 'auto']
    cp: Any = None
    rsc: Any = None

    @property
    def use_gpu(self) -> bool:
        return self.kind == 'gpu' and self.rapids.available


@lru_cache(maxsize=1)
def check_rapids() -> RapidsCheck:
    try:
        ver = metadata.version('rapids-singlecell')
    except metadata.PackageNotFoundError:
        return RapidsCheck(
            status=RapidsStatus.NOT_INSTALLED,
            detail="Distribution 'rapids-singlecell' not found. Install with `uv sync --extra rapids`.",
        )

    try:
        import rapids_singlecell as _  # noqa: F401
    except Exception as e:
        return RapidsCheck(
            status=RapidsStatus.IMPORT_FAILED,
            rsc_version=ver,
            detail=f'Import failed: {type(e).__name__}: {e}',
        )

    try:
        import cupy as cp

        n = int(cp.cuda.runtime.getDeviceCount())
        if n <= 0:
            return RapidsCheck(
                status=RapidsStatus.NO_GPU,
                rsc_version=ver,
                detail=f'No CUDA GPUs detected (device count = {n}).',
            )
    except Exception as e:
        return RapidsCheck(
            status=RapidsStatus.NO_GPU,
            rsc_version=ver,
            detail=f'GPU check failed: {type(e).__name__}: {e}',
        )

    return RapidsCheck(
        status=RapidsStatus.OK,
        rsc_version=ver,
        detail=f'{n} CUDA GPU(s) detected.',
    )


@lru_cache(maxsize=None)
def get_backend(which: Literal['cpu', 'gpu', 'auto'] = 'auto') -> ComputeBackend:
    if which == 'cpu':
        rapids = RapidsCheck(
            status=RapidsStatus.NOT_CHECKED,
            detail='CPU backend requested. Rapids availability check was skipped.',
        )
        return ComputeBackend(kind='cpu', rapids=rapids, preference=which)

    rapids = check_rapids()

    if rapids.available:
        import cupy as cp
        import rapids_singlecell as rsc

        return ComputeBackend(kind='gpu', rapids=rapids, preference=which, cp=cp, rsc=rsc)

    if which == 'gpu':
        raise RuntimeError(f'GPU backend requested but unavailable: {rapids.detail}')

    return ComputeBackend(kind='cpu', rapids=rapids, preference=which)
