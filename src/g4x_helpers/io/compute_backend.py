from __future__ import annotations

import warnings
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

    def __repr__(self):
        gap = 11
        repr_str = f'{"status":<{gap}} : {self.status[0:]}\n'
        if self.rsc_version is not None:
            repr_str += f'{"rsc_version":<{gap}} : {self.rsc_version}\n'
        repr_str += f'{"detail":<{gap}} : {self.detail}'

        return repr_str


@dataclass(frozen=True)
class ComputeBackend:
    kind: str
    rapids: RapidsCheck
    preference: str
    cp: Any = None
    rsc: Any = None

    @property
    def use_gpu(self) -> bool:
        return self.kind == 'gpu' and self.rapids.available


@lru_cache(maxsize=1)
def get_backend(preference: Literal['cpu', 'gpu', 'auto'] = 'auto') -> ComputeBackend:

    if preference == 'cpu':
        rapids = RapidsCheck(
            status=RapidsStatus.NOT_CHECKED,
            detail='CPU backend requested. Rapids availability check was skipped.',
        )

    if preference in ['gpu', 'auto']:
        rapids = check_rapids()
        if rapids.available:
            import cupy as cp
            import rapids_singlecell as rsc

            return ComputeBackend(kind='gpu', rapids=rapids, preference=preference, cp=cp, rsc=rsc)
        else:
            warnings.warn('GPU backend requested but unavailable, falling back to CPU')

    return ComputeBackend(kind='cpu', rapids=rapids, preference=preference)


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


# previous implementation
# class RapidsStatus(str, Enum):
#     s0 = 'not_installed'  # metadata missing
#     s1 = 'import_failed'  # dist present but module import fails
#     s2 = 'no_gpu'  # import ok but no usable CUDA GPU
#     s3 = 'ok'  # ready to use


# @dataclass(frozen=True)
# class RapidsCheck:
#     status: RapidsStatus
#     rsc_version: Optional[str] = None
#     detail: Optional[str] = None

#     def __repr__(self):
#         gap = 11
#         repr_str = f'{"status":<{gap}} : {self.status[0:]}\n'
#         if self.rsc_version is not None:
#             repr_str += f'{"rsc_version":<{gap}} : {self.rsc_version}\n'
#         repr_str += f'{"detail":<{gap}} : {self.detail}'

#         return repr_str


# def check_rapids() -> bool:
#     # --- Stage 1: metadata presence (fast) ---
#     def report_current():
#         res = RapidsCheck(status=status, rsc_version=ver, detail=detail)
#         print('--- Checking for availability of rapids-singlecell for GPU acceleration ---')
#         print(res)

#     ver = None
#     n = 0
#     try:
#         ver = metadata.version('rapids-singlecell')
#     except metadata.PackageNotFoundError:
#         status = RapidsStatus.s0
#         detail = "Distribution 'rapids-singlecell' not found. Install with `uv sync --extra rapids`. (linux only)"
#         report_current()
#         return False

#     # --- Stage 2: import module (realistic) ---
#     try:
#         import rapids_singlecell as _  # noqa: F401
#     except Exception as e:
#         status = RapidsStatus.s1
#         detail = f'Import failed: {type(e).__name__}: {e}'
#         report_current()
#         return False

#     # --- Stage 3: verify a CUDA GPU is usable ---
#     # You can choose your preferred “GPU probe”. Cupy is a common, lightweight-ish probe.
#     try:
#         import cupy as cp  # many RAPIDS stacks already include it; if not, this will fail clearly

#         n = int(cp.cuda.runtime.getDeviceCount())
#         if n <= 0:
#             status = RapidsStatus.s2
#             detail = f'No CUDA GPUs detected (device count = {n}).'
#             report_current()
#             return False
#     except Exception as e:
#         status = RapidsStatus.s2
#         detail = f'GPU check failed: {type(e).__name__}: {e}'
#         report_current()
#         return False

#     status = RapidsStatus.s3
#     detail = f'{n} CUDA GPU(s) detected.'
#     report_current()
#     return True
