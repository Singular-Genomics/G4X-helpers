from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from importlib import metadata
from typing import Optional


class RapidsStatus(str, Enum):
    s0 = 'not_installed'  # metadata missing
    s1 = 'import_failed'  # dist present but module import fails
    s2 = 'no_gpu'  # import ok but no usable CUDA GPU
    s3 = 'ok'  # ready to use


@dataclass(frozen=True)
class RapidsCheck:
    status: RapidsStatus
    rsc_version: Optional[str] = None
    detail: Optional[str] = None

    def __repr__(self):
        gap = 11
        repr_str = f'{"status":<{gap}} : {self.status[0:]}\n'
        repr_str += f'{"rsc_version":<{gap}} : {self.rsc_version}\n'
        repr_str += f'{"detail":<{gap}} : {self.detail}'

        return repr_str


def check_rapids() -> RapidsCheck:
    # --- Stage 1: metadata presence (fast) ---
    try:
        ver = metadata.version('rapids-singlecell')
    except metadata.PackageNotFoundError:
        return RapidsCheck(
            status=RapidsStatus.s0,
            detail="Distribution 'rapids-singlecell' not found. Install with `uv sync --extra rapids`. (linux only)",
        )

    # --- Stage 2: import module (realistic) ---
    try:
        import rapids_singlecell as _  # noqa: F401
    except Exception as e:
        return RapidsCheck(
            status=RapidsStatus.s1,
            rsc_version=ver,
            detail=f'Import failed: {type(e).__name__}: {e}',
        )

    # --- Stage 3: verify a CUDA GPU is usable ---
    # You can choose your preferred “GPU probe”. Cupy is a common, lightweight-ish probe.
    try:
        import cupy as cp  # many RAPIDS stacks already include it; if not, this will fail clearly

        n = int(cp.cuda.runtime.getDeviceCount())
        if n <= 0:
            return RapidsCheck(
                status=RapidsStatus.s2,
                rsc_version=ver,
                detail=f'No CUDA GPUs detected (device count = {n}).',
            )
    except Exception as e:
        return RapidsCheck(
            status=RapidsStatus.s2,
            rsc_version=ver,
            detail=f'GPU check failed: {type(e).__name__}: {e}',
        )

    return RapidsCheck(status=RapidsStatus.s3, rsc_version=ver, detail=f'{n} CUDA GPU(s) detected.')


### >>> usage
# res = check_rapids()

# if res.stage == RapidsStage.OK:
#     # use RAPIDS backend
