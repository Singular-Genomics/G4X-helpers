import os

DEFAULT_THREADS = max(1, (os.cpu_count() // 2 or 4))

PX_SIZE_um = 0.32

PROBE_PATTERN = r'^(.*?)-([ACGT]{2,30})-([^-]+)$'

primer_read_map = {
    'SP1': 1,
    'm7a': 2,
    'm9a': 3,
    'm6a': 4,
    'm8a': 5,
    'm3a': 6,
}
