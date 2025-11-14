import importlib.resources as resources
from pathlib import Path

from pathschema import validate

# register fonts supplied with package
path = resources.files('g4x_helpers.schemas')
schemas = {}
for cont in path.iterdir():
    schemas[cont.name.removesuffix('.txt')] = cont


def _print_details(path, result, errors_only=True):
    gap = 21
    for err_path, errors in result.errors_by_path.items():
        err_path = Path(err_path)

        relative = err_path.relative_to(path)

        if not errors_only:
            if len(errors) == 0:
                relative = 'root' if str(relative) == '.' else relative
                print(f'{"path valid":<{gap}} - {relative}')

        for err in errors:
            if err.startswith('Missing required file'):
                err_full = err
                err = 'Missing required file'
                relative = relative / err_full.removeprefix('Missing required file ')

            print(f'{err:<{gap}} - {relative}')


def validate_g4x_data(
    path,
    schema_name: str,
    formats: dict | None = None,
    report: str | None = 'short',
):
    path = Path(path)

    if formats is None:
        formats = {'sample_id': 'A01'}

    with open(schemas[schema_name], 'r') as f:
        schema = f.read()

    schema = schema.format(**formats)
    result = validate(path, schema)

    ok = not result.has_error()

    # store callables, not results
    reports = {
        True: {
            'short': lambda: _print_details(path, result, errors_only=True),
            'long': lambda: _print_details(path, result, errors_only=False),
        },
        False: {
            'short': lambda: _print_details(path, result, errors_only=True),
            'long': lambda: _print_details(path, result, errors_only=False),
        },
    }

    if report:
        try:
            reports[ok][report]()  # <-- call the function
        except KeyError:
            raise ValueError(f'Unknown report type: {report!r}')

    return ok
