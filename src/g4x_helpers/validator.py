import importlib.resources as resources
from pathlib import Path

from pathschema import validate

# register fonts supplied with package
path = resources.files('g4x_helpers.out_schemas')
schemas = {}
for cont in path.iterdir():
    schemas[cont.name.removesuffix('.txt')] = cont


def _print_details(path, result):
    gap = 16
    for err_path, errors in result.errors_by_path.items():
        err_path = Path(err_path)

        relative = err_path.relative_to(path)

        if len(errors) == 0:
            relative = 'root' if str(relative) == '.' else relative
            print(f'{"path valid":<{gap}} - {relative}')

        for err in errors:
            print(f'{err:<{gap}} - {relative}')


def validate_g4x_data(path, schema_name: str, formats: dict = {'sample_id': 'B01'}, report: str = 'short'):
    path = Path(path)

    with open(schemas[schema_name], 'r') as f:
        schema = f.read()

    schema = schema.format(**formats)

    result = validate(path, schema)

    if not result.has_error():
        if report == 'short':
            print('Validation passed with no errors.')
        elif report == 'long':
            _print_details(path, result)

    else:
        if report == 'short':
            print('Validation failed.')
        elif report == 'long':
            _print_details(path, result)
