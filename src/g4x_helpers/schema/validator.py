from pathlib import Path


class DataValidator:
    DEFAULT_TARGET_PATH = '.'
    TYPE = 'file'
    VALIDATION_TESTS = ()

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        tests = []
        for base in reversed(cls.__mro__[1:]):
            tests.extend(getattr(base, 'VALIDATION_TESTS', ()))

        for name, value in cls.__dict__.items():
            if getattr(value, '_is_validation_test', False):
                tests.append(name)

        cls.VALIDATION_TESTS = tuple(dict.fromkeys(tests))

    def __init__(self, smp_dir, target_path: str | None = None, validate_absence: bool = False):
        self.smp_dir = Path(smp_dir)
        self._target_path = Path(target_path) if target_path is not None else Path(type(self).DEFAULT_TARGET_PATH)
        self.name = type(self).__name__  # + 'Validator'
        self.validate_absence = validate_absence
        self.backup_dir = self.smp_dir / 'g4x_helpers' / 'migration_backup'

    @property
    def target_path(self):
        return self._target_path

    @target_path.setter
    def target_path(self, value):
        self._target_path = Path(value)

    @property
    def target_path_resolved(self):
        rel_path = self.target_path if not self.has_wildcard else self.resolve_wildcard()
        return (self.smp_dir / rel_path).resolve()

    @property
    def path(self):
        return self.target_path_resolved

    @property
    def has_wildcard(self):
        parts = self.target_path.parts
        return any('*' in part or '?' in part or '[' in part for part in parts)

    @property
    def path_exists(self):
        return self.target_path_resolved.exists()

    def resolve_wildcard(self):
        if self.has_wildcard:
            matches = sorted(self.smp_dir.glob(str(self.target_path)))

            if len(matches) != 1:
                # print(
                #     f'{self.name}: Could not resolve wildcard, expected exactly 1 match for {self.target_path!s}, found {len(matches)}: {matches}'
                # )
                return self.target_path
            return matches[0].relative_to(self.smp_dir)
        else:
            return self.target_path

    @property
    def is_valid(self):
        return all(self.validation().values())

    def validation(self):
        validation_results = {
            # 'resolved_path': self.target_path,  #
            'path_exists': self.path_exists,
        }

        if self.VALIDATION_TESTS and self.path_exists:
            for test_name in self.VALIDATION_TESTS:
                validation_results[test_name] = getattr(self, test_name)()

        if self.validate_absence:
            validation_results['absent'] = not self.path_exists

        return validation_results

    def report_validation(self):
        val_name = self.name.removesuffix('Validator')
        resolved_target_name = self.target_path_resolved.relative_to(self.smp_dir)
        if not self.is_valid:
            code = '[!invalid]'
            msg = f'{code} {val_name} {self.TYPE}: {resolved_target_name}\n'
            msg += f'{len(code) * " "} reason: {self.validation()}'
            print(msg)
        else:
            print(f'[valid] {val_name} {self.TYPE}: {resolved_target_name}')


def validation_test(func):
    func._is_validation_test = True
    return func
