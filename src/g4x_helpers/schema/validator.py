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


class BaseValidator:
    DEFAULT_TARGET_PATH = '.'
    VALIDATION_TESTS = ()

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.register_tests()

    def __init__(
        self,
        target_path: str | None = None,
        root: str | None = None,
        resolve: bool = False,
        validate_absence: bool = False,
    ):
        self.name = type(self).__name__  # + 'Validator'
        self.root = Path(root) if root is not None else None
        self._target_path = Path(target_path) if target_path is not None else Path(type(self).DEFAULT_TARGET_PATH)

        if self.root is not None:
            self._target_path = self.root / self._target_path

        if resolve:
            self._target_path = self._target_path.resolve(strict=False)

        self.validate_absence = validate_absence

    @classmethod
    def register_tests(cls):
        tests = []
        for base in reversed(cls.__mro__[1:]):
            tests.extend(getattr(base, 'VALIDATION_TESTS', ()))

        for name, value in cls.__dict__.items():
            if getattr(value, '_is_validation_test', False):
                tests.append(name)

        cls.VALIDATION_TESTS = tuple(dict.fromkeys(tests))

    @property
    def target_name(self):
        return self.target_path.name

    @property
    def target_path(self):
        return self._target_path.expanduser()

    @property
    def p(self):
        return self.target_path

    @property
    def n(self):
        return self.target_name

    @target_path.setter
    def target_path(self, value):
        self._target_path = Path(value)

    @property
    def target_rel(self):
        return self.target_path.relative_to(self.root) if self.root is not None else self.target_path.name

    @property
    def target_type(self):
        if self.target_path.is_file():
            return 'file'
        elif self.target_path.is_dir():
            return 'directory'
        else:
            return 'unknown'

    @property
    def is_valid(self):
        return all(self.validation().values())

    @property
    def path_exists(self):
        return self.target_path.exists()

    def validation(self):
        validation_results = {}

        if self.validate_absence:
            validation_results['absent'] = not self.path_exists
            return validation_results

        validation_results['path_exists'] = self.path_exists

        if self.VALIDATION_TESTS and self.path_exists:
            for test_name in self.VALIDATION_TESTS:
                validation_results[test_name] = getattr(self, test_name)()

        return validation_results

    def report_validation(self):
        if not self.is_valid:
            code = '[!invalid]'
            msg = f'{code} {self.name} {self.target_type}: {self.target_rel}\n'
            msg += f'{len(code) * " "} reason: {self.validation()}'
            return msg
        else:
            return f'[valid] {self.name} {self.target_type}: {self.target_rel}'
