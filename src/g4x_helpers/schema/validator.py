from pathlib import Path


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
        format: dict = {'sample_id': 'g4x_sample'},
        validate_absence: bool = False,
    ):
        self.name = type(self).__name__  # + 'Validator'
        self.root = Path(root) if root is not None else None
        self.resolve = resolve

        tpath_str = str(target_path) if target_path is not None else str(type(self).DEFAULT_TARGET_PATH)
        tpath_str = tpath_str.format(**format)

        self._target_path = Path(tpath_str)

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
        if self.root is not None:
            complete_path = self.root / self._target_path
        else:
            complete_path = self._target_path

        if self.resolve:
            complete_path = complete_path.resolve(strict=False)

        return complete_path.expanduser()

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
            return ''

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
