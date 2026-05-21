from __future__ import annotations

from pathlib import Path

from .. import c, io


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
        # format: dict | None = None,
        validate_absence: bool = False,
    ):
        self.name = type(self).__name__
        self.root = Path(root) if root is not None else None
        self.resolve = resolve

        # self.format = {'sample_id': 'g4x_sample'} if format is None else format
        # tpath_str = str(target_path) if target_path is not None else str(type(self).DEFAULT_TARGET_PATH)
        # tpath_str = tpath_str.format(**self.format)
        target_path = Path(target_path) if target_path is not None else Path(type(self).DEFAULT_TARGET_PATH)

        self._target_path = target_path
        self.validate_absence = validate_absence

    @classmethod
    def register_tests(cls):
        tests = []

        # Scan the whole MRO so decorated tests on BaseValidator are included too
        for base in reversed(cls.__mro__):
            for name, value in base.__dict__.items():
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

        if '*' in complete_path.name:
            candidates = list(complete_path.parent.glob(complete_path.name))
            if len(candidates) > 1:
                raise ValueError(f'Multiple files found matching pattern {complete_path.name}: {candidates}')
            elif len(candidates) == 0:
                pass
                # raise ValueError(f'No files found matching pattern {complete_path.name}')
            else:
                complete_path = candidates[0]

        if self.resolve:
            complete_path = complete_path.resolve(strict=False)

        return complete_path.expanduser()

    @target_path.setter
    def target_path(self, value):
        self._target_path = Path(value)

    @property
    def p(self):
        return self.target_path

    @property
    def n(self):
        return self.target_name

    @property
    def target_rel(self):
        return self.target_path.relative_to(self.root) if self.root is not None else self.target_path.name

    @property
    def target_type(self):
        if self.target_path.is_file():
            return 'file'
        elif self.target_path.is_dir():
            return 'directory'
        return ''

    @validation_test
    def path_exists(self):
        return self.target_path.exists()

    @property
    def is_valid(self):
        return all(self.validation().values())

    def validation(self):
        validation_results = {}

        if self.validate_absence:
            validation_results['absent'] = not self.path_exists()
            return validation_results

        validation_results['path_exists'] = self.path_exists()

        if self.VALIDATION_TESTS and self.path_exists():
            for test_name in self.VALIDATION_TESTS:
                validation_results[test_name] = getattr(self, test_name)()

        return validation_results

    def report_validation(self):
        if not self.is_valid:
            code = '[!invalid]'
            msg = f'{code} {self.name} {self.target_type}: {self.target_rel}\n'
            msg += f'{" " * len(code)} reason: {self.validation()}'
            return msg
        return f'[valid] {self.name} {self.target_type}: {self.target_rel}'


class FileValidator(BaseValidator):
    def load(self, *args, **kwargs):
        return self._try_load(*args, **kwargs)

    def _try_load(self, *args, **kwargs):
        load_method = getattr(self, '_load_method', False)
        if not load_method:
            raise NotImplementedError(f'{type(self).__name__} does not have a _load_method defined.')
        else:
            try:
                return load_method(*args, **kwargs)
            except Exception as e:
                raise ImportError(f'Error loading {type(self).__name__} from {self.target_rel}!\nreason: {e}') from None


class TableValidator(FileValidator):
    SCHEMA = {}

    def _load_method(self, lazy: bool = False, use_cache: bool = False):
        return io.import_table(self.target_path, lazy=lazy, use_cache=use_cache)

    @validation_test
    def is_table(self):
        suffix = self.target_path.name.removesuffix('.gz').split('.')[-1]
        return suffix in {'csv', 'tsv', 'parquet'}

    @validation_test
    def correct_schema(self):
        lf = self.load(lazy=True)
        lf_names = lf.collect_schema().names()
        return set(self.SCHEMA.keys()).issubset(lf_names)


class FolderValidator(BaseValidator):
    EXPECTED_FILES = {}
    ALT_FILES = {}
    EXPECTED_DIRS = {}

    def __init__(self, root=None, target_path=None, **kwargs):
        super().__init__(root=root, target_path=target_path, **kwargs)

    @property
    def expected_present(self):
        expected_present = True
        for f in self.EXPECTED_FILES:
            if f not in self.existing_file_names():
                expected_present = False
        return expected_present

    @property
    def alt_present(self):
        if self.ALT_FILES:
            alt_present = True
            for f in self.ALT_FILES:
                if f not in self.existing_file_names():
                    alt_present = False
        else:
            alt_present = False
        return alt_present

    @validation_test
    def is_directory(self):
        return self.target_type == 'directory'

    @validation_test
    def dirs_present(self):
        for d in self.EXPECTED_DIRS:
            if d not in self.existing_directory_names():
                return False
        return True

    @validation_test
    def files_present(self):
        return self.expected_present or self.alt_present

    def existing_files(self):
        return [p for p in self.p.iterdir() if p.is_file()]

    def existing_directories(self):
        return [p for p in self.p.iterdir() if p.is_dir()]

    def existing_file_names(self):
        return [p.name for p in self.existing_files()]

    def existing_directory_names(self):
        return [p.name for p in self.existing_directories()]


class DirectoryValidator(BaseValidator):
    EXPECTED_FILES = {}
    EXPECTED_DIRS = {}
    FILE_MAP = {}

    def __init__(self, root=None, target_path=None, **kwargs):
        super().__init__(root=root, target_path=target_path, **kwargs)

    def existing_files(self):
        return [p for p in self.p.iterdir() if p.is_file()]

    def existing_directories(self):
        return [p for p in self.p.iterdir() if p.is_dir()]

    def existing_file_names(self):
        return [p.name for p in self.existing_files()]

    def existing_directory_names(self):
        return [p.name for p in self.existing_directories()]

    @property
    def requested_files(self):
        mapped_expected = {f: [f] for f in self.EXPECTED_FILES}
        return self.FILE_MAP | mapped_expected

    @property
    def mapped_files(self):
        existing_files = {}
        for file_name, file_paths in self.requested_files.items():
            for f in file_paths:
                query = self.p / f
                if query.exists():
                    existing_files[file_name] = query
                    break
        return existing_files

    @validation_test
    def is_directory(self):
        return self.target_type == 'directory'

    @validation_test
    def dirs_present(self):
        for d in self.EXPECTED_DIRS:
            if d not in self.existing_directory_names():
                return False
        return True

    @validation_test
    def files_present(self):
        return set(self.mapped_files.keys()).issuperset(set(self.requested_files.keys()))


class ImgDirectoryValidator(DirectoryValidator):
    VALID_IMG_TYPES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]

    @property
    def img_type(self):
        ome_files = sum([f.endswith(c.PREFERRED_IMG_SUFFIX) for f in self.existing_file_names()])
        jp2_files = sum([f.endswith(c.ALT_IMG_SUFFIX) for f in self.existing_file_names()])

        img_type = 'unknown'
        if ome_files > 0 and jp2_files > 0:
            img_type = 'mixed'
        elif ome_files > 0 and jp2_files == 0:
            img_type = '.ome.tiff'
        elif jp2_files > 0 and ome_files == 0:
            img_type = '.jp2'
        return img_type

    @validation_test
    def img_type_valid(self):
        return self.img_type in self.VALID_IMG_TYPES

    def get_img(self, query: str):
        return self.mapped_files[query]
