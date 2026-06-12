import logging
import re
import shutil

from ... import io
from ..validator import BaseValidator, DirectoryValidator, FileValidator

LOGGER = logging.getLogger(__name__)


class MigrationError(Exception):
    pass


class DataMigrator(BaseValidator):
    IS_OPTIONAL = False
    VERSION_VALIDATORS = {}
    VERSION_CLASS_PATTERN = re.compile(r'.*_V\d+$')
    COPY_CURRENT = True

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        auto_validators = {}
        for name, obj in cls.__dict__.items():
            if isinstance(obj, type) and issubclass(obj, BaseValidator) and cls.VERSION_CLASS_PATTERN.match(name):
                key = name[0] + name[1:]

                auto_validators[key] = obj

        explicit_validators = getattr(cls, 'VERSION_VALIDATORS', {})

        cls.VERSION_VALIDATORS = {
            **auto_validators,
            **explicit_validators,
        }

    def _current_validator(self):
        for cls in type(self).__mro__:
            if cls is DataMigrator:
                continue

            if issubclass(cls, BaseValidator):
                return cls(root=self.root)

        raise TypeError(f'Could not infer "current" validator for {self._name}')

    def versions(self):
        versions = {'current': self._current_validator()}
        versions.update(
            {name: validator_cls(root=self.root) for name, validator_cls in type(self).VERSION_VALIDATORS.items()}
        )
        return versions

    @property
    def _name(self):
        return type(self).__name__.removesuffix('_Migrator')

    @property
    def target_validator(self):
        return type(self).__bases__[1]

    @property
    def valid_versions(self):
        return list({k for k, v in self.versions().items() if v.is_valid})

    @property
    def mig_version(self):
        if len(self.valid_versions) == 0:
            if self.IS_OPTIONAL:
                return 'missing_optional'
            else:
                raise ValueError(f'No valid versions found for {self._name} in {self.root}')
        if 'current' in self.valid_versions:
            mig_version = 'current'
        else:
            sorted_legacy_versions = sorted(self.valid_versions, key=lambda x: int(x.split('_V')[-1]))
            mig_version = sorted_legacy_versions[-1]
        return mig_version

    @property
    def migrator(self):
        return self.versions()[self.mig_version]

    @property
    def is_migratable(self):
        return self.migration_status[0]

    @property
    def migration_status(self):
        if self.IS_OPTIONAL:
            return True, f'{self._name} is optional, migration not required.'

        has_legacy = len(self.valid_versions) > 0
        if not has_legacy:
            msg = f'{self._name} has no valid versions detected.'
        else:
            msg = f'{self._name} has migratable versions: {self.valid_versions}.'

        return has_legacy, msg

    @property
    def is_file(self):
        return isinstance(self, FileValidator)

    @property
    def is_folder(self):
        return isinstance(self, DirectoryValidator)

    def copy_if_current(self, out_path):
        if self.is_file:
            file_out = io.pathval.ensure_parent_dir(out_path / self.DEFAULT_TARGET_PATH)
            shutil.copy(self.migrator.p, file_out)
        elif self.is_folder:
            shutil.copytree(self.migrator.p, out_path / self.DEFAULT_TARGET_PATH)
        else:
            raise TypeError(f'Cannot copy {self._name} because it is neither a file nor a folder.')

    def migrate(self, out_path, logger: logging.Logger | None = None, *args, **kwargs):
        log = logger or LOGGER
        kwargs['logger'] = log
        log.info(f'Initializing migration of {self._name}')

        if not self.valid_versions and self.IS_OPTIONAL:
            log.warning(f'No valid versions of optional {self._name} found. Skipping migration.')
            return

        if not self.is_migratable:
            raise ValueError(f'{self._name} is not migratable. {self.migration_status[1]}')

        out_path = io.pathval.validate_dir_path(out_path)

        if self.COPY_CURRENT and 'current' in self.valid_versions:
            log.debug(f'{self._name} has correct schema, copying data without migration.')
            self.copy_if_current(out_path)
            return

        log.debug(f'Detected legacy version: "{self.mig_version}"')

        migrate_method = getattr(self, '_migrate_method', False)
        if not migrate_method:
            raise NotImplementedError(f'{self._name} does not have a _migrate_method defined.')
        else:
            try:
                result = migrate_method(out_path, *args, **kwargs)
                target = self.target_validator(root=out_path)
                if target.is_valid:
                    log.debug(f'✓ Migration successful for {self._name}!')
                else:
                    raise ImportError(f'Migrated {self._name} did not pass validation! {target.validation()}.')
                return result

            except Exception as e:
                raise ImportError(f'Error migrating {self._name} from {self.target_path}!\nreason: {e}') from None
