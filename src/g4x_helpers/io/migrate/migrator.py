import re
import shutil
from pathlib import Path


class MigrationError(Exception):
    pass


class DataMigrator:
    ABSENT_SENTINEL = '__ABSENT__'

    def __init__(self, smp_dir, probes: dict, name: str | None = None):
        self.name = name or '_unnamed-migrator_'
        self.smp_dir = Path(smp_dir)

        self.probes = self._validate_probes(probes)
        self.backup_dir = self.smp_dir / 'g4x-helpers' / 'migration_backup'

    @property
    def target_path(self):
        return self.probes['current']

    @property
    def target_path_full(self):
        return self.smp_dir / self.target_path

    @property
    def detected_path(self):
        version = self.get_version()
        if version is None:
            return None
        return self.probes[version]

    @property
    def detected_path_full(self):
        return self.smp_dir / self.detected_path

    @property
    def backup_path(self):
        if self.backup_path_full is not None:
            return self.backup_path_full.relative_to(self.smp_dir)
        return None

    @property
    def backup_path_full(self):
        b_ver = self.get_version(backup=True)
        if b_ver is None:
            return None
        return self.backup_dir / self.probes[b_ver]

    @property
    def backup_exists(self):
        if self.backup_path_full is None:
            return False
        return self.backup_path_full.exists()

    @property
    def has_absent_target(self):
        return self.target_path == self.ABSENT_SENTINEL

    @property
    def matched_versions(self):
        return {version: meta for version, meta in self.probes.items() if (self.smp_dir / meta).exists()}

    @property
    def matched_backup_versions(self):
        return {version: meta for version, meta in self.probes.items() if (self.backup_dir / meta).exists()}

    def get_version(self, backup: bool = False):
        matched = self.matched_backup_versions if backup else self.matched_versions

        if 'current' in matched:
            # TODO: if both current and legacy versions coexist, optionally run cleanup logic for stale legacy paths.
            return 'current'

        legacy_versions = [v for v in matched if v != 'current']
        if len(legacy_versions) == 1:
            return legacy_versions[0]

        if len(legacy_versions) > 1:
            # print(f'Multiple legacy files found: {legacy_versions}. Using latest.')
            return max(legacy_versions, key=lambda v: int(v.rsplit('_', 1)[-1]))

        return None

    @property
    def is_valid(self):
        if self.has_absent_target:
            return self.get_version() is None
        return self.get_version() == 'current'

    def relocate(self, src, dst, how):
        src = Path(src)
        dst = Path(dst)

        if not dst.parent.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)

        if how == 'move':
            shutil.move(src, dst)

        elif how == 'copy':
            if src.is_dir():
                # TODO: dirs_exist_ok=True merges into existing backup dirs; consider strict backups to avoid stale files.
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dst)

        else:
            raise ValueError("how must be 'copy' or 'move'")

    def backup(self):
        print(f'{self.name}: creating backup for file: {self.detected_path}')
        if self.detected_path is None:
            raise MigrationError(f'{self.name}: No matching file found to backup.')

        if self.backup_exists:
            print(f'{self.name}: Backup file already exists: {self.backup_path}')
            return

        src = self.detected_path_full
        dst = self.backup_dir / self.detected_path
        self.relocate(src, dst, how='copy')

    def revert_backup(self):
        print(f'{self.name}: Reverting backup for file: {self.detected_path}')
        if not self.backup_exists:
            print(f'{self.name}: No backup file found: {self.backup_path}. cannot revert.')
            return

        if not self.has_absent_target and self.target_path_full.exists():
            if self.target_path_full.is_dir():
                shutil.rmtree(self.target_path_full)
            else:
                self.target_path_full.unlink()

        src = self.backup_path_full
        dst = self.smp_dir / self.backup_path_full.relative_to(self.backup_dir)
        self.relocate(src, dst, how='move')

    def migrate(self, backup: bool = True):
        print(f'{self.name}: Starting migration for path: {self.detected_path}')

        version = self.get_version()

        if version is None:
            raise MigrationError(f'{self.name}: No matching file found.')

        if version == 'current':
            print(f'{self.name}: path is valid, skipping migration.')
            return True

        # check for migration method first to the backup does not run if this fails
        migrate_method = self._get_migrate_method(version)

        if backup:
            self.backup()

        migrate_method()

        print(f'{self.name}: Finished migration for: {self.detected_path}')
        print(f'{self.name}: Path is now valid: {self.is_valid}')
        return True

    def _validate_probes(self, probes):
        if 'current' not in probes:
            raise ValueError(f'{self.name}: The "current" key must be present in the probes dictionary.')

        pattern = r'^legacy_version_\d+$'

        def is_valid(s):
            return bool(re.match(pattern, s))

        invalid_keys = [k for k in probes.keys() if k != 'current' and not is_valid(k)]

        if invalid_keys:
            raise ValueError(
                f'{self.name}: The following legacy keys are invalid: {invalid_keys}, please use the format "legacy_version_<number>".'
            )
        return probes

    def _get_migrate_method(self, version):
        migrate_method = getattr(self, f'_migrate_{version}', None)
        if migrate_method is None:
            raise MigrationError(
                f'{self.name}: No migration implemented for {version}.\nPlease create a subclass with a migration method named "_migrate_{version}"'
            )
        return migrate_method
