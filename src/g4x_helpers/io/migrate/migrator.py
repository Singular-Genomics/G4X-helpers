import shutil
from pathlib import Path


class LegacyFileMigrator:
    def __init__(self, smp_dir, probes: dict):
        self.smp_dir = Path(smp_dir)
        if 'current' not in probes:
            raise ValueError()
        self.probes = probes
        self.backup_dir = self.smp_dir / 'g4x-helpers' / 'migration_backup'

    @property
    def target_file(self):
        return self.probes['current']

    @property
    def target_file_full(self):
        return self.smp_dir / self.target_file

    @property
    def detected_file(self):
        version = self.get_version()
        if version is None:
            return None
        return self.probes[version]

    @property
    def detected_file_full(self):
        return self.smp_dir / self.detected_file

    @property
    def backup_file(self):
        if self.backup_file_full is not None:
            return self.backup_file_full.relative_to(self.smp_dir)
        return None

    @property
    def backup_file_full(self):
        b_ver = self.get_version(backup=True)
        if b_ver is None:
            return None
        return self.backup_dir / self.probes[b_ver]

    @property
    def backup_exists(self):
        if self.backup_file_full is None:
            return False
        return self.backup_file_full.exists()

    @property
    def matched_versions(self):
        return {version: meta for version, meta in self.probes.items() if (self.smp_dir / meta).exists()}

    @property
    def matched_backup_versions(self):
        return {version: meta for version, meta in self.probes.items() if (self.backup_dir / meta).exists()}

    def get_version(self, backup: bool = False):
        matched = self.matched_backup_versions if backup else self.matched_versions

        if 'current' in matched:
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
        return self.get_version() == 'current'

    @property
    def is_legacy(self):
        return self.get_version() not in (None, 'current')

    # def copy_or_move(self, src, dst, how):
    #     if not dst.parent.exists():
    #         dst.parent.mkdir(parents=True, exist_ok=True)

    #     if how == 'move':
    #         src.rename(dst)
    #     elif how == 'copy':
    #         shutil.copy2(src, dst)
    #     else:
    #         raise ValueError("how must be 'copy' or 'move'")

    def copy_or_move(self, src, dst, how):
        src = Path(src)
        dst = Path(dst)

        if not dst.parent.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)

        if how == 'move':
            shutil.move(src, dst)

        elif how == 'copy':
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dst)

        else:
            raise ValueError("how must be 'copy' or 'move'")

    def backup(self):
        print(f'Creating backup for file: {self.detected_file}')
        if self.backup_exists:
            print(f'Backup file already exists: {self.backup_file}')
            return

        src = self.detected_file_full
        dst = self.backup_dir / self.detected_file
        self.copy_or_move(src, dst, how='copy')

    def revert_backup(self):
        print(f'Reverting backup for file: {self.detected_file}')
        if not self.backup_exists:
            print(f'No backup file found: {self.backup_file}. cannot revert.')
            return

        if self.probes['current'] != '__ABSENT__':
            self.detected_file_full.unlink()

        src = self.backup_file_full
        dst = self.smp_dir / self.backup_file_full.relative_to(self.backup_dir)
        self.copy_or_move(src, dst, how='move')

    def migrate(self, keep_detected: bool = False, backup: bool = True):
        print(f'Starting migration for file: {self.detected_file}')

        version = self.get_version()

        if version is None:
            raise FileNotFoundError('No matching file found.')

        if version == 'current':
            print('File is valid, skipping migration.')
            return True

        if backup:
            self.backup()

        migrate_method = getattr(self, f'_migrate_{version}', None)
        if migrate_method is None:
            raise NotImplementedError(f'No migration implemented for {version}')

        migrate_method(keep_detected=keep_detected)

        print(f'Finished migration for file: {self.detected_file}')
        print(f'File is now valid: {self.is_valid}')
