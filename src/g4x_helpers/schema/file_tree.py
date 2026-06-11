import copy
from pathlib import Path
from typing import Literal

from . import definition as sd
from . import utils as ut

MAIN_VALIDATOR = sd.SampleG4X

ASSAY_AGNOSTIC_VALIDATORS = (
    sd.SampleSheet,
    sd.QCSummary,
    sd.HnEDir,
    sd.Segmentation,
    sd.BeadMask,
    sd.SingleCellFolder,
    sd.CellMetadata,
    sd.AdataH5,
    sd.ViewerZarr,
)

TX_VALIDATORS = (
    sd.Manifest,
    sd.RawFeatures,
    sd.TxTable,
    sd.CellxGene,
    sd.ClusteringUmap,
    sd.Dgex,
)

PR_VALIDATORS = (
    sd.ProteinPanel,
    sd.ProteinDir,
    sd.CellxProt,
)


class FlatTree:
    def __init__(self, sample_dir: str):

        self.smp_dir = Path(sample_dir)
        validators = self.fetch_validators([MAIN_VALIDATOR])
        validators += self.fetch_validators(ASSAY_AGNOSTIC_VALIDATORS)
        validators += self.fetch_validators(TX_VALIDATORS)
        validators += self.fetch_validators(PR_VALIDATORS)

        self.validators = validators

        for v in self.validators:
            setattr(self, v.name, v)

    def fetch_validators(self, validators: tuple):
        fetched = []

        for validator_cls in validators:
            fetched.append(validator_cls(root=self.smp_dir))

        return fetched


class FileTree:
    def __init__(self, sample_dir: str, alt_source: str | None = None):

        self.smp_dir = Path(sample_dir)
        self.alt_source = Path(alt_source) if alt_source else None

        meta_validator = MAIN_VALIDATOR(root=self.smp_dir)
        if not meta_validator.path_exists():
            msg = 'Missing sample.g4x\n'
            msg += 'G4X-helpers 4 requires that G4X-data must contain a metadata file named "sample.g4x"\n\n'
            msg += 'If this data was generated with a software version prior to 26.1, you can convert it to the lastest schema using "g4x-helpers migrate"'
            raise ValidationError(msg)

        if not meta_validator.is_valid:
            raise ValidationError(f'sample.g4x is not valid\nCaused by: {meta_validator.validation()}')

        assay_type = ut.detect_assay_type(meta_validator.load())
        if assay_type == 'undefined':
            raise ValidationError('Could not detect assay type from sample.g4x')

        self.assay_type = assay_type
        self.tx_detected = True if assay_type in ['combined', 'tx_only'] else False
        self.pr_detected = True if assay_type in ['combined', 'pr_only'] else False

        validators = [meta_validator]
        validators += self.fetch_validators(ASSAY_AGNOSTIC_VALIDATORS)
        if self.tx_detected:
            validators += self.fetch_validators(TX_VALIDATORS)
        if self.pr_detected:
            validators += self.fetch_validators(PR_VALIDATORS)

        self.validators = validators
        self.raw_validators = [v for v in self.validators if v.PRIMARY]
        self.secondary_validators = [v for v in self.validators if not v.PRIMARY]

        for v in self.validators:
            setattr(self, v.name, v)

    def fetch_validators(self, validators: tuple):
        fetched = []

        for validator_cls in validators:
            if self.alt_source:
                validator = validator_cls(root=self.alt_source)

                if validator.is_valid:
                    fetched.append(validator)
                    continue

            fetched.append(validator_cls(root=self.smp_dir))

        return fetched

    def copy(self):
        return copy.deepcopy(self)

    @property
    def is_valid_raw(self):
        return all([v.is_valid for v in self.raw_validators])

    @property
    def is_valid_all(self):
        return all([v.is_valid for v in self.validators])

    @property
    def errors(self):
        return {v.name: v.validation() for v in self.validators if not v.is_valid}

    @property
    def reports(self):
        return {v.name: v.validation() for v in self.validators}

    def _val_report_minimal(self, raw_only: bool = True, report_pass: bool = True, raise_exception: bool = True):
        gate = self.is_valid_raw if raw_only else self.is_valid_all
        what = 'raw' if raw_only else 'all'
        validation_title = f'G4X-[{what} data] validation'
        if not gate:
            msg = f'{validation_title} failed for:\n'
            msg += f'{self.smp_dir}\n'

            errs = pretty_dict_str(self.errors)
            msg += f'\nErrors:\n{errs}'
            if raise_exception:
                raise ValidationError(msg)
            else:
                return msg
        else:
            if report_pass:
                return f'{validation_title} passed for:\n{self.smp_dir}'

    def _val_report_verbose(self, raw_only: bool = False, raise_exception: bool = True):

        if self.SampleG4X.path_exists():
            msg = f'Detected G4X-metadata file:\n{self.SampleG4X.target_path.resolve()}'
            msg += f'\nassay type: {self.assay_type}'
            msg += '\n\n> Validating required raw data ...'
        else:
            return self.SampleG4X.report_validation()

        for v in self.raw_validators:
            msg += f'\n{v.report_validation()}'

        if not raw_only:
            msg += '\n\n> Validating secondary data ...'
            for v in self.secondary_validators:
                msg += f'\n{v.report_validation()}'

        msg += '\n\n'
        msg += self._val_report_minimal(raw_only=raw_only, raise_exception=raise_exception)
        return msg

    def validation_report(
        self,
        format: Literal['verbose', 'minimal'] = 'verbose',
        raw_only: bool = False,
        report_pass: bool = True,
        raise_exception: bool = True,
    ):
        if format == 'verbose':
            msg = self._val_report_verbose(raw_only=raw_only, raise_exception=raise_exception)
        elif format == 'minimal':
            msg = self._val_report_minimal(raw_only=raw_only, report_pass=report_pass, raise_exception=raise_exception)
        else:
            raise ValueError(f"Invalid format: {format}. Expected 'verbose' or 'minimal'.")

        if msg is None:
            return None

        return msg


class ValidationError(Exception):
    pass


def kv_line_gap(key, value, gap=2):
    value = '<undefined>' if not value else value
    line = f'{key:<{gap}}'
    line += ' - '
    line += f'{value}'

    return line


def pretty_dict_str(d):
    max_len = max([len(k) for k in d.keys()])
    msg = ''
    for k, v in d.items():
        msg += kv_line_gap(k, v, gap=max_len) + '\n'
    return msg
