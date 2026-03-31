import copy
import json
from pathlib import Path
from typing import Literal

from .definition import (
    AdataH5,
    BeadMask,
    CellMetadata,
    CellxGene,
    CellxProt,
    ClusteringUmap,
    Dgex,
    HnEDir,
    Manifest,
    ProteinDir,
    ProteinPanel,
    QCSummary,
    RawFeatures,
    SampleMetadata,
    SampleSheet,
    Segmentation,
    SingleCellFolder,
    TxTable,
    ViewerZarr,
)


class FileTree:
    def __init__(self, sample_dir):
        self.smp_dir = Path(sample_dir)

        meta_validator = SampleMetadata(root=self.smp_dir)

        self.raw_validators = [meta_validator]
        self.secondary_validators = []

        self.tx_detected = False
        self.pr_detected = False

        if meta_validator.is_valid:
            sample_id = meta_validator.load()['sample_id']
            self.detect_assay_type(meta_validator.target_path)

            self.raw_validators.extend(
                [
                    SampleSheet(root=self.smp_dir),
                    Segmentation(root=self.smp_dir),
                    BeadMask(root=self.smp_dir),
                    HnEDir(root=self.smp_dir),
                ]
            )

            self.secondary_validators.extend(
                [
                    QCSummary(root=self.smp_dir, format={'sample_id': sample_id}),
                    ViewerZarr(root=self.smp_dir),
                    SingleCellFolder(root=self.smp_dir),
                    CellMetadata(root=self.smp_dir),
                    AdataH5(root=self.smp_dir),
                    ClusteringUmap(root=self.smp_dir),
                    Dgex(root=self.smp_dir),
                ]
            )

            if self.tx_detected:
                self.raw_validators.extend(
                    [
                        Manifest(root=self.smp_dir),
                        RawFeatures(root=self.smp_dir),
                    ]
                )
                self.secondary_validators.extend(
                    [
                        TxTable(root=self.smp_dir),
                        CellxGene(root=self.smp_dir),
                    ]
                )

            if self.pr_detected:
                self.raw_validators.extend(
                    [
                        ProteinPanel(root=self.smp_dir),
                        ProteinDir(root=self.smp_dir),
                    ]
                )
                self.secondary_validators.extend(
                    [
                        CellxProt(root=self.smp_dir),
                    ]
                )

        self.validators = self.raw_validators + self.secondary_validators

        for v in self.validators:
            setattr(self, v.name, v)

    def detect_assay_type(self, meta_path):
        with open(meta_path, 'r') as f:
            smp_meta = json.load(f)

        # TODO ensure these hooks are a good choice
        # inspect tx and pr presence from valid metadata
        self.tx_detected = 'transcript_panel' in smp_meta
        self.pr_detected = 'protein_panel' in smp_meta

        self.assay_type = (
            'combined'
            if self.tx_detected and self.pr_detected
            else 'tx_only'
            if self.tx_detected
            else 'pr_only'
            if self.pr_detected
            else 'undefined'
        )

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

        if self.SampleMetadata.path_exists:
            msg = f'Detected G4X-metadata file:\n{self.SampleMetadata.target_path}'
            msg += f'\nassay type: {self.assay_type}'
            msg += '\n\n> Validating required raw data ...'
        else:
            return self.SampleMetadata.report_validation()

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
