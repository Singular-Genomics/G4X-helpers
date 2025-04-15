#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import io
import re
import os
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup
import importlib.util

mydir = os.path.dirname(os.path.realpath(__file__))
version_py = os.path.join(mydir, "src/pySing/version.py")
spec = importlib.util.spec_from_file_location("version", version_py)
version = importlib.util.module_from_spec(spec)
spec.loader.exec_module(version)
__version__ = version.__version__


def read(*names, **kwargs):
    with io.open(join(dirname(__file__), *names), encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()

def check_module(module_name, verbose: bool=True):
    """
    Check if a Python module is installed and available for import.

    Parameters:
    ----------
    module_name : str
        The name of the module to check.
    verbose : bool, optional
        If True, prints a message indicating whether the module is available. Defaults to True.

    Returns:
    -------
    bool
        True if the module is available, False otherwise.
    """
    
    if importlib.util.find_spec(module_name) is not None:
        if verbose:
            print(f"{module_name} is installed and available.")
        return True
    else:
        if verbose:
            print(f"{module_name} is not available.")
        return False
    
extra_seq = [
    'biopython',
    'pysam',
    'gget ~= 0.27.9',
    'pyfaidx == 0.8.1.3'
]

extra_scseq = [
    'scanpy == 1.10.2',
    'celltypist ~= 1.6.3',
    'statannot',
    'pybiomart',
    'scrublet',
    'scikit-misc',
    'leidenalg ~= 0.10.2',
    'igraph ~= 0.11.8',
]

extra_gpu = [
    'rapids-singlecell[rapids12] ~= 0.10.11',
    'polars[gpu] ~= 1.8.2'
]

extra_spatial = [
    'spatialdata ~=0.3.0',
    'spatialdata-io ~= 0.1.7'
]

extra_dev = [
    *extra_seq, 
    *extra_scseq,
    *extra_spatial
]


setup(
    name='pySing',
    version=__version__,
    license='MIT',
    description='Python helpers for Singular Genomics.',
    author='Kenneth Gouin III',
    author_email='kennethg@singulargenomics.com',
    url='https://bitbucket.org/singulargenomics/pysing/src/master/',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Private :: Do Not Upload',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.10, <3.11',
    extras_require={
        'seq': extra_seq,
        'scseq': extra_scseq,
        'spatial': extra_spatial,
        'gpu': extra_gpu,
        'dev': extra_dev
    },
    entry_points={
        'console_scripts': [
            'pySing = pySing.cli:main',
            'filter_pair = pySing.utils.filter_pair:launch_filter_pair',
            'aggregate_metrics = pySing.utils.metric_aggregator:launch_aggregator',
            'parse_secondary = pySing.utils.parse_secondary:launch_parser',
            'flowcell_plot = pySing.utils.flowcell_plotter:launch_fc_plotter',
            'brightphase_runner = pySing.utils.brightphase_runner:launch_brightphase'
        ]
    },
)
