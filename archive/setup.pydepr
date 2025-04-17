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
version_py = os.path.join(mydir, "src/g4x_helpers/version.py")
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

setup(
    name='g4x_helpers',
    version=__version__,
    license='MIT',
    description='Python helpers for G4X.',
    author='Kenneth Gouin III',
    author_email='kennethg@singulargenomics.com',
    url='https://github.com/Singular-Genomics/G4X-helpers',
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
)
