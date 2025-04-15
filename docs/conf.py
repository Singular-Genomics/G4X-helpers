# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../src'))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
]

## get version information
import importlib.util
version_py = os.path.join(os.path.abspath('../'), "src/g4x_helpers/version.py")
spec = importlib.util.spec_from_file_location("version", version_py)
version = importlib.util.module_from_spec(spec)
spec.loader.exec_module(version)
__version__ = version.__version__

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

source_suffix = '.rst'
master_doc = 'index'
project = 'g4x-helpers'
repository_url = 'https://github.com/Singular-Genomics/G4X-helpers'
copyright = '2025, Singular Genomics'
author = 'Kenneth Gouin III'
version = release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['.']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'
html_title = "G4X-helpers"
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
    '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_ivar = True
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    bbknn=("https://bbknn.readthedocs.io/en/latest/", None),
    cuml=("https://docs.rapids.ai/api/cuml/stable/", None),
    cycler=("https://matplotlib.org/cycler/", None),
    dask=("https://docs.dask.org/en/stable/", None),
    dask_ml=("https://ml.dask.org/", None),
    h5py=("https://docs.h5py.org/en/stable/", None),
    ipython=("https://ipython.readthedocs.io/en/stable/", None),
    igraph=("https://python.igraph.org/en/stable/api/", None),
    leidenalg=("https://leidenalg.readthedocs.io/en/latest/", None),
    louvain=("https://louvain-igraph.readthedocs.io/en/latest/", None),
    matplotlib=("https://matplotlib.org/stable/", None),
    networkx=("https://networkx.org/documentation/stable/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    pynndescent=("https://pynndescent.readthedocs.io/en/latest/", None),
    pytest=("https://docs.pytest.org/en/latest/", None),
    python=("https://docs.python.org/3", None),
    rapids_singlecell=("https://rapids-singlecell.readthedocs.io/en/latest/", None),
    scanpy=("https://scanpy.readthedocs.io/en/stable/", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    tutorials=("https://scanpy-tutorials.readthedocs.io/en/latest/", None),
)