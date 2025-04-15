========
Overview
========

Python helpers at Singular Genomics.

Installation
************
.. code-block:: bash

    ## Clone the repo
    git clone git@bitbucket.org:singulargenomics/pysing.git

    ## Navigate to repo
    cd pySing

    ## To install basic dependencies (covers brightphase and utils)
    pip install .

    ## To install dependencies for seq module
    pip install .[seq]

    ## To install dependencies for scseq module
    pip install .[scseq]

    ## To install all dependencies
    pip install .[seq,scseq]


Documentation
*************

Full documentation located `here <http://hpc-05.sg.local:3000/pysing/main/latest/sphinx_docs/_build/html/index.html>`_


The ``utils`` module contains helper functions, dataclasses, and standalone modules, some of which have CLI entrypoints as shown below.


The ``seq`` module contains functions related to fastq processing and other common NGS applications.


The ``scseq`` module contains functions related to single-cell-sequencing processing.


The ``brightphase`` module was lifted from `smithii <https://bitbucket.org/singulargenomics/smithii/src/master/>`_ on 20230712 and modified to make it useable outside of the smithii framework/environment.


``utils``
=========

- ``plotting.py`` (module containing helpful plotting functions, especially fix_plot)

- ``filter_pair.py`` (CLI: ``filter_pair -h`` - module to filter and pair fastqs, suitable for parallelization)

- ``metric_aggregator.py`` (CLI: ``aggregate_metrics -h`` - module that attempts to abstract aggregation of metrics from granular level to coarse level)

- ``flowcell_plotter.py`` (CLI: ``flowcell_plot -h`` - module that creates flow cell montages from various input data types with optional metrics)

- ``parse_secondary.py`` (CLI: ``parse_secondary -h`` - deprecated since Acc2, useful for scraping Acc1 results and making merged-run plots)

- ``brightphase_runner.py`` (CLI: ``brightphase_runner -h`` - launch brightphase on a directory of h5s)

``seq``
=======

- ``rawdata.py`` (fastq and h5 processing)

- ``post_aln.py`` (bam and vcf processing)

``scseq``
=========

- ``generalfunctions.py``

- ``qcfunctions.py``

- ``dgexfunctions.py``

- ``populationfunctions.py`` (work in progress, no guarantees)

- ``airrfunctions.py`` (work in progress, no guarantees)

``brightphase``
===============

- see CLI entrypoint in ``utils`` above