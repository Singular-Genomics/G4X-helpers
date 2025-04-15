Installation
============

Docker
------

Docker images that contain ``pySing`` and all associated dependencies are available through AWS ECR.

AWS ECR

   ``docker pull 448313343180.dkr.ecr.us-west-1.amazonaws.com/pysing_full:<version>``
   ``docker pull 448313343180.dkr.ecr.us-west-1.amazonaws.com/pysing_basic:<version>``

Source installation / CLI usage
-------------------------------

``pySing`` can be installed and run directly as a python package

1. Prepare a python environment

   - Install ``conda``, ``miniconda`` or ``mamba``
   - Create the environment: ``conda create -n pysing_env python=3.10``
   - Activate the environment: ``conda activate pysing_env``

2. Clone and install ``pySing``

   - Clone `pySing <git@bitbucket.org:singulargenomics/pysing.git>`_
   - ``cd pySing``
   - ``pip install .`` or ``pip install -e .`` for an editable install
   - Extra modules that can be installed include: ``seq``, ``scseq``, ``spatial``, ``gpu``
   - ``pip install .[seq,scseq,spatial]``
   - If installing the ``gpu`` module, run like below:
   - ``pip install --extra-index-url=https://pypi.nvidia.com .[gpu]``