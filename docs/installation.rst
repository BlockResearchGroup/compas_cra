********************************************************************************
Installation
********************************************************************************


Install with conda
==================

Create a new conda environment

.. code-block:: bash:

    conda config --add channels conda-forge
    conda create -n cra
    conda activate cra

Clone compas_cra:

.. code-block:: bash

    git clone git@github.com:BlockResearchGroup/compas_cra.git
    cd compas_cra

Install compas_cra with all dependencies:

.. code-block:: bash

    pip install -r requirements-dev.txt

Install ipopt solvers for pyomo:

.. code-block:: bash

    conda install -c conda-forge ipopt

Known issue (Windows)

pyomo cannot find ipopt location: ``pyomo.common.errors.ApplicationError: No executable found for solver 'ipopt'``.
Please refer this thread for solution: https://groups.google.com/g/open-dsopf/c/wYPbZp-HLCw?pli=1
