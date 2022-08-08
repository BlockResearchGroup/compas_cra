********************************************************************************
Installation
********************************************************************************


Install with conda
==================


Set up conda channels

.. code-block:: bash

    conda config --add channels conda-forge


Clone `compas_cra`:

.. code-block:: bash

    git clone git@github.com:BlockResearchGroup/compas_cra.git
    cd compas_cra

Install compas_cra with all dependencies in a new conda environment:

.. code-block:: bash

    conda env create -f env_osx.yml  # (Mac)
    conda env create -f env_win.yml  # (Windows)
    conda env create -f env_linux.yml  # (Linux)

    conda activate cra



Update conda packages

.. code-block:: bash

    conda env update cra --file env_osx.yml --prune  # (Mac)
    conda env update cra --file env_win.yml --prune  # (Windows)
    conda env update cra --file env_linux.yml --prune  # (Linux)

Known issues (Windows)

- pyomo cannot find ipopt location: ``pyomo.common.errors.ApplicationError: No executable found for solver 'ipopt'``. Please refer this thread for solution: https://groups.google.com/g/open-dsopf/c/wYPbZp-HLCw?pli=1

- `compas_cra` uses IPOPT solver, so it might not work for PC with AMD processor.
