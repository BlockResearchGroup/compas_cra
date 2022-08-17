********************************************************************************
Installation
********************************************************************************


Install with conda
==================
.. _Installation:

Set up conda channels

.. code-block:: bash

    conda config --add channels conda-forge


Clone :mod:`compas_cra`:

.. code-block:: bash

    git clone git@github.com:BlockResearchGroup/compas_cra.git
    cd compas_cra

Install :mod:`compas_cra` with all dependencies in a new conda environment:

.. code-block:: bash

    conda env create -f env_osx.yml  # (Mac)
    conda env create -f env_win.yml  # (Windows)
    conda env create -f env_linux.yml  # (Linux)

    conda activate cra


Verify installation
===================

After running:

.. code-block:: bash

    invoke test

You should see something like:

.. code-block:: bash

    ============================= test session starts ==============================
    platform darwin -- Python 3.8.13, pytest-7.0.1, pluggy-1.0.0
    rootdir: ~/compas-dev/compas_cra, configfile: setup.cfg, testpaths: tests
    collected 4 items

    tests/test_cra.py .                                                      [ 25%]
    tests/test_cra_penalty.py .                                              [ 50%]
    tests/test_ipopt.py .                                                    [ 75%]
    tests/test_rbe.py .                                                      [100%]

    =============================== warnings summary ===============================
    .
    .
    .
    .
    .
    -- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
    ======================== 4 passed, 5 warnings in 2.41s =========================

Update conda packages
=====================

.. code-block:: bash

    conda env update cra --file env_osx.yml --prune  # (Mac)
    conda env update cra --file env_win.yml --prune  # (Windows)
    conda env update cra --file env_linux.yml --prune  # (Linux)

Known issues (Windows)
======================

- pyomo cannot find ipopt location: ``pyomo.common.errors.ApplicationError: No executable found for solver 'ipopt'``. Please refer this thread for solution: https://groups.google.com/g/open-dsopf/c/wYPbZp-HLCw?pli=1

- :mod:`compas_cra` uses `IPOPT <https://coin-or.github.io/Ipopt/>`_ solver, so it might not work for PC with AMD processor.
