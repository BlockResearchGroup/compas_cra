.. _Installation:

********************************************************************************
Installation
********************************************************************************

Install with conda and pip (recommended)
========================================

Create an ``conda`` environment ``cra``.

.. code-block:: bash

    conda create -n cra python=3.8

Activate the environment.

.. code-block:: bash

    conda activate cra

Install ``compas_cra`` with ``pip``.

.. code-block:: bash

    pip install compas_cra

Install `IPOPT <https://coin-or.github.io/Ipopt/>`_ solver.

.. code-block:: bash

     conda install ipopt  # For Windows: conda install ipopt=3.11.1

Install `compas_view2 <https://compas.dev/compas_view2/>`_ for visualisation.

.. code-block:: bash

     conda install -c conda-forge compas_view2=0.7.0

Verify that the installation was successful.

.. code-block:: bash

    python -m compas_cra

.. code-block:: bash

    Yay! COMPAS CRA is installed correctly!

Try the :ref:`Tutorial` script.

.. code-block:: bash

    python scripts/tutorial_cubes.py


Developer Guide
===============

You can also install COMPAS CRA manually from `source <https://github.com/BlockResearchGroup/compas_cra>`_.

Install from source
-------------------

Create a virtual environment using your tool of choice (e.g. virtualenv, conda, etc), optional.

.. code-block:: bash

    conda create -n cra
    conda activate cra

Go to your directory and clone :mod:`compas_cra`:

.. code-block:: bash

    git clone git@github.com:BlockResearchGroup/compas_cra.git
    cd compas_cra

Install requirements

.. code-block:: bash

    pip install -r requirements.txt
    pip install -r requirements-dev.txt

In the `requirements-dev.txt <https://github.com/BlockResearchGroup/compas_cra/blob/main/requirements-dev.txt>`_, we also installed COMPAS CRA as an editable version from local source.

.. code-block:: bash

    pip install -e .

Install `IPOPT <https://coin-or.github.io/Ipopt/>`_ solver.

.. code-block:: bash

    conda install ipopt  # For Windows: conda install ipopt=3.11.1

Install `compas_view2 <https://compas.dev/compas_view2/>`_ for visualisation.

.. code-block:: bash

    conda install -c conda-forge compas_view2=0.7.0


A quicker way - from `.yml` file
--------------------------------

Set up conda channels

.. code-block:: bash

    conda config --add channels conda-forge


Clone :mod:`compas_cra`:

.. code-block:: bash

    git clone git@github.com:BlockResearchGroup/compas_cra.git
    cd compas_cra

Install COMPAS CRA with all dependencies in a new conda environment:

.. code-block:: bash

    conda env create -f env_osx.yml  # (Mac)
    conda env create -f env_win.yml  # (Windows)
    conda env create -f env_linux.yml  # (Linux)

    conda activate cra  # you can change the environment name in .yml file

Update conda packages
---------------------

.. code-block:: bash

    conda env update cra --file env_osx.yml --prune  # (Mac)
    conda env update cra --file env_win.yml --prune  # (Windows)
    conda env update cra --file env_linux.yml --prune  # (Linux)


Verify installation
-------------------

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
    -- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
    ======================== 4 passed, 5 warnings in 2.41s =========================


Rhino Installation
==================

:mod:`compas_cra` is developed independent of the functionality of CAD software.
However, CAD software can be useful to create geometrical objects.
For a more detailed information on how to install COMPAS and its packages for Rhino,
please refer to `Working in Rhino <https://compas.dev/compas/latest/gettingstarted/rhino.html>`_ page of the COMPAS documentation.

In order to install COMPAS CRA for Rhino, do

::

    $ python -m compas_rhino.uninstall
    $ python -m compas_rhino.install
    $ python -m compas_rhino.install -p compas_cra

Every time a new file is opened in Rhino, be sure to restart Rhino or reset the Python Script Engine before running scripts.


Verify Conda Environment, COMPAS CRA, and COMPAS in Rhino
---------------------------------------------------------

.. literalinclude:: ../scripts/rhinoenv.py
    :language: python

File can also be found in `scripts/rhinoenv.py <https://github.com/BlockResearchGroup/compas_cra/blob/main/scripts/rhinoenv.py>`_

Known issues (Windows)
======================

- pyomo cannot find ipopt location: ``pyomo.common.errors.ApplicationError: No executable found for solver 'ipopt'``. Please refer this thread for solution: https://groups.google.com/g/open-dsopf/c/wYPbZp-HLCw?pli=1

- :mod:`compas_cra` uses `IPOPT <https://coin-or.github.io/Ipopt/>`_ solver, so it might not work for PC with AMD processor.
