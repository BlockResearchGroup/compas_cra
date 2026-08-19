.. _Installation:

********************************************************************************
Installation
********************************************************************************

Stable
======

.. code-block:: bash

    pip install compas_cra

That is all that is needed, on Windows, macOS (Apple Silicon and Intel) and Linux.
The wheels bundle a statically linked `IPOPT <https://coin-or.github.io/Ipopt/>`_
executable, so no conda environment, no homebrew and no manual download of solver
binaries are required.

To also install the viewers:

.. code-block:: bash

    pip install compas_cra[viz]

Verify the solver is available with:

.. code-block:: bash

    compas-cra-ipopt --version


Latest
======

The latest version can be installed from local source.

.. code-block:: bash

    git clone https://github.com/blockresearchgroup/compas_cra.git
    cd compas_cra
    pip install -e ".[dev]"

A source install contains no IPOPT binary, since that is added to the wheels at build
time. Either build one with ``packaging/build_ipopt.sh`` (see ``packaging/README.md``),
or provide ``ipopt`` yourself on the ``PATH`` - which is what the conda development
environment does:

.. code-block:: bash

    conda env create -f environment.yml
    conda activate cra-dev


Conda
=====

Using conda for the solver is still supported: if no bundled binary is present,
``compas_cra`` falls back to whatever ``ipopt`` is on the ``PATH``.

.. code-block:: bash

    conda create -n cra -c conda-forge python=3.10 ipopt compas compas_viewer
    conda activate cra
    pip install compas_assembly compas_cra
