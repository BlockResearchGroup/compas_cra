.. _Installation:

********************************************************************************
Installation
********************************************************************************

Stable
======

Stable releases can be installed using a combination of conda and pip.

.. code-block:: bash

    conda create -n cra python=3.10 pyomo=6.4.2 ipopt=3.14.9 compas

.. note::

    On Windows, you may have to install ``ipopt=3.11.1``.

.. code-block:: bash

    conda activate cra
    pip install compas_assembly compas_cra

To use the CRA viewer, you should also install :mod:`compas_viewer`.

.. code-block:: bash

    conda install compas_viewer


Latest
======

The latest version can be installed from local source using an environment file.
Please use the correct environment file for your system
(``env_linux.yml``, ``env_osx.yml``, ``env_win.yml``)

.. code-block:: bash

    git clone https://github.com/blockresearchgroup/compas_cra.git
    cd compas_cra
    conda env create -f env_osx.yml
    conda activate cra-dev

Note that this will automatically create an editable install that can be used for development.
