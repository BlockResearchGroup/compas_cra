.. _Installation:

********************************************************************************
Installation
********************************************************************************

Stable
======

Stable releases can be installed using a combination of conda and pip.

.. code-block:: bash

    conda create -n cra python=3.9 pyomo=6.4.2 ipopt compas

.. note::

    On Windows, you should install ``ipopt=3.11.1``.

.. code-block:: bash

    conda activate cra
    pip install compas_cra

To use the CRA viewer, you should also install :mod:`compas_view2`
and the COMPAS 2 migration from the compatibility branch of the github repo.
(we will replace this by :mod:`compas_viewer` soon).

.. code-block:: bash

    conda install matplotlib compas_view2
    pip install git+https://github.com/compas-dev/compas_view2.git@compas2


Latest
======

The latest version can be installed from local source using an environment file.
Please use the correct environment file for your system
(``env_linux.yml``, ``env_osx.yml``, ``env_win.yml``)

.. code-block:: bash

    git clone https://github.com/blockresearchgroup/compas_cra.git
    cd compas_cra
    conda env create -f env_osx.yml

Note that this will automatically create an editable install that can be used for development.
