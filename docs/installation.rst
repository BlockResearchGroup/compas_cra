.. _Installation:

********************************************************************************
Installation
********************************************************************************

Stable
======

Stable releases can be installed on ``osx`` and ``linux`` using a combination of conda and pip.

.. code-block:: bash

    conda create -n cra -c conda-forge python=3.10 pyomo=6.4.2 ipopt=3.14.9 compas compas_viewer
    conda activate cra
    pip install compas_assembly compas_cra


Latest
======

The latest version can be installed on ``osx`` and ``linux`` from local source using an environment file.

.. code-block:: bash

    git clone https://github.com/blockresearchgroup/compas_cra.git
    cd compas_cra
    conda env create -f environment.yml
    conda activate cra-dev

Note that this will automatically create an editable install that can be used for development.


Windows
=======

On Windows, the first part of the installation is the same.

However, once you start using CRA, you will very likely run into a problem with missing ``ipopt`` binares.
To solve this, you can download the binaries from the ``ipopt`` github release and place it in in the correct folder of you environment.
(Thank you Petras for figuring this out :)

The bnary of release ``ipopt=3.14.9`` is available here:

* <https://github.com/coin-or/Ipopt/releases/download/releases%2F3.14.9/Ipopt-3.14.9-win64-msvs2019-md.zip>

Copy the contents of the ``bin`` folder of the archive into the ``Library\bin`` of your environment.
For example, into ``C:\Users\You\anaconda3\envs\cra\Library\bin``.

If this still doesn't fix the problem, please let us know.