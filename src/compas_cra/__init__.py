"""
********************************************************************************
compas_cra
********************************************************************************

.. currentmodule:: compas_cra


.. toctree::
    :maxdepth: 1

    compas_cra.datastructures
    compas_cra.algorithms
    compas_cra.equilibrium
    compas_cra.geometry
    compas_cra.viewers

"""

from __future__ import print_function

import os

__author__ = ["Gene Ting-Chun Kao"]
__copyright__ = "Gene Ting-Chun Kao"
__license__ = "MIT License"
__email__ = "kao@arch.ethz.ch, kao.gene@gmail.com"
__version__ = "0.2.0"


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.abspath(os.path.join(HOME, "data"))
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))

SRC = os.path.abspath(os.path.join(HOME, "src"))
CRA = os.path.abspath(os.path.join(SRC, "compas_cra"))
SAMPLE = os.path.abspath(os.path.join(CRA, "data", "samples"))

__all__ = ["HOME", "DATA", "DOCS", "TEMP"]
