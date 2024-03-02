from __future__ import print_function

import os

__author__ = ["Gene Ting-Chun Kao"]
__copyright__ = "Gene Ting-Chun Kao"
__license__ = "MIT License"
__email__ = "kao@arch.ethz.ch, kao.gene@gmail.com"
__version__ = "0.3.0"


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.abspath(os.path.join(HOME, "data"))
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))

SRC = os.path.abspath(os.path.join(HOME, "src"))
CRA = os.path.abspath(os.path.join(SRC, "compas_cra"))
SAMPLE = os.path.abspath(os.path.join(CRA, "data", "samples"))

__all__ = ["HOME", "DATA", "DOCS", "TEMP"]
