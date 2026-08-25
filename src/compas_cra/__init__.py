import os

__author__ = ["Gene Ting-Chun Kao"]
__copyright__ = "Gene Ting-Chun Kao"
__license__ = "MIT License"
__email__ = "kao@arch.ethz.ch, kao.gene@gmail.com"
__version__ = "0.7.2"


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.abspath(os.path.join(HOME, "data"))
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))

SRC = os.path.abspath(os.path.join(HOME, "src"))
# Relative to the package itself, so the samples are found in an installed wheel too,
# not only in a source checkout.
CRA = os.path.abspath(HERE)
SAMPLE = os.path.abspath(os.path.join(CRA, "data", "samples"))

__all__ = ["HOME", "DATA", "DOCS", "TEMP"]
