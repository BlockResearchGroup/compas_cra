
"""
********************************************************************************
datastructures
********************************************************************************

.. currentmodule:: compas_cra.datastructures


Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    CRA_Assembly


Functions
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    assembly_interfaces_numpy


"""


from __future__ import absolute_import

from .cra_assembly import *  # noqa: F401 F403
from .interfaces_numpy import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
