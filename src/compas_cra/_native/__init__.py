"""In-process IPOPT solver, compiled into compas_cra.

The heavy lifting happens in the ``_core`` extension module, which links IPOPT and the
MUMPS linear solver statically. It is built from ``native/src/ipopt_nb.cpp`` by the
CMake build at the repository root, so it is present in the platform wheels and absent
from a plain source checkout. See :mod:`compas_cra.nlp` for the high-level API.
"""

from ._core import IPOPT_VERSION
from ._core import solve_nlp

__all__ = ["solve_nlp", "IPOPT_VERSION"]
