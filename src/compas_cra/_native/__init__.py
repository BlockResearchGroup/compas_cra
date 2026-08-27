"""In-process IPOPT solver, compiled into compas_cra.

The heavy lifting happens in the ``_core`` extension module, which links IPOPT and the
MUMPS linear solver statically. It is built from ``native/src/ipopt_nb.cpp`` by the
CMake build at the repository root, so it is present in the platform wheels and absent
from a plain source checkout. See :mod:`compas_cra.nlp` for the high-level API.
"""

import os

try:
    from ._core import IPOPT_VERSION
    from ._core import solve_nlp
except ImportError:
    # A platform wheel never gets here: delvewheel grafts the runtime DLLs into it. An
    # editable install on Windows does - the extension links MSYS2's shared runtimes
    # (libgfortran, libopenblas, libstdc++, ...) out of ucrt64\bin, and Python does not
    # search PATH for extension DLLs. Point the loader at them and try once more.
    if os.name != "nt":
        raise
    _candidates = [
        os.path.join(_root, "ucrt64", "bin")
        for _root in (os.environ.get("MSYS2_ROOT"), "C:\\msys64", "D:\\msys64")
        if _root
    ]
    _found = [_d for _d in _candidates if os.path.isdir(_d)]
    if not _found:
        raise
    for _d in _found:
        os.add_dll_directory(_d)
    from ._core import IPOPT_VERSION
    from ._core import solve_nlp

__all__ = ["solve_nlp", "IPOPT_VERSION"]
