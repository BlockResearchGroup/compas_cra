"""Locate the IPOPT executable used to solve CRA equilibrium problems.

COMPAS CRA solves its models with :mod:`pyomo`, which drives IPOPT through the AMPL
``.nl`` file interface: it shells out to the ``ipopt`` command line executable rather
than linking against the library. There are no IPOPT wheels on PyPI (``cyipopt``
publishes an sdist only, and needs a system IPOPT to compile against), so the platform
wheels of :mod:`compas_cra` ship a statically linked ``ipopt`` binary inside this
package. See ``packaging/build_ipopt.sh`` for how it is built.

The source distribution contains no binary; there :func:`executable` falls back to an
``ipopt`` on ``PATH``, e.g. one installed with conda.
"""

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

__all__ = ["IPOPT_VERSION", "bundled", "executable", "ipopt_version"]

#: Version of IPOPT that the bundled binaries are built from.
IPOPT_VERSION = "3.14.19"


def _binary_name():
    return "ipopt.exe" if sys.platform == "win32" else "ipopt"


def bundled():
    """Return the path of the IPOPT executable shipped in this package.

    Returns
    -------
    pathlib.Path or None
        None if this is an installation without a bundled binary (e.g. from the sdist,
        or on a platform we do not build wheels for).

    """
    path = Path(__file__).parent / "bin" / _binary_name()
    if not path.is_file():
        return None
    if os.name == "posix" and not os.access(path, os.X_OK):
        # Wheels do not carry the executable bit on package data, only on scripts.
        try:
            path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        except OSError:
            return None
    return path


def executable():
    """Return the path of the IPOPT executable to solve with.

    The binary bundled in this package is preferred over one on ``PATH``, so that the
    solver is the version this release was tested against regardless of what else the
    environment provides.

    Returns
    -------
    pathlib.Path or None
        None if no IPOPT executable can be found at all.

    """
    path = bundled()
    if path is not None:
        return path
    found = shutil.which("ipopt")
    return Path(found) if found else None


def ipopt_version():
    """Return the version reported by the IPOPT executable itself.

    Returns
    -------
    str or None
        None if no executable is available or it cannot be run.

    """
    exe = executable()
    if exe is None:
        return None
    try:
        out = subprocess.run([str(exe), "--version"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    line = (out.stdout or out.stderr).strip().splitlines()
    return line[0].strip() if line else None
