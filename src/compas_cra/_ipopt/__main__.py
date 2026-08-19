"""Run the IPOPT executable that :mod:`compas_cra` would solve with.

Installed as the ``compas-cra-ipopt`` console script, and reachable as
``python -m compas_cra._ipopt``. Useful for checking which solver is being used.

The command is deliberately not called ``ipopt``: pip would install it over a conda
environment's own ``bin/ipopt``, and :func:`compas_cra._ipopt.executable` would then
find this wrapper instead of a real solver.
"""

import os
import subprocess
import sys

from compas_cra._ipopt import executable


def main():
    """Forward all arguments to the IPOPT executable."""
    exe = executable()
    if exe is None:
        sys.stderr.write(
            "no IPOPT executable is bundled with this installation of compas_cra, and none was found on PATH\n"
        )
        return 1
    args = [str(exe)] + sys.argv[1:]
    if os.name == "posix":
        os.execv(str(exe), args)  # replaces this process, never returns
    return subprocess.call(args)


if __name__ == "__main__":
    sys.exit(main())
