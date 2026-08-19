"""Run the bundled IPOPT executable.

Installed as the ``ipopt`` console script, so that ``SolverFactory("ipopt")`` finds a
solver on ``PATH`` in any environment where :mod:`compas_cra` is installed, and also
reachable as ``python -m compas_cra._ipopt``.
"""

import os
import subprocess
import sys

from compas_cra._ipopt import bundled


def main():
    """Forward all arguments to the bundled IPOPT executable."""
    exe = bundled()
    if exe is None:
        sys.stderr.write("no IPOPT executable is bundled with this installation of compas_cra\n")
        return 1
    args = [str(exe)] + sys.argv[1:]
    if os.name == "posix":
        os.execv(str(exe), args)  # replaces this process, never returns
    return subprocess.call(args)


if __name__ == "__main__":
    sys.exit(main())
