"""Construction of the IPOPT solver used by all CRA formulations."""

import pyomo.environ as pyo

from compas_cra._ipopt import executable


def ipopt_solver(**options):
    """Create a pyomo IPOPT solver.

    The IPOPT executable bundled with :mod:`compas_cra` is used when it is available,
    so that no conda environment or manually installed solver is needed. Otherwise the
    solver is resolved from ``PATH`` as before.

    Parameters
    ----------
    **options : dict, optional
        IPOPT options, e.g. ``tol=1e-10``.
        See https://coin-or.github.io/Ipopt/OPTIONS.html

    Returns
    -------
    :class:`pyomo.opt.base.solvers.OptSolver`

    Raises
    ------
    RuntimeError
        If no IPOPT executable can be found.

    """
    exe = executable()
    if exe is None:
        raise RuntimeError(
            "No IPOPT executable found. This installation of compas_cra does not "
            "bundle one (source installs and unsupported platforms do not), and none "
            "is available on PATH. Install a binary wheel with `pip install compas_cra` "
            "or provide ipopt yourself, e.g. `conda install -c conda-forge ipopt`."
        )
    solver = pyo.SolverFactory("ipopt", executable=str(exe))
    for key, value in options.items():
        solver.options[key] = value
    return solver
