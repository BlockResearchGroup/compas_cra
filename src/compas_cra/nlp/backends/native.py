"""In-process IPOPT backend, provided by the :mod:`compas_cra._native` extension.

The extension links IPOPT (with the MUMPS linear solver) statically into a Python
module, so solving happens in-process with numpy-array callbacks: no subprocess, no
``.nl`` files, no bundled executable.
"""

import numpy as np

from ..problem import NLPResult

__all__ = ["is_available", "solve"]

# IPOPT ApplicationReturnStatus -> NLPResult status. Feasible_Point_Found counts as
# acceptable, matching the pyomo path (its result check passes both optimal and
# feasible termination conditions).
_STATUS = {
    0: NLPResult.OPTIMAL,  # Solve_Succeeded
    1: NLPResult.ACCEPTABLE,  # Solved_To_Acceptable_Level
    2: NLPResult.INFEASIBLE,  # Infeasible_Problem_Detected
    6: NLPResult.ACCEPTABLE,  # Feasible_Point_Found
}

_DEFAULT_OPTIONS = {
    "sb": "yes",  # suppress the ipopt banner
}


def is_available():
    try:
        from ... import _native  # noqa: F401
    except ImportError:
        return False
    return True


def solve(problem, options=None, verbose=False):
    """Solve an :class:`~compas_cra.nlp.problem.NLPProblem` with in-process IPOPT."""
    from ... import _native

    opts = dict(_DEFAULT_OPTIONS)
    opts.update(options or {})
    if "print_level" not in opts:
        opts["print_level"] = 5 if verbose else 0
    if not problem.has_hessian:
        opts.setdefault("hessian_approximation", "limited-memory")

    raw = _native.solve_nlp(
        n=problem.n,
        m=problem.m,
        x_l=problem.x_l,
        x_u=problem.x_u,
        g_l=problem.g_l,
        g_u=problem.g_u,
        x0=problem.x0,
        jac_rows=problem.jac_rows,
        jac_cols=problem.jac_cols,
        eval_f=lambda x: float(problem.objective(x)),
        eval_grad=lambda x: np.ascontiguousarray(problem.gradient(x), dtype=float),
        eval_g=lambda x: np.ascontiguousarray(problem.constraints(x), dtype=float),
        eval_jac=lambda x: np.ascontiguousarray(problem.jacobian(x), dtype=float),
        hess_rows=problem.hess_rows,
        hess_cols=problem.hess_cols,
        eval_h=(
            None
            if not problem.has_hessian
            else lambda x, sigma, lam: np.ascontiguousarray(problem.hessian(x, sigma, lam), dtype=float)
        ),
        options=opts,
    )

    status = _STATUS.get(raw["status"], NLPResult.FAILED)
    message = raw["status_name"]

    # On degenerate problems (the CRA complementarity constraints are exactly that)
    # IPOPT's restoration phase can converge to a fully feasible, essentially optimal
    # point and still return Restoration_Failed because the filter refuses the step
    # back into regular mode. Judge the returned point on its merits: if it satisfies
    # the acceptable-level tolerances, it is a solution. Maximum_Iterations_Exceeded is
    # deliberately NOT in this list: an iteration-capped point is wherever the solver
    # happened to stop, feasible or not, and the executable path fails there too.
    if status is NLPResult.FAILED and message in (
        "Restoration_Failed",
        "Search_Direction_Becomes_Too_Small",
    ):
        viol = _point_violation(problem, raw["x"], raw["g"])
        tol = float(opts.get("acceptable_constr_viol_tol", 1e-6))
        if viol <= tol:
            status = NLPResult.ACCEPTABLE
            message = "{} (feasible point accepted, violation {:.2e})".format(message, viol)

    return NLPResult(
        x=raw["x"],
        obj=raw["obj"],
        status=status,
        status_message=message,
        g=raw["g"],
        mult_g=raw["mult_g"],
        iterations=raw.get("iterations"),
    )


def _point_violation(problem, x, g):
    """Worst constraint/bound violation of a candidate point."""
    x = np.asarray(x, dtype=float)
    g = np.asarray(g, dtype=float) if g is not None and len(g) else problem.constraints(x)
    viol = 0.0
    if problem.m:
        viol = max(viol, float(np.max(problem.g_l - g, initial=0.0)))
        viol = max(viol, float(np.max(g - problem.g_u, initial=0.0)))
    viol = max(viol, float(np.max(problem.x_l - x, initial=0.0)))
    viol = max(viol, float(np.max(x - problem.x_u, initial=0.0)))
    return viol
