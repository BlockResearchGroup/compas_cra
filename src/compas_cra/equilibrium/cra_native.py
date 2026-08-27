"""CRA, CRA-penalty and RBE solvers on the in-process IPOPT binding.

No pyomo, no subprocess: each solver builds its problem with
:mod:`~compas_cra.equilibrium.cra_nlp` and solves through the
``compas_cra._native`` extension, writing the same results to the assembly as the
historical pyomo + executable path (per-interface ``interface.forces`` and a
``displacement`` attribute per free node).
"""

import time

from compas_assembly.datastructures import Assembly

from ..nlp import solve_nlp
from .cra_nlp import cra_penalty_problem
from .cra_nlp import cra_problem
from .cra_nlp import rbe_problem

__all__ = ["cra_solve_native", "cra_penalty_solve_native", "rbe_solve_native"]

# the tolerances the CRA solver has always used. They look brutally tight for double
# precision, and they are, but they are load bearing: the solve reaches its answer
# through IPOPT's restoration phase, and relaxing them keeps it out of restoration
# without getting it any closer to convergence. Measured on the 20-block arch,
# tol=1e-8 / constr_viol_tol=1e-9 turns a solve that stops at a feasible point into
# Maximum_Iterations_Exceeded at the 3000-iteration cap.
_CRA_OPTIONS = {
    "tol": 1e-10,
    "constr_viol_tol": 1e-12,
    "compl_inf_tol": 1e-12,
    "acceptable_tol": 1e-8,
    "acceptable_constr_viol_tol": 1e-8,
    "acceptable_compl_inf_tol": 1e-8,
    # CRA's complementarity constraints are degenerate, and the monotone barrier update
    # crawls on them: the 20-block arch takes 1964 of its 3000 permitted iterations,
    # which leaves almost no headroom before a wheel that rounds slightly differently
    # runs out of iterations instead -- which is exactly what users on macOS hit, where
    # the wheel links Accelerate rather than OpenBLAS. The adaptive update finishes the
    # same solve in 668 iterations, on the same answer: interface resultants agree to
    # 1e-10, and a 40-block arch goes from 2757 iterations to 307.
    "mu_strategy": "adaptive",
}

# the tolerances the CRA penalty solver has always used
_CRA_PENALTY_OPTIONS = {
    "tol": 1e-8,
    "constr_viol_tol": 1e-7,
    "acceptable_tol": 1e-6,
    "acceptable_constr_viol_tol": 1e-5,
}


def cra_solve_native(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.0,
    d_bnd: float = 1e-3,
    eps: float = 1e-4,
    verbose: bool = False,
    timer: bool = False,
) -> Assembly:
    """CRA solver (paper Eq. 11) using the in-process IPOPT binding."""
    problem, layout = _timed(
        timer, lambda: cra_problem(assembly, mu=mu, density=density, d_bnd=d_bnd, eps=eps), "set up"
    )
    result = _timed(
        timer, lambda: solve_nlp(problem, backend="native", options=_CRA_OPTIONS, verbose=verbose), "solving"
    )
    _check(result)
    _result_to_assembly(result.x, layout, assembly, shift=3, verbose=verbose)
    return assembly


def cra_penalty_solve_native(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.0,
    d_bnd: float = 1e-3,
    eps: float = 1e-4,
    verbose: bool = False,
    timer: bool = False,
) -> Assembly:
    """CRA penalty solver (paper Eq. 14) using the in-process IPOPT binding."""
    problem, layout = _timed(
        timer, lambda: cra_penalty_problem(assembly, mu=mu, density=density, d_bnd=d_bnd, eps=eps), "set up"
    )
    result = _timed(
        timer, lambda: solve_nlp(problem, backend="native", options=_CRA_PENALTY_OPTIONS, verbose=verbose), "solving"
    )
    _check(result)
    _result_to_assembly(result.x, layout, assembly, shift=4, verbose=verbose)
    return assembly


def rbe_solve_native(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.0,
    verbose: bool = False,
    timer: bool = False,
) -> Assembly:
    """RBE solver (paper Eq. 6) using the in-process IPOPT binding."""
    problem, layout = _timed(timer, lambda: rbe_problem(assembly, mu=mu, density=density), "set up")
    result = _timed(timer, lambda: solve_nlp(problem, backend="native", verbose=verbose), "solving")
    _check(result)
    _result_to_assembly(result.x, layout, assembly, shift=4, verbose=verbose)
    return assembly


def _timed(timer, fn, label):
    start = time.time()
    out = fn()
    if timer:
        print("--- %s time: %s seconds ---" % (label, time.time() - start))
    return out


def _check(result):
    if not result.success:
        raise ValueError("solve failed: {} ({})".format(result.status, result.status_message))
    print("result: ", result.status)


def _result_to_assembly(x, layout, assembly, shift, verbose=False):
    """Write the solution vector to the assembly, mirroring pyomo_result_assembly."""
    f = x[layout["f"]]

    offset = 0
    for edge in assembly.graph.edges():
        for interface in assembly.graph.edge_attribute(edge, "interfaces"):
            interface.forces = []
            for i in range(len(interface.points)):
                base = offset + shift * i
                interface.forces.append(
                    {
                        "c_np": float(f[base + 0]),
                        "c_nn": float(f[base + 1]) if shift == 4 else 0,
                        "c_u": float(f[base + shift - 2]),
                        "c_v": float(f[base + shift - 1]),
                    }
                )
            offset += shift * len(interface.points)

    if layout["q"] is None:
        return
    q = x[layout["q"]]
    if verbose:
        print("q:", list(q))
    offset = 0
    for node in assembly.graph.nodes():
        if assembly.graph.node_attribute(node, "is_support"):
            continue
        assembly.graph.node_attribute(node, "displacement", [float(v) for v in q[offset : offset + 6]])
        offset += 6
