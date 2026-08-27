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
# The CRA system is square (as many equations as unknowns), and that made IPOPT's
# version the hidden variable of this solver's history. Until IPOPT 3.14.11, square
# problems were special-cased: the convergence check ignored dual feasibility and bound
# complementarity, so IPOPT stopped at the first primal-feasible point of the barrier
# path - an interior point, where every contact face still carries force. Every result
# in the CRA paper (Kao et al. 2022, doi:10.1016/j.cad.2022.103216) and every published
# example screenshot is such a point: on the curved-interface examples the load spreads
# over all 72 sub-faces. IPOPT 3.14.12 removed the special case, and from there the
# solver drives complementarity to zero and lands on a degenerate vertex instead - the
# same examples concentrate the whole load on 8 faces and zero the rest. Verified by
# building IPOPT 3.14.9, 3.14.11 and 3.14.14 from source against the same MUMPS 5.9:
# 3.14.11 still spreads, 3.14.14 concentrates - the flip is exactly at 3.14.12, per
# its changelog; not MUMPS, and not this repository - every code generation from the
# screenshot-day commit (15e5edc) to today produces both solutions depending only on
# the binary.
#
# mu_target restores the published behavior as an explicit setting rather than an
# accident of an old binary: the barrier stops at mu = 5e-6 instead of 0, which is the
# interior point the paper's figures show. Calibrated against IPOPT 3.14.9 on the
# curved examples (cube-curve-short: max/median resultant 0.00057/0.00050 vs the era's
# 0.00066/0.00041; cube-curve-tall: 0.00191/0.00050 vs 0.00179/0.00058; all 72 faces
# loaded in every case). Equilibrium stays exact - constr_viol_tol keeps the force
# balance at 1e-8 - only the force *distribution* on statically indeterminate contacts
# is selected. The tolerances sit just above the barrier target, which mu_target
# requires, and they are load-bearing: the degenerate examples are knife-edged in both
# mu_target and tol (curve-3-blocks solves at this pair and exhausts even 9000
# iterations one small step away), so change them together or not at all. The arch
# converges to the same resultants as before (1.9570 / 0.8430) with iteration headroom
# intact (810 of 9000).
_CRA_OPTIONS = {
    "tol": 1e-5,
    "dual_inf_tol": 1e-5,
    "compl_inf_tol": 1e-5,
    "constr_viol_tol": 1e-8,
    "acceptable_tol": 1e-4,
    "acceptable_constr_viol_tol": 1e-7,
    "acceptable_compl_inf_tol": 1e-4,
    # CRA's complementarity constraints are degenerate, and the monotone barrier update
    # crawls on them; the adaptive update reaches the same points in a fraction of the
    # iterations. See the note above for why mu_target is the load-distribution knob.
    "mu_strategy": "adaptive",
    "mu_target": 5e-6,
    # the degenerate curved examples legitimately need a few thousand iterations in
    # this regime; the default 3000 cap is calibrated for strict optimization
    "max_iter": 9000,
}

# the tolerances the CRA penalty solver has always used
_CRA_PENALTY_OPTIONS = {
    "tol": 1e-8,
    "constr_viol_tol": 1e-7,
    "acceptable_tol": 1e-6,
    "acceptable_constr_viol_tol": 1e-5,
    # same degenerate complementarity constraints as the CRA solver, same crawl under
    # the monotone update: the 20-block arch exhausted the 3000-iteration cap at 5e-2
    # violation. See _CRA_OPTIONS.
    "mu_strategy": "adaptive",
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
