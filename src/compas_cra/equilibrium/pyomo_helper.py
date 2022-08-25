"""Some functions to help building pyomo optimisation problems"""

from typing import Literal
from typing import Callable

import pyomo.environ as pyo
from numpy import zeros
from pyomo.core.base.matrix_constraint import MatrixConstraint


def initialisations(
    variable: Literal["f_tilde"],
) -> Callable:
    """Variable initialisations for pyomo.

    Parameters
    ----------
    variable : str
        * f_tilde: force, :math:`f ̃ = ({f_n}^+, {f_n}^-, f_u, f_v)`

    Returns
    -------
    Callable
        initialisations function for pyomo

    """

    def f_tilde_init(model, i):
        """initialise f ̃ with [1, 0, 1, 1]"""
        if i % 4 == 1:
            return 0.0
        return 1.0

    if variable == "f_tilde":
        return f_tilde_init


def bounds(
    variable: Literal["d", "f", "f_tilde"],
    d_bnd: float = 1e-3,
) -> Callable:
    r"""Variable bounds for pyomo.

    Parameters
    ----------
    variable : str
        * d: virtual displacement :math:`\delta d`
        * f: force, :math:`f = (f_n, f_u, f_v)`
        * f_tilde: force, :math:`f ̃ = ({f_n}^+, {f_n}^-, f_u, f_v)`
    d_bnd : float, optional
        displacement bounds, -d_bnd <= d <= d_bnd

    Returns
    -------
    Callable
        bounds constraint/domain function for pyomo

    """

    def f_tilde_bnds(model, i):
        """bounds of f ̃, f ̃ include [fn+, fn-, fu, fv]"""
        if i % 4 == 0 or i % 4 == 1:
            return pyo.NonNegativeReals
        return pyo.Reals

    def f_bnds(model, i):
        """bounds of f, f include [fn, fu, fv]"""
        if i % 3 == 0:
            return pyo.NonNegativeReals
        return pyo.Reals

    def d_bnds(model, i):
        return (-d_bnd, model.d[i], d_bnd)

    if variable == "f":
        return f_bnds
    if variable == "f_tilde":
        return f_tilde_bnds
    if variable == "d":
        return d_bnds


def objectives(
    solver: Literal["cra", "cra_penalty", "rbe"],
    weights: tuple = (1e0, 1e0, 1e6, 1e0),
) -> Callable:
    """Objective functions for pyomo.

    Parameters
    ----------
    solver : str
        * cra: CRA objective, :math:`W_{compression} * ||f_n||_2^2 + W_{α} * ||α||_2^2`
        * cra_penalty: CRA penalty objective, :math:`W_{compression} * ||{f_n}^+||_2^2 + W_{tension} * ||{f_n}^-||_2^2 + W_{α} * ||α||_2^2`
        * rbe: RBE objective, :math:`W_{compression} * ||{f_n}^+||_2^2 + W_{tension} * ||{f_n}^-||_2^2 + W_{friction} * ||{f_u}||_2^2 + W_{friction} * ||{f_v}||_2^2`
    weights : tuple, optional
        weighting factors, :math:`(W_{α}, W_{compression}, W_{tension}, W_{friction})`

    Returns
    -------
    Callable
        objective function for pyomo

    Notes
    -----

    dsodijsio jdois jois

    """

    def obj_rbe(model):
        """RBE objective function"""
        return _obj_weights(model)

    def obj_cra(model):
        """CRA objective function"""
        alpha_sum = pyo.dot_product(model.alpha, model.alpha)
        f_sum = 0
        for i in model.f_id:
            if i % 3 == 0:
                f_sum = f_sum + (model.f[i] * model.f[i])
        return f_sum + alpha_sum

    def obj_cra_penalty(model):
        """CRA penalty objective function"""
        alpha_sum = pyo.dot_product(model.alpha, model.alpha) * weights[0]  # alpha
        f_sum = _obj_weights(model)
        return alpha_sum + f_sum

    def _obj_weights(model):
        f_sum = 0
        for i in model.f_id:
            if i % 4 == 1:
                f_sum = f_sum + (model.f[i] * model.f[i] * weights[2])  # tension
            elif i % 4 == 0:
                f_sum = f_sum + (model.f[i] * model.f[i] * weights[1])  # compression
            elif i % 4 == 2 or i % 3 == 0:
                f_sum = f_sum + (model.f[i] * model.f[i] * weights[3])  # friction
        return f_sum

    if solver == "cra":
        return obj_cra
    if solver == "cra_penalty":
        return obj_cra_penalty
    if solver == "rbe":
        return obj_rbe


def constraints(
    name: Literal[
        "contact",
        "penalty_contact",
        "fn_np",
        "no_penetration",
        "ft_dt",
        "penalty_ft_dt",
    ],
    eps: float = 1e-4,
) -> Callable:
    r"""Constraint functions for pyomo.

    Parameters
    ----------
    name : str
        * contact: contact constraint, :math:`{f_{jkn}^i}\: ({\delta d_{jkn}^i} + eps) = 0`
        * penalty_contact: penalty formulation contact constraint, :math:`{f_{jkn}^{i+}}\:({\delta d_{jkn}^i} + eps) = 0`
        * fn_np: fn+ and fn- cannot coexist, :math:`{f_{jkn}^{i+}} \: {f_{jkn}^{i-}} = 0`
        * no_penetration: no penetration constraint, :math:`{f_{jkn}^{i+}}\:({\delta d_{jkn}^i} + eps) = 0`
        * ft_dt: friction and virtual sliding alignment, :math:`f_{jkt}^{i} = -{α_{jk}^i} \: \delta{d}_{jkt}^{i}`
        * penalty_ft_dt: penalty formulation friction and virtual sliding alignment, :math:`f_{jkt}^{i} = -{α_{jk}^i} \: \delta{d}_{jkt}^{i}`
    eps : float, optional
        epsilon, overlapping parameter

    Returns
    -------
    Callable
         constraint function for pyomo

    """

    def contact_con(model, i):
        """contact constraint"""
        dn = model.d[i * 3]
        fn = model.f[i * 3]
        return ((dn + eps) * fn, 0)

    def penalty_contact_con(model, i):
        """penalty formulation contact constraint"""
        dn = model.d[i * 3]
        fn = model.f[i * 4]
        return ((dn + eps) * fn, 0)

    def fn_np_con(model, i):
        """fn+ and fn- cannot coexist constraints"""
        return (model.f[i * 4] * model.f[i * 4 + 1], 0)

    def no_penetration_con(m, t):
        """no penetration constraint"""
        return (0, m.d[t * 3] + eps, None)

    def ft_dt_con(model, i, xyz):
        """friction and virtual sliding alignment"""
        d_t = model.displs[i * 3 + 1] + model.displs[i * 3 + 2]
        f_t = model.forces[i * 3 + 1] + model.forces[i * 3 + 2]
        return (f_t[xyz], -d_t[xyz] * model.alpha[i])

    def penalty_ft_dt_con(model, i, xyz):
        """penalty formulation friction and virtual sliding alignment"""
        d_t = model.displs[i * 3 + 1] + model.displs[i * 3 + 2]
        f_t = model.forces[i * 4 + 2] + model.forces[i * 4 + 3]
        return (f_t[xyz], -d_t[xyz] * model.alpha[i])

    if name == "contact":
        return contact_con
    if name == "penalty_contact":
        return penalty_contact_con
    if name == "fn_np":
        return fn_np_con
    if name == "no_penetration":
        return no_penetration_con
    if name == "ft_dt":
        return ft_dt_con
    if name == "penalty_ft_dt":
        return penalty_ft_dt_con


def static_equilibrium_constraints(model, aeq, afr, p) -> Callable:
    """Create equilibrium and friction constraints.

    Parameters
    ----------
    model : model, optional
        Pyomo model object
    aeq : :class:`~scipy.sparse.csr_matrix`
        Aeq matrix for equilibrium equation.
    afr : :class:`~scipy.sparse.csr_matrix`
        Afr matrix for friction equation,
    p : :class:`~numpy.ndarray`
        External force p.

    Returns
    -------
    Callable, Callable
        equilibrium and friction constraint functions for pyomo

    """

    equilibrium_constraints = MatrixConstraint(
        aeq.data, aeq.indices, aeq.indptr, -p.flatten(), -p.flatten(), model.array_f
    )

    friction_constraint = MatrixConstraint(
        afr.data,
        afr.indices,
        afr.indptr,
        [None for i in range(afr.shape[0])],
        zeros(afr.shape[0]),
        model.array_f,
    )
    return equilibrium_constraints, friction_constraint


def pyomo_result_check(result):
    """Check if pyomo optimisation result, raise error if the problem is infeasible."""
    if (
        result.solver.termination_condition is not pyo.TerminationCondition.optimal
        and result.solver.termination_condition is not pyo.TerminationCondition.feasible
    ):
        raise ValueError(result.solver.termination_condition)

    print("result: ", result.solver.termination_condition)


def pyomo_result_assembly(model, assembly, penalty=False, verbose=False):
    """Save pyomo optimisation results to assembly."""

    shift = 3
    if penalty:
        shift = 4  # for cra_penalty and rbe shift number is 4

    # save force to assembly
    offset = 0
    for edge in assembly.graph.edges():
        interfaces = assembly.graph.edge_attribute(edge, "interfaces")
        for interface in interfaces:
            interface.forces = []
            n = len(interface.points)
            for i in range(n):
                interface.forces.append(
                    {
                        "c_np": model.f[offset + shift * i + 0].value,
                        "c_nn": model.f[offset + shift * i + 1].value if penalty else 0,
                        "c_u": model.f[
                            offset + shift * i + 1 + (1 if penalty else 0)
                        ].value,
                        "c_v": model.f[
                            offset + shift * i + 2 + (1 if penalty else 0)
                        ].value,
                    }
                )
            offset += shift * n

    # save displacement to assembly
    if model.find_component("q") is not None:
        q = [model.q[i].value for i in range(len(model.q))]
        if verbose:
            print("q:", q)
        offset = 0
        for node in assembly.graph.nodes():
            if assembly.graph.node_attribute(node, "is_support"):
                continue
            displacement = q[offset : offset + 6]
            assembly.graph.node_attribute(node, "displacement", displacement)
            offset += 6
