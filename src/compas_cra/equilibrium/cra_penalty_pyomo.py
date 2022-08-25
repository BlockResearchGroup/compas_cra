"""Nonlinear formulation with penalty to solve Coupled Rigid-block Equilibrium Using Pyomo + IPOPT"""

import time
import numpy as np
import pyomo.environ as pyo

from compas_assembly.datastructures import Assembly
from .cra_helper import num_vertices, num_free
from .cra_helper import unit_basis
from .cra_helper import equilibrium_setup, friction_setup, external_force_setup
from .pyomo_helper import bounds, objectives, constraints
from .pyomo_helper import static_equilibrium_constraints
from .pyomo_helper import pyomo_result_check, pyomo_result_assembly


def cra_penalty_solve(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.0,
    d_bnd: float = 1e-3,
    eps: float = 1e-4,
    verbose: bool = False,
    timer: bool = False,
) -> Assembly:
    r"""CRA solver with penalty formulation using Pyomo + IPOPT.

    Parameters
    ----------
    assembly : :class:`~compas_assembly.datastructures.Assembly`
        The rigid block assembly.
    mu : float, optional
        Friction coefficient value.
    density : float, optional
        Density of the block material.
    d_bnd : float, optional
        The bound of virtual displacement d.
    eps : float, optional
        Epsilon, contact overlapping parameter.
    verbose : bool, optional
        Print information during the execution of the algorithm.
    timer : bool, optional
        Time the solving time.

    Returns
    -------
    :class:`~compas_assembly.datastructures.Assembly`
        The assembly is updated in place, also return Assembly for compas.rpc and compas.cloud


    Notes
    -----
     This function solves the following optimisation problem, `Eq.(14) <https://www.sciencedirect.com/science/article/pii/S0010448522000161?via%3Dihub#fd14>`_ :

    .. math::

        \begin{align}
            \begin{split}
                \min_{{\bf{\tilde{f}}},\, \delta{\bf{q}},\, \pmb\alpha} \quad & \left\| {\bf{f}}_n^{+} \right\|_2^2 +
                \gamma \left\| {\bf{f}}_n^{-} \right\|_2^2 + \left\| \pmb\alpha \right\|_2^2 \\
                \textrm{s.t.} \quad & {{\bf{A}}_{eq}}\:{\bf{B}}\:{\bf{\tilde{f}}} = -{\bf{p}} \\
                & {\bf{A}}_{fr}\:{\bf{B}}\:{\bf{\tilde{f}}} \le {\bf{0}} \\
                & {\bf{A}}_{eq}^\intercal\:{\delta\bf{q}} = \delta{\bf{d}} \\
                & {f_{jkn}^{i+}}\:({\delta d_{jkn}^i} + \varepsilon) = 0 \\
                & {f_{jkn}^{i+}} \: {f_{jkn}^{i-}} = 0 \\
                & {\bf{f}}_{jkt}^{i} = -{\alpha_{jk}^i} \: \delta{\bf{d}}_{jkt}^{i} \\
                & \left\lvert\, \delta {\bf{d}}_{jk\cdot}^{i} \,\right\lvert \le \eta \\
                & \gamma \, , \varepsilon \, , \eta  \ge 0 \\
                & f_{jkn}^{i+}\, ,f_{jkn}^{i-} \, , {\alpha_{jk}^i} \, , ({\delta d_{jkn}^i} + \varepsilon) \ge 0 \;, \quad \forall i,j,k \;,
            \end{split}
        \end{align}

    For more information please check our research paper:
    `Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies <https://doi.org/10.1016/j.cad.2022.103216>`_

    """

    if timer:
        start_time = time.time()

    model = pyo.ConcreteModel()

    v_num = num_vertices(assembly)  # number of vertices
    free_num = num_free(assembly)  # number of free blocks
    f_basis = unit_basis(assembly, True)
    d_basis = unit_basis(assembly)

    model.v_id = pyo.Set(initialize=range(v_num))  # vertex indices
    model.f_id = pyo.Set(initialize=range(v_num * 4))  # force indices
    model.d_id = pyo.Set(initialize=range(v_num * 3))  # displacement indices
    model.q_id = pyo.Set(initialize=range(free_num * 6))  # q indices

    model.f = pyo.Var(model.f_id, initialize=0, domain=bounds("f_tilde"))
    model.q = pyo.Var(model.q_id, initialize=0)
    model.alpha = pyo.Var(model.v_id, initialize=0, within=pyo.NonNegativeReals)

    model.array_f = np.array([model.f[i] for i in model.f_id])
    model.array_q = np.array([model.q[i] for i in model.q_id])

    aeq = equilibrium_setup(assembly, penalty=False)
    aeq_b = equilibrium_setup(assembly, penalty=True)
    afr_b = friction_setup(assembly, mu, penalty=True)
    p = external_force_setup(assembly, density)

    model.d = aeq.toarray().T @ model.array_q
    model.forces = (
        f_basis * model.array_f[:, np.newaxis]
    )  # force x in global coordinate
    model.displs = (
        d_basis * model.d[:, np.newaxis]
    )  # displacement d in global coordinate

    obj_cra_penalty = objectives("cra_penalty")
    bound_d = bounds("d", d_bnd)
    constraint_contact = constraints("penalty_contact", eps)
    constraint_no_penetration = constraints("no_penetration", eps)
    constraint_penalty_ft_dt = constraints("penalty_ft_dt")
    constraint_fn_np = constraints("fn_np")

    eq_con, fr_con = static_equilibrium_constraints(model, aeq_b, afr_b, p)

    model.obj = pyo.Objective(rule=obj_cra_penalty, sense=pyo.minimize)
    model.ceq = eq_con
    model.cfr = fr_con
    model.d_bnd = pyo.Constraint(model.d_id, rule=bound_d)
    model.c_con = pyo.Constraint(model.v_id, rule=constraint_contact)
    model.p_con = pyo.Constraint(model.v_id, rule=constraint_no_penetration)
    model.fn_np = pyo.Constraint(model.v_id, rule=constraint_fn_np)
    model.ft_dt = pyo.Constraint(
        model.v_id, [i for i in range(3)], rule=constraint_penalty_ft_dt
    )

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory("ipopt")
    solver.options["tol"] = 1e-8  # same as default tolerance
    solver.options["constr_viol_tol"] = 1e-7  # constraint tolerance
    solver.options["acceptable_tol"] = 1e-6
    solver.options["acceptable_constr_viol_tol"] = 1e-5
    # https://coin-or.github.io/Ipopt/OPTIONS.html
    result = solver.solve(model, tee=verbose)

    if timer:
        print("--- solving time: %s seconds ---" % (time.time() - start_time))

    if verbose:
        model.f.display()
        model.q.display()
        model.alpha.display()
        print("objective value: ")
        model.obj.display()

    pyomo_result_check(result)
    pyomo_result_assembly(model, assembly, penalty=True, verbose=verbose)

    return assembly
