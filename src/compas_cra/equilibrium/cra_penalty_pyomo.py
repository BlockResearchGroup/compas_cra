#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Nonlinear formulation with penalty to solve Coupled Rigid-block Equilibrium
Using Pyomo + IPOPT
"""

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


__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['cra_penalty_solve']


def cra_penalty_solve(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.,
    d_bnd: float = 1e-3,
    eps: float = 1e-4,
    verbose: bool = False,
    timer: bool = False
):
    """CRA solver with penalty formulation using Pyomo + IPOPT. """
    # TODO: docs

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

    model.f = pyo.Var(model.f_id, initialize=0, domain=bounds('f_tilde'))
    model.q = pyo.Var(model.q_id, initialize=0)
    model.alpha = pyo.Var(model.v_id, initialize=0, within=pyo.NonNegativeReals)

    model.array_f = np.array([model.f[i] for i in model.f_id])
    model.array_q = np.array([model.q[i] for i in model.q_id])

    aeq = equilibrium_setup(assembly, penalty=False)
    aeq_b = equilibrium_setup(assembly, penalty=True)
    afr_b = friction_setup(assembly, mu, penalty=True)
    p = external_force_setup(assembly, density)

    model.d = aeq.toarray().T @ model.array_q
    model.forces = f_basis * model.array_f[:, np.newaxis]  # force x in global coordinate
    model.displs = d_basis * model.d[:, np.newaxis]  # displacement d in global coordinate

    obj_cra_penalty = objectives('cra_penalty')
    bound_d = bounds('d', d_bnd)
    constraint_contact = constraints('penalty_contact', eps)
    constraint_no_penetration = constraints('no_penetration', eps)
    constraint_penalty_ft_dt = constraints('penalty_ft_dt')
    constraint_fn_np = constraints('fn_np')

    eq_con, fr_con = static_equilibrium_constraints(model, aeq_b, afr_b, p)

    model.obj = pyo.Objective(rule=obj_cra_penalty, sense=pyo.minimize)
    model.ceq = eq_con
    model.cfr = fr_con
    model.d_bnd = pyo.Constraint(model.d_id, rule=bound_d)
    model.c_con = pyo.Constraint(model.v_id, rule=constraint_contact)
    model.p_con = pyo.Constraint(model.v_id, rule=constraint_no_penetration)
    model.fn_np = pyo.Constraint(model.v_id, rule=constraint_fn_np)
    model.ft_dt = pyo.Constraint(model.v_id, [i for i in range(3)], rule=constraint_penalty_ft_dt)

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-8  # same as default tolerance
    solver.options['constr_viol_tol'] = 1e-7  # constraint tolerance
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
