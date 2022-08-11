#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Nonlinear formulation to solve Coupled Rigid-block Equilibrium
Using Pyomo + IPOPT
"""

import numpy as np
import pyomo.environ as pyo
import time

from pyomo.core.base.matrix_constraint import MatrixConstraint
from compas_assembly.datastructures import Assembly
from compas_cra.equilibrium.cra_helper import make_aeq, make_afr, unit_basis
from compas_cra.equilibrium.pyomo_helper import bounds, objectives, constraints
from compas_cra.equilibrium.pyomo_helper import pyomo_result_assembly

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['cra_solve']


def cra_solve(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.,
    d_bnd: float = 1e-3,
    eps: float = 1e-4,
    verbose: bool = False,
    timer: bool = False
):
    """CRA solver using Pyomo + IPOPT. """

    n = assembly.graph.number_of_nodes()
    key_index = {key: index for index, key in enumerate(assembly.graph.nodes())}

    fixed = [key for key in assembly.graph.nodes_where({'is_support': True})]
    fixed = [key_index[key] for key in fixed]
    free = list(set(range(n)) - set(fixed))

    aeqcsr, vcount = make_aeq(assembly)
    aeqcsr = aeqcsr[[index * 6 + i for index in free for i in range(6)], :]
    aeq = aeqcsr.toarray()
    if verbose:
        print("Aeq: ", aeq.shape)

    p = [[0, 0, 0, 0, 0, 0] for i in range(n)]
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, 'block')
        index = key_index[node]
        p[index][2] = -block.volume() * density

    p = np.array(p, dtype=float)
    p = p[free, :].reshape((-1, 1), order='C')

    afrcsr = make_afr(vcount, fcon_number=8, mu=mu)
    afr = afrcsr.toarray()
    if verbose:
        print("Afr: ", afr.shape)

    basis = unit_basis(assembly)
    if verbose:
        print("Unit basis: ", basis.shape)

    model = pyo.ConcreteModel()
    if timer:
        start_time = time.time()

    v_num = vcount  # number of vertices
    free_num = len(free)  # number of free blocks
    v_index = list(range(v_num))  # vertex indices
    f_index = list(range(v_num * 3))  # force indices
    d_index = f_index  # displacement indices
    q_index = list(range(free_num * 6))  # q indices

    bound_f = bounds('f')

    model.fid = pyo.Set(initialize=f_index)
    model.f = pyo.Var(model.fid, initialize=1, domain=bound_f)
    model.q = pyo.Var(q_index, initialize=0)
    model.alpha = pyo.Var(v_index, initialize=0, within=pyo.NonNegativeReals)

    f = np.array([model.f[i] for i in f_index])
    q = np.array([model.q[i] for i in q_index])
    d = aeq.T @ q

    model.d = d
    model.forces = basis * f[:, np.newaxis]  # force x in global coordinate
    model.displs = basis * d[:, np.newaxis]  # displacement d in global coordinate

    obj_cra = objectives('cra')
    bound_d = bounds('d', d_bnd)
    constraint_contact = constraints('contact', eps)
    constraint_no_penetration = constraints('no_penetration', eps)
    constraint_ft_dt = constraints('ft_dt')

    model.obj = pyo.Objective(rule=obj_cra, sense=pyo.minimize)
    model.ceq = MatrixConstraint(aeqcsr.data, aeqcsr.indices, aeqcsr.indptr,
                                 -p.flatten(), -p.flatten(), f)
    model.cfr = MatrixConstraint(afrcsr.data, afrcsr.indices, afrcsr.indptr,
                                 [None for i in range(afr.shape[0])],
                                 np.zeros(afr.shape[0]), f)
    model.d_bnd = pyo.Constraint(d_index, rule=bound_d)
    model.c_con = pyo.Constraint(v_index, rule=constraint_contact)
    model.p_con = pyo.Constraint(v_index, rule=constraint_no_penetration)
    model.ft_dt = pyo.Constraint(v_index, [i for i in range(3)], rule=constraint_ft_dt)

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-8  # same as default tolerance
    solver.options['constr_viol_tol'] = 1e-7  # constraint tolerance
    # https://coin-or.github.io/Ipopt/OPTIONS.html

    results = solver.solve(model, tee=verbose)

    if timer:
        print("--- solving time: %s seconds ---" % (time.time() - start_time))

    if verbose:
        model.f.display()
        model.q.display()
        model.alpha.display()

    if results.solver.termination_condition is not \
       pyo.TerminationCondition.optimal and \
       results.solver.termination_condition is not \
       pyo.TerminationCondition.feasible:
        raise ValueError(results.solver.termination_condition)

    print("result: ", results.solver.termination_condition)

    pyomo_result_assembly(model, assembly, penalty=False, verbose=verbose)

    return assembly


if __name__ == '__main__':
    pass
