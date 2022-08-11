#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Rigid-block Equilibrium
Using Pyomo + MOSEK
"""

import numpy as np
import pyomo.environ as pyo
import time

from pyomo.core.base.matrix_constraint import MatrixConstraint
from compas_assembly.datastructures import Assembly
from compas_cra.equilibrium.cra_penalty_helper import make_aeq_b, make_afr_b
from compas_cra.equilibrium.pyomo_helper import bounds, objectives
from compas_cra.equilibrium.pyomo_helper import pyomo_result_assembly

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['rbe_solve']


def rbe_solve(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.,
    verbose: bool = False,
    timer: bool = False
):
    """RBE solver with penalty formulation using Pyomo + MOSEK. """

    n = assembly.graph.number_of_nodes()
    key_index = {key: index for index, key in enumerate(assembly.graph.nodes())}

    fixed = [key for key in assembly.graph.nodes_where({'is_support': True})]
    fixed = [key_index[key] for key in fixed]
    free = list(set(range(n)) - set(fixed))

    aeq_b_csr, vcount = make_aeq_b(assembly)
    aeq_b_csr = aeq_b_csr[[index * 6 + i for index in free for i in range(6)], :]

    p = [[0, 0, 0, 0, 0, 0] for i in range(n)]
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, 'block')
        index = key_index[node]
        p[index][2] = -1 * block.volume() * density

    p = np.array(p, dtype=float)
    p = p[free, :].reshape((-1, 1), order='C')

    afr_b_csr = make_afr_b(vcount, mu=mu, friction_net=False)
    # afr_b = afr_b_csr.toarray()

    model = pyo.ConcreteModel()
    if timer:
        start_time = time.time()

    v_num = vcount  # number of vertices
    f_index = list(range(v_num * 4))  # force indices

    model.f_id = pyo.Set(initialize=f_index)
    # free_num = len(free)  # number
    # eq_index = [i for i in range(6 * free_num)]
    # fr_index = [i for i in range(v_num * 8)]  # friction constraint indices

    bound_f_tilde = bounds('f_tilde')
    # model.f = pyo.Var(f_index, initialize=0, domain=f_tilde_bnds)
    model.f = pyo.Var(model.f_id, initialize=0, domain=bound_f_tilde)
    # model.f = pyo.Var(f_index, initialize=f_init, domain=f_tilde_bnds)

    f = np.array([model.f[i] for i in f_index])

    # def eq_con(m, t):
    #     return (sum(aeq_b_csr[t, i] * m.f[i] for i in f_index), -p[t][0])
    #
    # def fr_con(m, t):
    #     return (None, sum(afr_b[t, i] * m.f[i] for i in f_index), 0)

    obj_rbe = objectives('rbe')

    model.obj = pyo.Objective(rule=obj_rbe, sense=pyo.minimize)
    model.ceq = MatrixConstraint(aeq_b_csr.data, aeq_b_csr.indices, aeq_b_csr.indptr,
                                 -p.flatten(), -p.flatten(), f)
    model.cfr = MatrixConstraint(afr_b_csr.data, afr_b_csr.indices, afr_b_csr.indptr,
                                 [None for i in range(afr_b_csr.shape[0])],
                                 np.zeros(afr_b_csr.shape[0]), f)
    # model.ceq = pyo.Constraint(eq_index, rule=eq_con)
    # model.cfr = pyo.Constraint(fr_index, rule=fr_con)

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory('ipopt')
    results = solver.solve(model, tee=verbose)

    if timer:
        print("--- solving time: %s seconds ---" % (time.time() - start_time))

    if verbose:
        model.f.display()

    if results.solver.termination_condition is not \
       pyo.TerminationCondition.optimal and \
       results.solver.termination_condition is not \
       pyo.TerminationCondition.feasible:
        raise ValueError(results.solver.termination_condition)

    print("result: ", results.solver.termination_condition)
    print("obj", model.obj.display())

    model.f.display()

    pyomo_result_assembly(model, assembly, penalty=True, verbose=verbose)

    return assembly
