#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Nonlinear formulation with penalty to solve Coupled Rigid-block Equilibrium
Using Pyomo + IPOPT
"""

import numpy as np
import pyomo.environ as pyo
import time

from pyomo.core.base.matrix_constraint import MatrixConstraint
from compas_assembly.datastructures import Assembly
from compas_cra.equilibrium.cra_helper import make_aeq, unit_basis
from compas_cra.equilibrium.cra_penalty_helper import make_aeq_b, make_afr_b, unit_basis_penalty
from compas_cra.equilibrium.pyomo_helper import bounds, objectives


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

    n = assembly.graph.number_of_nodes()
    key_index = {key: index for index, key in enumerate(assembly.graph.nodes())}

    fixed = [key for key in assembly.graph.nodes_where({'is_support': True})]
    fixed = [key_index[key] for key in fixed]
    free = list(set(range(n)) - set(fixed))

    aeq_csr, vcount = make_aeq(assembly)
    aeq_csr = aeq_csr[[index * 6 + i for index in free for i in range(6)], :]
    aeq = aeq_csr.toarray()

    aeq_b_csr, _ = make_aeq_b(assembly)
    aeq_b_csr = aeq_b_csr[[index * 6 + i for index in free for i in range(6)], :]

    p = [[0, 0, 0, 0, 0, 0] for i in range(n)]
    for node in assembly.graph.nodes():
        block = assembly.graph.node_attribute(node, 'block')
        index = key_index[node]
        p[index][2] = -1 * block.volume() * density

    p = np.array(p, dtype=float)
    p = p[free, :].reshape((-1, 1), order='C')

    afr_b_csr = make_afr_b(vcount, mu=mu, friction_net=False)
    afr_b = afr_b_csr.toarray()

    f_basis = unit_basis_penalty(assembly)
    d_basis = unit_basis(assembly)

    model = pyo.ConcreteModel()
    if timer:
        start_time = time.time()

    v_num = vcount  # number of vertices
    free_num = len(free)  # number of free blocks
    v_index = list(range(v_num))  # vertex indices
    f_index = list(range(v_num * 4))  # force indices
    d_index = list(range(v_num * 3))  # displacement indices
    q_index = list(range(free_num * 6))  # q indices

    bound_f_tilde = bounds('f_tilde')

    model.fid = pyo.Set(initialize=f_index)
    # model.f = pyo.Var(f_index, initialize=f_init, domain=f_bnds)
    model.f = pyo.Var(model.fid, initialize=0, domain=bound_f_tilde)
    model.q = pyo.Var(q_index, initialize=0)
    model.alpha = pyo.Var(v_index, initialize=0, within=pyo.NonNegativeReals)

    f = np.array([model.f[i] for i in f_index])
    q = np.array([model.q[i] for i in q_index])
    d = aeq.T @ q

    forces = f_basis * f[:, np.newaxis]  # force x in global coordinate
    displs = d_basis * d[:, np.newaxis]  # displacement d in global coordinate

    model.d = d

    def contact_con(m, t):
        dn = m.d[t * 3]
        fn = m.f[t * 4]
        return ((dn + eps) * fn, 0)

    def fnc_con(m, t):  # fn+ and fn- cannot coexist
        return (m.f[t * 4] * m.f[t * 4 + 1], 0)

    def nonpen_con(m, t):
        return (0, m.d[t * 3] + eps, None)

    def ftdt_con(m, t, xyz):
        dt = displs[t * 3 + 1] + displs[t * 3 + 2]
        ft = forces[t * 4 + 2] + forces[t * 4 + 3]
        return (ft[xyz], -dt[xyz] * m.alpha[t])

    obj_cra_penalty = objectives(solver='cra_penalty')
    bound_d = bounds(variable='d', d_bnd=d_bnd)

    model.obj = pyo.Objective(rule=obj_cra_penalty, sense=pyo.minimize)
    model.ceq = MatrixConstraint(aeq_b_csr.data, aeq_b_csr.indices, aeq_b_csr.indptr,
                                 -p.flatten(), -p.flatten(), f)
    model.cfr = MatrixConstraint(afr_b_csr.data, afr_b_csr.indices, afr_b_csr.indptr,
                                 [None for i in range(afr_b.shape[0])],
                                 np.zeros(afr_b.shape[0]), f)
    model.dbnd = pyo.Constraint(d_index, rule=bound_d)
    model.ccon = pyo.Constraint(v_index, rule=contact_con)
    model.pcon = pyo.Constraint(v_index, rule=nonpen_con)

    model.fncon = pyo.Constraint(v_index, rule=fnc_con)

    model.cftdt = pyo.Constraint(v_index, [i for i in range(3)], rule=ftdt_con)

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-8
    solver.options['constr_viol_tol'] = 1e-7
    results = solver.solve(model, tee=verbose)

    if timer:
        print("--- solving time: %s seconds ---" % (time.time() - start_time))

    # f_star = np.array([model.f[i].value for i in f_index])
    # q_star = np.array([model.q[i].value for i in q_index])
    # d_star = aeq.T @ q_star
    # for v in v_index:
    #     fnp = f_star[v * 4]
    #     fnn = f_star[v * 4 + 1]
    #     dn = d_star[v * 3]
    #     print("fn: ", fnp, ", dn: ", dn, ", fn * dn: ", (fnp - fnn) * dn)

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
    print("obj", model.obj.display())

    # =========================================================================
    # assign forces and object displacement y
    offset = 0
    for edge in assembly.graph.edges():
        interfaces = assembly.graph.edge_attribute(edge, 'interfaces')
        for interface in interfaces:
            interface.forces = []
            n = len(interface.points)
            for i in range(n):
                interface.forces.append({
                    'c_np': model.f[offset + 4 * i + 0].value,
                    'c_nn': model.f[offset + 4 * i + 1].value,
                    'c_u': model.f[offset + 4 * i + 2].value,
                    'c_v': model.f[offset + 4 * i + 3].value
                })

            offset += 4 * n

    q = [model.q[i].value * 1 for i in range(6 * free_num)]
    if verbose:
        print("q:", q)

    offset = 0
    for node in assembly.graph.nodes():
        if assembly.graph.node_attribute(node, 'is_support'):
            continue
        displacement = q[offset:offset+6]
        assembly.graph.node_attribute(node, 'displacement', displacement)
        offset += 6
    # =========================================================================
