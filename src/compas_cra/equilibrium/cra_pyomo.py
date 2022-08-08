#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Nonlinear formulation to solve Coupled Rigid-block Equilibrium
Using Pyomo + IPOPT
"""

import numpy as np
import pyomo.environ as pyo
import time

from compas_cra.equilibrium.cra_helper import make_aeq, make_afr, unit_basis
from compas_cra.equilibrium.pyomo_helper import f_bnds
from compas_cra.equilibrium.pyomo_helper import obj_cra
from pyomo.core.base.matrix_constraint import MatrixConstraint
from compas_assembly.datastructures import Assembly

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
    #
    # assembly.graph.node_attribute(0, "is_support", True)

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
    v_index = [i for i in range(v_num)]  # vertex indices
    f_index = [i for i in range(v_num * 3)]  # force indices
    d_index = f_index  # displacement indices
    q_index = [i for i in range(free_num * 6)]  # q indices

    model.fid = pyo.Set(initialize=f_index)
    model.f = pyo.Var(model.fid, initialize=1, domain=f_bnds)
    model.q = pyo.Var(q_index, initialize=0)
    model.alpha = pyo.Var(v_index, initialize=0, within=pyo.NonNegativeReals)

    f = np.array([model.f[i] for i in f_index])
    q = np.array([model.q[i] for i in q_index])
    d = aeq.T @ q

    forces = basis * f[:, np.newaxis]  # force x in global coordinate
    displs = basis * d[:, np.newaxis]  # displacement d in global coordinate

    model.d = d

    def d_bnds(m, t):
        return (-d_bnd, d[t], d_bnd)

    def contact_con(m, t):
        dn = m.d[t * 3]
        fn = m.f[t * 3]
        return ((dn + eps) * fn, 0)

    def nonpen_con(m, t):
        return (0, m.d[t * 3] + eps, None)

    def ftdt_con(m, t, xyz):
        dt = displs[t * 3 + 1] + displs[t * 3 + 2]
        ft = forces[t * 3 + 1] + forces[t * 3 + 2]
        return (ft[xyz], -dt[xyz] * m.alpha[t])

    model.obj = pyo.Objective(rule=obj_cra, sense=pyo.minimize)
    model.ceq = MatrixConstraint(aeqcsr.data, aeqcsr.indices, aeqcsr.indptr,
                                 -p.flatten(), -p.flatten(), f)
    model.cfr = MatrixConstraint(afrcsr.data, afrcsr.indices, afrcsr.indptr,
                                 [None for i in range(afr.shape[0])],
                                 np.zeros(afr.shape[0]), f)
    model.dbnd = pyo.Constraint(d_index, rule=d_bnds)
    model.ccon = pyo.Constraint(v_index, rule=contact_con)
    model.pcon = pyo.Constraint(v_index, rule=nonpen_con)
    model.cftdt = pyo.Constraint(v_index, [i for i in range(3)], rule=ftdt_con)

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

    # =========================================================================
    # save forces to assembly
    offset = 0
    for edge in assembly.graph.edges():
        interfaces = assembly.graph.edge_attribute(edge, 'interfaces')
        for interface in interfaces:
            interface.forces = []
            n = len(interface.points)
            for i in range(n):
                interface.forces.append({
                    'c_np': model.f[offset + 3 * i + 0].value,
                    'c_nn': 0,
                    'c_u': model.f[offset + 3 * i + 1].value,
                    'c_v': model.f[offset + 3 * i + 2].value
                })
            offset += 3 * n

        # interface = assembly.graph.edge_attribute(edge, 'interface')
        # interface.forces = []
        # n = len(interface.points)
        # for i in range(n):
        #     interface.forces.append({
        #         'c_np': model.f[offset + 3 * i + 0].value,
        #         'c_nn': 0,
        #         'c_u': model.f[offset + 3 * i + 1].value,
        #         'c_v': model.f[offset + 3 * i + 2].value
        #     })
        # offset += 3 * n

    # save displacements to assembly
    q = [model.q[i].value * 1 for i in range(6 * free_num)]
    # d = aeq.T @ q
    # f = [model.f[i].value for i in f_index]
    # for v in v_index:
    #     dn = d[v * 3]
    #     fn = f[v * 3]
    #     print("contact_con ", (dn + eps) * fn)

    if verbose:
        print("q:", q)
    offset = 0
    for node in assembly.graph.nodes():
        if assembly.graph.node_attribute(node, 'is_support'):
            continue
        displacement = q[offset:offset + 6]
        assembly.graph.node_attribute(node, 'displacement', displacement)
        offset += 6

    return assembly


if __name__ == '__main__':
    pass
