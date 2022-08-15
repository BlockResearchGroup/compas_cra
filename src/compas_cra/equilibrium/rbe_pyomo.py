#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Rigid-block Equilibrium
Using Pyomo + IPOPT
"""

import time
import numpy as np
import pyomo.environ as pyo

from compas_assembly.datastructures import Assembly
from .cra_helper import num_vertices
from .cra_helper import equilibrium_setup, friction_setup, external_force_setup
from .pyomo_helper import bounds, objectives
from .pyomo_helper import static_equilibrium_constraints
from .pyomo_helper import pyomo_result_check, pyomo_result_assembly

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
    """RBE solver with penalty formulation using Pyomo + IPOPT. """

    model = pyo.ConcreteModel()

    if timer:
        start_time = time.time()

    v_num = num_vertices(assembly)  # number of vertices

    model.f_id = pyo.Set(initialize=range(v_num * 4))  # force indices
    model.f = pyo.Var(model.f_id, initialize=0, domain=bounds('f_tilde'))
    model.array_f = np.array([model.f[i] for i in model.f_id])

    aeq_b = equilibrium_setup(assembly, penalty=True)
    afr_b = friction_setup(assembly, mu, penalty=True)
    p = external_force_setup(assembly, density)

    obj_rbe = objectives('rbe', (0, 0, 1e+6))
    eq_con, fr_con = static_equilibrium_constraints(model, aeq_b, afr_b, p)

    model.obj = pyo.Objective(rule=obj_rbe, sense=pyo.minimize)
    model.ceq = eq_con
    model.cfr = fr_con

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory('ipopt')
    result = solver.solve(model, tee=verbose)

    if timer:
        print("--- solving time: %s seconds ---" % (time.time() - start_time))

    if verbose:
        model.f.display()
        print("objective value: ")
        model.obj.display()

    pyomo_result_check(result)
    pyomo_result_assembly(model, assembly, penalty=True, verbose=verbose)

    return assembly
