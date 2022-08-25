"""Rigid-block Equilibrium Using Pyomo + IPOPT"""

import time
import numpy as np
import pyomo.environ as pyo

from compas_assembly.datastructures import Assembly
from .cra_helper import num_vertices
from .cra_helper import equilibrium_setup, friction_setup, external_force_setup
from .pyomo_helper import bounds, objectives
from .pyomo_helper import static_equilibrium_constraints
from .pyomo_helper import pyomo_result_check, pyomo_result_assembly


def rbe_solve(
    assembly: Assembly,
    mu: float = 0.84,
    density: float = 1.0,
    verbose: bool = False,
    timer: bool = False,
) -> Assembly:
    r"""RBE solver with penalty formulation using Pyomo + IPOPT.

    Parameters
    ----------
    assembly : :class:`~compas_assembly.datastructures.Assembly`
        The rigid block assembly.
    mu : float, optional
        Friction coefficient value.
    density : float, optional
        Density of the block material.
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
    This function solves the following optimisation problem, `Eq.(6) <https://www.sciencedirect.com/science/article/pii/S0010448522000161?via%3Dihub#fd6>`_ :

    .. math::

        \begin{align}
            \begin{split}
                \min_{\bf{\tilde{f}}} \quad & \frac{1}{2}\:{\bf{\tilde{f}}^\intercal}\:{\bf{H}}\:{\bf{\tilde{f}}} \\
                \textrm{s.t.} \quad & {{\bf{A}}_{eq}}\:{\bf{B}}\:{\bf{\tilde{f}}} = -{\bf{p}} \\
                & {\bf{A}}_{fr}\:{\bf{B}}\:{\bf{\tilde{f}}} \le {\bf{0}} \\
                & f_{jkn}^{i+}\, ,f_{jkn}^{i-} \ge 0 \;, \quad \forall i,j,k \;,
            \end{split}
        \end{align}

    For more information please check our research paper:
    `Coupled Rigid-Block Analysis: Stability-Aware Design of Complex Discrete-Element Assemblies <https://doi.org/10.1016/j.cad.2022.103216>`_

    """

    model = pyo.ConcreteModel()

    if timer:
        start_time = time.time()

    v_num = num_vertices(assembly)  # number of vertices

    model.f_id = pyo.Set(initialize=range(v_num * 4))  # force indices
    model.f = pyo.Var(model.f_id, initialize=0, domain=bounds("f_tilde"))
    model.array_f = np.array([model.f[i] for i in model.f_id])

    aeq_b = equilibrium_setup(assembly, penalty=True)
    afr_b = friction_setup(assembly, mu, penalty=True)
    p = external_force_setup(assembly, density)

    obj_rbe = objectives("rbe", (0, 1e0, 1e6, 1e0))
    eq_con, fr_con = static_equilibrium_constraints(model, aeq_b, afr_b, p)

    model.obj = pyo.Objective(rule=obj_rbe, sense=pyo.minimize)
    model.ceq = eq_con
    model.cfr = fr_con

    if timer:
        print("--- set up time: %s seconds ---" % (time.time() - start_time))
    print("finished setup... now trying to solve it...")
    if timer:
        start_time = time.time()

    solver = pyo.SolverFactory("ipopt")
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
