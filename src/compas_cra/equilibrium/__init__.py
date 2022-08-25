"""
********************************************************************************
equilibrium
********************************************************************************

.. currentmodule:: compas_cra.equilibrium


Solvers
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    cra_solve
    cra_penalty_solve
    rbe_solve

--------------------------

The following helper functions can be useful if you're developing your own formulation.

Equilibrium Helper Functions
============================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    equilibrium_setup
    friction_setup
    external_force_setup
    density_setup
    make_aeq
    make_afr
    unit_basis
    num_vertices
    num_free
    free_nodes

Pyomo Helper Functions
======================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    initialisations
    bounds
    objectives
    constraints
    static_equilibrium_constraints
    pyomo_result_check
    pyomo_result_assembly

"""

from .cra_pyomo import cra_solve
from .cra_penalty_pyomo import cra_penalty_solve
from .rbe_pyomo import rbe_solve
from .cra_helper import (
    equilibrium_setup,
    friction_setup,
    external_force_setup,
    density_setup,
    make_aeq,
    make_afr,
    unit_basis,
    num_vertices,
    num_free,
    free_nodes,
)
from .pyomo_helper import (
    initialisations,
    bounds,
    objectives,
    constraints,
    static_equilibrium_constraints,
    pyomo_result_check,
    pyomo_result_assembly,
)

__all__ = [
    "cra_solve",
    "cra_penalty_solve",
    "rbe_solve",
    "equilibrium_setup",
    "friction_setup",
    "external_force_setup",
    "density_setup",
    "make_aeq",
    "make_afr",
    "unit_basis",
    "num_vertices",
    "num_free",
    "free_nodes",
    "initialisations",
    "bounds",
    "objectives",
    "constraints",
    "static_equilibrium_constraints",
    "pyomo_result_check",
    "pyomo_result_assembly",
]
