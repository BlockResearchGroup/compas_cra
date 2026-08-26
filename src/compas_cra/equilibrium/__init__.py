from .cra_native import cra_solve_native
from .cra_native import cra_penalty_solve_native
from .cra_native import rbe_solve_native

# The solvers run on the in-process IPOPT binding (compas_cra._native); the
# historical names are the canonical API, and are aliases of the same functions.
from .cra_native import cra_solve_native as cra_solve
from .cra_native import cra_penalty_solve_native as cra_penalty_solve
from .cra_native import rbe_solve_native as rbe_solve
from .cra_nlp import cra_problem
from .cra_nlp import cra_penalty_problem
from .cra_nlp import rbe_problem
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

__all__ = [
    "cra_solve",
    "cra_penalty_solve",
    "rbe_solve",
    "cra_solve_native",
    "cra_penalty_solve_native",
    "rbe_solve_native",
    "cra_problem",
    "cra_penalty_problem",
    "rbe_problem",
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
]
