"""Solver-agnostic sparse NLP layer.

:class:`~compas_cra.nlp.problem.NLPProblem` describes a sparse nonlinear program;
:func:`~compas_cra.nlp.solve_nlp` solves one with the best available backend
(currently the in-process IPOPT binding in :mod:`compas_cra._native`).
"""

from .backends import available_backends
from .backends import solve_nlp
from .problem import NLPProblem
from .problem import NLPResult

__all__ = ["NLPProblem", "NLPResult", "solve_nlp", "available_backends"]
