"""Solver-agnostic description of a sparse nonlinear program.

The problem form is the standard NLP

.. math::

    \\min_x f(x) \\quad \\textrm{s.t.} \\quad g_L \\le g(x) \\le g_U, \\; x_L \\le x \\le x_U

with first derivatives given as a sparse Jacobian in fixed triplet structure, and an
optional exact Hessian of the Lagrangian (lower triangle, fixed triplet structure,
duplicate entries are summed). This mirrors what interior point solvers such as IPOPT
consume, without depending on any particular solver package.
"""

import numpy as np

__all__ = ["NLPProblem", "NLPResult"]


class NLPProblem:
    """A sparse nonlinear program with callback-based evaluation.

    Parameters
    ----------
    n : int
        Number of variables.
    x_l, x_u : ndarray of shape (n,)
        Variable bounds. Use -inf/inf for free variables.
    g_l, g_u : ndarray of shape (m,)
        Constraint bounds. Equal entries make an equality constraint.
    x0 : ndarray of shape (n,)
        Starting point.
    objective : callable
        ``objective(x) -> float``
    gradient : callable
        ``gradient(x) -> ndarray (n,)``
    constraints : callable
        ``constraints(x) -> ndarray (m,)``
    jac_rows, jac_cols : ndarray of int
        Fixed sparsity structure of the constraint Jacobian (0-based).
    jacobian : callable
        ``jacobian(x) -> ndarray (nnz,)`` values matching ``jac_rows``/``jac_cols``.
    hess_rows, hess_cols : ndarray of int, optional
        Fixed lower-triangle sparsity of the Lagrangian Hessian (0-based, row >= col).
        Duplicate entries are summed by the solver.
    hessian : callable, optional
        ``hessian(x, sigma, lam) -> ndarray (nnz_h,)`` values of
        ``sigma * H(f) + sum_i lam[i] * H(g_i)``. When omitted, solvers fall back to a
        quasi-Newton approximation.
    """

    def __init__(
        self,
        n,
        x_l,
        x_u,
        g_l,
        g_u,
        x0,
        objective,
        gradient,
        constraints,
        jac_rows,
        jac_cols,
        jacobian,
        hess_rows=None,
        hess_cols=None,
        hessian=None,
    ):
        self.n = int(n)
        self.x_l = np.asarray(x_l, dtype=float)
        self.x_u = np.asarray(x_u, dtype=float)
        self.g_l = np.asarray(g_l, dtype=float)
        self.g_u = np.asarray(g_u, dtype=float)
        self.m = self.g_l.shape[0]
        self.x0 = np.asarray(x0, dtype=float)
        self.objective = objective
        self.gradient = gradient
        self.constraints = constraints
        self.jac_rows = np.asarray(jac_rows, dtype=np.int32)
        self.jac_cols = np.asarray(jac_cols, dtype=np.int32)
        self.jacobian = jacobian
        self.hess_rows = None if hess_rows is None else np.asarray(hess_rows, dtype=np.int32)
        self.hess_cols = None if hess_cols is None else np.asarray(hess_cols, dtype=np.int32)
        self.hessian = hessian
        self._validate()

    def _validate(self):
        if self.x_l.shape != (self.n,) or self.x_u.shape != (self.n,) or self.x0.shape != (self.n,):
            raise ValueError("variable bounds and starting point must have shape (n,)")
        if self.g_u.shape != (self.m,):
            raise ValueError("constraint bounds must have equal shape (m,)")
        if self.jac_rows.shape != self.jac_cols.shape:
            raise ValueError("jac_rows and jac_cols must have the same shape")
        if (self.hessian is None) != (self.hess_rows is None) or (self.hessian is None) != (self.hess_cols is None):
            raise ValueError("hessian, hess_rows and hess_cols must be given together")

    @property
    def has_hessian(self):
        return self.hessian is not None

    # ------------------------------------------------------------------
    # verification helpers (used by the test suite, handy for debugging)
    # ------------------------------------------------------------------

    def check_gradient(self, x=None, step=1e-6):
        """Max abs difference between the analytic gradient and finite differences."""
        x = self.x0 if x is None else np.asarray(x, dtype=float)
        analytic = np.asarray(self.gradient(x), dtype=float)
        fd = np.empty_like(analytic)
        for i in range(self.n):
            e = np.zeros(self.n)
            e[i] = step
            fd[i] = (self.objective(x + e) - self.objective(x - e)) / (2 * step)
        return float(np.max(np.abs(analytic - fd))) if self.n else 0.0

    def check_jacobian(self, x=None, step=1e-6):
        """Max abs difference between the analytic Jacobian and finite differences."""
        from scipy.sparse import coo_matrix

        x = self.x0 if x is None else np.asarray(x, dtype=float)
        values = np.asarray(self.jacobian(x), dtype=float)
        analytic = coo_matrix((values, (self.jac_rows, self.jac_cols)), shape=(self.m, self.n)).toarray()
        fd = np.empty((self.m, self.n))
        for i in range(self.n):
            e = np.zeros(self.n)
            e[i] = step
            fd[:, i] = (self.constraints(x + e) - self.constraints(x - e)) / (2 * step)
        return float(np.max(np.abs(analytic - fd))) if self.m and self.n else 0.0

    def check_hessian(self, x=None, sigma=1.0, lam=None, step=1e-5):
        """Max abs difference between the analytic Lagrangian Hessian and finite differences."""
        from scipy.sparse import coo_matrix

        if not self.has_hessian:
            raise ValueError("problem has no hessian callback")
        x = self.x0 if x is None else np.asarray(x, dtype=float)
        lam = np.ones(self.m) if lam is None else np.asarray(lam, dtype=float)
        values = np.asarray(self.hessian(x, sigma, lam), dtype=float)
        lower = coo_matrix((values, (self.hess_rows, self.hess_cols)), shape=(self.n, self.n)).toarray()
        analytic = lower + lower.T - np.diag(np.diag(lower))

        def lagrangian_grad(x):
            return sigma * np.asarray(self.gradient(x), dtype=float) + lam @ self._jac_dense(x)

        fd = np.empty((self.n, self.n))
        for i in range(self.n):
            e = np.zeros(self.n)
            e[i] = step
            fd[:, i] = (lagrangian_grad(x + e) - lagrangian_grad(x - e)) / (2 * step)
        fd = (fd + fd.T) / 2
        return float(np.max(np.abs(analytic - fd)))

    def _jac_dense(self, x):
        from scipy.sparse import coo_matrix

        values = np.asarray(self.jacobian(x), dtype=float)
        return coo_matrix((values, (self.jac_rows, self.jac_cols)), shape=(self.m, self.n)).toarray()


class NLPResult:
    """Result of an NLP solve."""

    OPTIMAL = "optimal"
    ACCEPTABLE = "acceptable"
    INFEASIBLE = "infeasible"
    FAILED = "failed"

    def __init__(self, x, obj, status, status_message, g=None, mult_g=None, iterations=None):
        self.x = np.asarray(x, dtype=float)
        self.obj = float(obj)
        self.status = status
        self.status_message = status_message
        self.g = g
        self.mult_g = mult_g
        self.iterations = iterations

    @property
    def success(self):
        return self.status in (self.OPTIMAL, self.ACCEPTABLE)

    def __repr__(self):
        return "NLPResult(status={!r}, obj={:.6g}, iterations={})".format(self.status, self.obj, self.iterations)
