"""Minimal smoke test for a built compas_cra wheel.

Solves min x^2 s.t. x >= 1 — enough to prove the extension loads, IPOPT runs, the
MUMPS linear solver works and callbacks round-trip. Needs only numpy, so it can run
in the bare test environments of the wheel-building CI, before the full test suite
(which needs compas and friends) gets a chance to run.
"""

import numpy as np

from compas_cra import _native as csn

print("compas_cra._native, IPOPT", csn.IPOPT_VERSION)

res = csn.solve_nlp(
    n=1,
    m=1,
    x_l=np.array([-10.0]),
    x_u=np.array([10.0]),
    g_l=np.array([1.0]),
    g_u=np.array([1e19]),
    x0=np.array([5.0]),
    jac_rows=np.array([0], dtype=np.int32),
    jac_cols=np.array([0], dtype=np.int32),
    eval_f=lambda x: float(x[0] ** 2),
    eval_grad=lambda x: np.array([2.0 * x[0]]),
    eval_g=lambda x: np.array([x[0]]),
    eval_jac=lambda x: np.array([1.0]),
    hess_rows=np.array([0], dtype=np.int32),
    hess_cols=np.array([0], dtype=np.int32),
    eval_h=lambda x, sigma, lam: np.array([2.0 * sigma]),
    options={"print_level": 0, "sb": "yes"},
)

assert res["status"] == 0, res["status_name"]
assert abs(res["x"][0] - 1.0) < 1e-6, res["x"]
assert abs(res["obj"] - 1.0) < 1e-6, res["obj"]
print("smoke test OK: x =", res["x"][0])
