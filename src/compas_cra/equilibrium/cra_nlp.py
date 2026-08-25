"""The CRA optimisation problem as a plain sparse NLP, independent of pyomo.

This builds the same problem that :func:`~compas_cra.equilibrium.cra_solve`
formulates through pyomo — Eq.(11) of the CRA paper — directly as matrices and
callbacks with exact first and second derivatives:

.. math::

    \\begin{align}
        \\min_{f, \\delta q, \\alpha} \\quad & \\|f_n\\|_2^2 + \\|\\alpha\\|_2^2 \\\\
        \\textrm{s.t.} \\quad & A_{eq} f = -p \\\\
        & A_{fr} f \\le 0 \\\\
        & -\\eta \\le A_{eq}^\\intercal \\delta q \\le \\eta \\\\
        & f_n^i (\\delta d_n^i + \\varepsilon) = 0 \\\\
        & f_t^i + \\alpha^i \\delta d_t^i = 0 \\\\
        & \\delta d_n^i + \\varepsilon \\ge 0, \\quad f_n \\ge 0, \\; \\alpha \\ge 0
    \\end{align}

Variable layout: ``x = [f (3V), q (6F), alpha (V)]`` where ``V`` is the number of
contact vertices and ``F`` the number of free blocks. All nonlinearities are bilinear
(``f_n · d_n`` and ``alpha · d_t``), so the Jacobian has fixed sparsity and the exact
Lagrangian Hessian has constant per-constraint blocks.
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix

from ..nlp.problem import NLPProblem
from .cra_helper import equilibrium_setup
from .cra_helper import external_force_setup
from .cra_helper import friction_setup
from .cra_helper import num_free
from .cra_helper import num_vertices
from .cra_helper import unit_basis

__all__ = ["cra_problem", "cra_penalty_problem", "rbe_problem"]

INF = 1e19  # ipopt's conventional infinity


def _tilde_weights(v_count, w_comp, w_tens, w_fric):
    """Objective weights per f-tilde component, replicating the pyomo index rule.

    The historical rule is kept verbatim, including its quirk: the fv component
    (i % 4 == 3) is only weighted when i % 3 == 0, mirroring ``_obj_weights``.
    """
    w = np.zeros(4 * v_count)
    for i in range(4 * v_count):
        if i % 4 == 1:
            w[i] = w_tens
        elif i % 4 == 0:
            w[i] = w_comp
        elif i % 4 == 2 or i % 3 == 0:
            w[i] = w_fric
    return w


def _keep_nonempty_rows(matrix):
    """Rows of a csr matrix that have entries (pyomo skips structurally empty rows)."""
    keep = np.flatnonzero(np.diff(matrix.indptr))
    return csr_matrix(matrix[keep]), keep


def _tangent_map(basis, v_count, shift):
    """K with K[3i+x, shift*i + off_u] = u_i[x] and off_v column = v_i[x].

    Maps a force/displacement vector to its global tangential component per
    (vertex, xyz). ``shift`` is 3 for [fn, fu, fv] layouts and 4 for f-tilde.
    """
    off_u = shift - 2
    rows, cols, data = [], [], []
    for i in range(v_count):
        u = basis[3 * i + 1]
        v = basis[3 * i + 2]
        for x in range(3):
            if u[x]:
                rows.append(3 * i + x)
                cols.append(shift * i + off_u)
                data.append(u[x])
            if v[x]:
                rows.append(3 * i + x)
                cols.append(shift * i + off_u + 1)
                data.append(v[x])
    return csr_matrix(coo_matrix((data, (rows, cols)), shape=(3 * v_count, shift * v_count)))


def cra_problem(assembly, mu=0.84, density=1.0, d_bnd=1e-3, eps=1e-4):
    """Build the CRA NLP for an assembly.

    Returns
    -------
    tuple(:class:`~compas_cra.nlp.problem.NLPProblem`, dict)
        The problem, and a layout dict with the index ranges needed to interpret the
        solution vector (``nv``, ``nfree``, ``f``, ``q``, ``alpha`` slices).
    """
    V = num_vertices(assembly)
    F = num_free(assembly)
    if V == 0 or F == 0:
        raise ValueError("assembly has no contact interfaces or no free blocks")

    aeq = csr_matrix(equilibrium_setup(assembly))  # (6F, 3V)
    afr = csr_matrix(friction_setup(assembly, mu))  # (8V, 3V)
    p = external_force_setup(assembly, density).flatten()  # (6F,)
    basis = unit_basis(assembly)  # (3V, 3) rows [w, u, v] per vertex

    if aeq.shape != (6 * F, 3 * V):
        # csr_matrix construction infers the shape from the largest index; pad
        aeq.resize((6 * F, 3 * V))
    if afr.shape != (8 * V, 3 * V):
        afr.resize((8 * V, 3 * V))

    B = csr_matrix(aeq.T)  # (3V, 6F): d = B q
    Bn = B[0::3]  # (V, 6F): normal components d_n

    # K maps f -> tangential force vector per (vertex, xyz): K[3i+x, 3i+1] = u_i[x],
    # K[3i+x, 3i+2] = v_i[x]. The same matrix maps d -> tangential displacement, so the
    # friction alignment constraint is K f + alpha ∘ (K B q) = 0.
    k_rows, k_cols, k_data = [], [], []
    for i in range(V):
        u = basis[3 * i + 1]
        v = basis[3 * i + 2]
        for x in range(3):
            if u[x]:
                k_rows.append(3 * i + x)
                k_cols.append(3 * i + 1)
                k_data.append(u[x])
            if v[x]:
                k_rows.append(3 * i + x)
                k_cols.append(3 * i + 2)
                k_data.append(v[x])
    K = csr_matrix(coo_matrix((k_data, (k_rows, k_cols)), shape=(3 * V, 3 * V)))
    KB = csr_matrix(K @ B)  # (3V, 6F)

    # ------------------------------------------------------------------
    # variables: x = [f (3V), q (6F), alpha (V)]
    # ------------------------------------------------------------------
    nf, nq, na = 3 * V, 6 * F, V
    n = nf + nq + na
    of, oq, oa = 0, nf, nf + nq  # column offsets

    x_l = np.full(n, -INF)
    x_u = np.full(n, INF)
    x_l[of : of + nf : 3] = 0.0  # f_n >= 0
    x_l[oa:] = 0.0  # alpha >= 0

    x0 = np.zeros(n)
    x0[of : of + nf] = 1.0  # pyomo initialises every f component with 1

    # ------------------------------------------------------------------
    # constraints, in blocks. pyomo skips structurally empty rows, mirror that.
    # ------------------------------------------------------------------
    aeq_keep = np.flatnonzero(np.diff(aeq.indptr))
    aeq_k = csr_matrix(aeq[aeq_keep])
    m_eq = aeq_k.shape[0]
    m_fr = afr.shape[0]

    m = m_eq + m_fr + 3 * V + V + V + 3 * V
    o_fr = m_eq
    o_db = o_fr + m_fr
    o_ct = o_db + 3 * V
    o_np = o_ct + V
    o_ft = o_np + V

    g_l = np.empty(m)
    g_u = np.empty(m)
    g_l[:m_eq] = g_u[:m_eq] = -p[aeq_keep]
    g_l[o_fr:o_db] = -INF
    g_u[o_fr:o_db] = 0.0
    g_l[o_db:o_ct] = -d_bnd
    g_u[o_db:o_ct] = d_bnd
    g_l[o_ct:o_np] = g_u[o_ct:o_np] = 0.0  # f_n (d_n + eps) = 0
    g_l[o_np:o_ft] = -eps  # d_n >= -eps
    g_u[o_np:o_ft] = INF
    g_l[o_ft:] = g_u[o_ft:] = 0.0  # K f + alpha ∘ (K B q) = 0

    # ------------------------------------------------------------------
    # fixed Jacobian structure: constant blocks first, then the varying ones
    # ------------------------------------------------------------------
    aeq_coo = aeq_k.tocoo()
    afr_coo = afr.tocoo()
    b_coo = B.tocoo()
    bn_csr = csr_matrix(Bn)
    bn_coo = bn_csr.tocoo()
    k_coo = K.tocoo()
    kb_csr = KB
    kb_coo = kb_csr.tocoo()

    jr, jc = [], []
    const_parts = []

    # 1: equilibrium  A_eq f = -p
    jr.append(aeq_coo.row)
    jc.append(of + aeq_coo.col)
    const_parts.append(aeq_coo.data)
    # 2: friction  A_fr f <= 0
    jr.append(o_fr + afr_coo.row)
    jc.append(of + afr_coo.col)
    const_parts.append(afr_coo.data)
    # 3: displacement bounds  B q
    jr.append(o_db + b_coo.row)
    jc.append(oq + b_coo.col)
    const_parts.append(b_coo.data)
    # 5: no penetration  (B q)_n  (constant block, placed before the varying ones)
    jr.append(o_np + bn_coo.row)
    jc.append(oq + bn_coo.col)
    const_parts.append(bn_coo.data)
    # 6a: friction alignment, d(g)/d(f) = K
    jr.append(o_ft + k_coo.row)
    jc.append(of + k_coo.col)
    const_parts.append(k_coo.data)

    n_const = sum(part.shape[0] for part in const_parts)

    # varying entries
    # 4a: contact, d/d(f_n_i) = d_n_i + eps            (V entries, order i)
    jr.append(o_ct + np.arange(V))
    jc.append(of + 3 * np.arange(V))
    # 4b: contact, d/d(q_j) = f_n_i * Bn[i, j]         (nnz(Bn), csr order)
    jr.append(o_ct + bn_coo.row)
    jc.append(oq + bn_coo.col)
    # 6b: friction alignment, d/d(q_j) = alpha_i * KB[r, j]   (nnz(KB), csr order)
    jr.append(o_ft + kb_coo.row)
    jc.append(oq + kb_coo.col)
    # 6c: friction alignment, d/d(alpha_i) = (KB q)_r  (3V entries, order r)
    jr.append(o_ft + np.arange(3 * V))
    jc.append(oa + np.arange(3 * V) // 3)

    jac_rows = np.concatenate(jr)
    jac_cols = np.concatenate(jc)
    const_values = np.concatenate(const_parts)

    bn_entry_row = bn_coo.row  # vertex index of each Bn entry
    kb_entry_row = kb_coo.row  # (3i+x) row index of each KB entry
    fn_idx = of + 3 * np.arange(V)
    alpha_rep = np.arange(3 * V) // 3  # vertex index per ftdt row

    def split(x):
        return x[of : of + nf], x[oq : oq + nq], x[oa : oa + na]

    def objective(x):
        f, _, alpha = split(x)
        fn = f[0::3]
        return float(fn @ fn + alpha @ alpha)

    def gradient(x):
        f, _, alpha = split(x)
        grad = np.zeros(n)
        grad[fn_idx] = 2.0 * f[0::3]
        grad[oa:] = 2.0 * alpha
        return grad

    def constraints(x):
        f, q, alpha = split(x)
        d = B @ q
        dn = d[0::3]
        kf = K @ f
        kbq = KB @ q
        g = np.empty(m)
        g[:m_eq] = aeq_k @ f
        g[o_fr:o_db] = afr @ f
        g[o_db:o_ct] = d
        g[o_ct:o_np] = f[0::3] * (dn + eps)
        g[o_np:o_ft] = dn
        g[o_ft:] = kf + alpha[alpha_rep] * kbq
        return g

    def jacobian(x):
        f, q, alpha = split(x)
        dn = bn_csr @ q
        kbq = kb_csr @ q
        fn = f[0::3]
        return np.concatenate(
            [
                const_values,
                dn + eps,  # 4a
                fn[bn_entry_row] * bn_coo.data,  # 4b
                alpha[alpha_rep][kb_entry_row] * kb_coo.data,  # 6b
                kbq,  # 6c
            ]
        )

    # ------------------------------------------------------------------
    # exact Lagrangian Hessian, lower triangle (row >= col; duplicates are summed)
    # ------------------------------------------------------------------
    # The ftdt cross terms produce one entry per xyz component, all landing on the same
    # (alpha_i, q_j) cell. IPOPT would sum such duplicates, but ASL (the pyomo path)
    # merges them structurally — do the same so both solvers factorise the exact same
    # KKT matrix.
    kb_pair = (kb_entry_row // 3).astype(np.int64) * nq + kb_coo.col
    kb_unique, kb_inverse = np.unique(kb_pair, return_inverse=True)

    h_rows = np.concatenate(
        [
            fn_idx,  # objective: 2 sigma on f_n diagonal
            oa + np.arange(V),  # objective: 2 sigma on alpha diagonal
            oq + bn_coo.col,  # contact cross term (q_j, f_n_i)
            oa + (kb_unique // nq).astype(np.int32),  # ftdt cross term (alpha_i, q_j), merged
        ]
    )
    h_cols = np.concatenate(
        [
            fn_idx,
            oa + np.arange(V),
            fn_idx[bn_entry_row],
            oq + (kb_unique % nq).astype(np.int32),
        ]
    )

    def hessian(x, sigma, lam):
        lam_ct = lam[o_ct:o_np]
        lam_ft = lam[o_ft:]
        ftdt = np.bincount(kb_inverse, weights=lam_ft[kb_entry_row] * kb_coo.data, minlength=kb_unique.shape[0])
        return np.concatenate(
            [
                np.full(V, 2.0 * sigma),
                np.full(V, 2.0 * sigma),
                lam_ct[bn_entry_row] * bn_coo.data,
                ftdt,
            ]
        )

    problem = NLPProblem(
        n=n,
        x_l=x_l,
        x_u=x_u,
        g_l=g_l,
        g_u=g_u,
        x0=x0,
        objective=objective,
        gradient=gradient,
        constraints=constraints,
        jac_rows=jac_rows,
        jac_cols=jac_cols,
        jacobian=jacobian,
        hess_rows=h_rows,
        hess_cols=h_cols,
        hessian=hessian,
    )
    assert problem.jac_rows.shape[0] == n_const + V + bn_coo.nnz + kb_coo.nnz + 3 * V

    layout = {
        "nv": V,
        "nfree": F,
        "f": slice(of, of + nf),
        "q": slice(oq, oq + nq),
        "alpha": slice(oa, oa + na),
    }
    return problem, layout


def rbe_problem(assembly, mu=0.84, density=1.0):
    """Build the RBE quadratic program (paper Eq. 6) for an assembly.

    Variables are the penalty forces ``f-tilde = [fn+, fn-, fu, fv]`` per contact
    vertex; the constraints are linear (equilibrium and friction cone), the objective
    a diagonal quadratic with the historical component weights.
    """
    V = num_vertices(assembly)
    if V == 0:
        raise ValueError("assembly has no contact interfaces")

    aeq_b = csr_matrix(equilibrium_setup(assembly, penalty=True))
    afr_b = csr_matrix(friction_setup(assembly, mu, penalty=True))
    p = external_force_setup(assembly, density).flatten()

    n = 4 * V
    if aeq_b.shape[1] != n:
        aeq_b.resize((aeq_b.shape[0], n))
    if afr_b.shape != (8 * V, n):
        afr_b.resize((8 * V, n))

    weights = _tilde_weights(V, w_comp=1e0, w_tens=1e6, w_fric=1e0)

    x_l = np.full(n, -INF)
    x_u = np.full(n, INF)
    fn_idx = np.flatnonzero((np.arange(n) % 4) <= 1)
    x_l[fn_idx] = 0.0  # fn+ and fn- are nonnegative
    x0 = np.zeros(n)

    aeq_k, aeq_keep = _keep_nonempty_rows(aeq_b)
    m_eq = aeq_k.shape[0]
    m = m_eq + afr_b.shape[0]

    g_l = np.empty(m)
    g_u = np.empty(m)
    g_l[:m_eq] = g_u[:m_eq] = -p[aeq_keep]
    g_l[m_eq:] = -INF
    g_u[m_eq:] = 0.0

    aeq_coo = aeq_k.tocoo()
    afr_coo = afr_b.tocoo()
    jac_rows = np.concatenate([aeq_coo.row, m_eq + afr_coo.row])
    jac_cols = np.concatenate([aeq_coo.col, afr_coo.col])
    jac_values = np.concatenate([aeq_coo.data, afr_coo.data])

    def objective(x):
        return float(weights @ (x * x))

    def gradient(x):
        return 2.0 * weights * x

    def constraints(x):
        g = np.empty(m)
        g[:m_eq] = aeq_k @ x
        g[m_eq:] = afr_b @ x
        return g

    def jacobian(x):
        return jac_values

    w_nz = np.flatnonzero(weights)

    def hessian(x, sigma, lam):
        return 2.0 * sigma * weights[w_nz]

    problem = NLPProblem(
        n=n,
        x_l=x_l,
        x_u=x_u,
        g_l=g_l,
        g_u=g_u,
        x0=x0,
        objective=objective,
        gradient=gradient,
        constraints=constraints,
        jac_rows=jac_rows,
        jac_cols=jac_cols,
        jacobian=jacobian,
        hess_rows=w_nz,
        hess_cols=w_nz,
        hessian=hessian,
    )
    layout = {"nv": V, "nfree": 0, "f": slice(0, n), "q": None, "alpha": None}
    return problem, layout


def cra_penalty_problem(assembly, mu=0.84, density=1.0, d_bnd=1e-3, eps=1e-4):
    """Build the CRA penalty problem (paper Eq. 14) for an assembly.

    Variable layout: ``x = [f-tilde (4V), q (6F), alpha (V)]`` with
    ``f-tilde = [fn+, fn-, fu, fv]`` per contact vertex.
    """
    V = num_vertices(assembly)
    F = num_free(assembly)
    if V == 0 or F == 0:
        raise ValueError("assembly has no contact interfaces or no free blocks")

    aeq = csr_matrix(equilibrium_setup(assembly, penalty=False))  # (6F, 3V), for d
    aeq_b = csr_matrix(equilibrium_setup(assembly, penalty=True))  # (6F, 4V)
    afr_b = csr_matrix(friction_setup(assembly, mu, penalty=True))  # (8V, 4V)
    p = external_force_setup(assembly, density).flatten()
    d_basis = unit_basis(assembly)  # (3V, 3) rows [w, u, v]

    if aeq.shape != (6 * F, 3 * V):
        aeq.resize((6 * F, 3 * V))
    if aeq_b.shape != (6 * F, 4 * V):
        aeq_b.resize((6 * F, 4 * V))
    if afr_b.shape != (8 * V, 4 * V):
        afr_b.resize((8 * V, 4 * V))

    B = csr_matrix(aeq.T)  # (3V, 6F): d = B q
    Bn = csr_matrix(B[0::3])  # (V, 6F)
    Kf = _tangent_map(d_basis, V, shift=4)  # (3V, 4V): tangential force
    Kd = _tangent_map(d_basis, V, shift=3)  # (3V, 3V): tangential displacement
    KB = csr_matrix(Kd @ B)  # (3V, 6F)

    nf, nq, na = 4 * V, 6 * F, V
    n = nf + nq + na
    of, oq, oa = 0, nf, nf + nq

    weights = _tilde_weights(V, w_comp=1e0, w_tens=1e6, w_fric=1e0)

    x_l = np.full(n, -INF)
    x_u = np.full(n, INF)
    fn_idx = of + np.flatnonzero((np.arange(nf) % 4) <= 1)
    x_l[fn_idx] = 0.0  # fn+ and fn- are nonnegative
    x_l[oa:] = 0.0  # alpha >= 0
    x0 = np.zeros(n)

    aeq_k, aeq_keep = _keep_nonempty_rows(aeq_b)
    m_eq = aeq_k.shape[0]
    m_fr = afr_b.shape[0]

    m = m_eq + m_fr + 3 * V + V + V + V + 3 * V
    o_fr = m_eq
    o_db = o_fr + m_fr
    o_ct = o_db + 3 * V
    o_np = o_ct + V
    o_pm = o_np + V  # fn+ * fn- = 0
    o_ft = o_pm + V

    g_l = np.empty(m)
    g_u = np.empty(m)
    g_l[:m_eq] = g_u[:m_eq] = -p[aeq_keep]
    g_l[o_fr:o_db] = -INF
    g_u[o_fr:o_db] = 0.0
    g_l[o_db:o_ct] = -d_bnd
    g_u[o_db:o_ct] = d_bnd
    g_l[o_ct:o_np] = g_u[o_ct:o_np] = 0.0  # fn+ (d_n + eps) = 0
    g_l[o_np:o_pm] = -eps  # d_n >= -eps
    g_u[o_np:o_pm] = INF
    g_l[o_pm:o_ft] = g_u[o_pm:o_ft] = 0.0  # fn+ fn- = 0
    g_l[o_ft:] = g_u[o_ft:] = 0.0  # Kf f + alpha ∘ (Kd B q) = 0

    aeq_coo = aeq_k.tocoo()
    afr_coo = afr_b.tocoo()
    b_coo = B.tocoo()
    bn_coo = Bn.tocoo()
    kf_coo = Kf.tocoo()
    kb_coo = KB.tocoo()

    fnp_idx = of + 4 * np.arange(V)  # fn+ variable indices
    fnm_idx = fnp_idx + 1  # fn- variable indices
    alpha_rep = np.arange(3 * V) // 3

    jr = [
        aeq_coo.row,
        o_fr + afr_coo.row,
        o_db + b_coo.row,
        o_np + bn_coo.row,
        o_ft + kf_coo.row,
    ]
    jc = [
        of + aeq_coo.col,
        of + afr_coo.col,
        oq + b_coo.col,
        oq + bn_coo.col,
        of + kf_coo.col,
    ]
    const_values = np.concatenate([aeq_coo.data, afr_coo.data, b_coo.data, bn_coo.data, kf_coo.data])

    # varying entries
    jr += [
        o_ct + np.arange(V),  # contact d/d(fn+)
        o_ct + bn_coo.row,  # contact d/d(q)
        o_pm + np.arange(V),  # fn+ fn- d/d(fn+)
        o_pm + np.arange(V),  # fn+ fn- d/d(fn-)
        o_ft + kb_coo.row,  # ftdt d/d(q)
        o_ft + np.arange(3 * V),  # ftdt d/d(alpha)
    ]
    jc += [
        fnp_idx,
        oq + bn_coo.col,
        fnp_idx,
        fnm_idx,
        oq + kb_coo.col,
        oa + alpha_rep,
    ]
    jac_rows = np.concatenate(jr)
    jac_cols = np.concatenate(jc)

    def split(x):
        return x[of:oq], x[oq:oa], x[oa:]

    def objective(x):
        f, _, alpha = split(x)
        return float(weights @ (f * f) + alpha @ alpha)

    def gradient(x):
        f, _, alpha = split(x)
        grad = np.zeros(n)
        grad[of:oq] = 2.0 * weights * f
        grad[oa:] = 2.0 * alpha
        return grad

    def constraints(x):
        f, q, alpha = split(x)
        d = B @ q
        dn = d[0::3]
        fnp = f[0::4]
        fnm = f[1::4]
        kbq = KB @ q
        g = np.empty(m)
        g[:m_eq] = aeq_k @ f
        g[o_fr:o_db] = afr_b @ f
        g[o_db:o_ct] = d
        g[o_ct:o_np] = fnp * (dn + eps)
        g[o_np:o_pm] = dn
        g[o_pm:o_ft] = fnp * fnm
        g[o_ft:] = Kf @ f + alpha[alpha_rep] * kbq
        return g

    def jacobian(x):
        f, q, alpha = split(x)
        dn = Bn @ q
        kbq = KB @ q
        fnp = f[0::4]
        fnm = f[1::4]
        return np.concatenate(
            [
                const_values,
                dn + eps,
                fnp[bn_coo.row] * bn_coo.data,
                fnm,  # d(fn+ fn-)/d(fn+)
                fnp,  # d(fn+ fn-)/d(fn-)
                alpha[alpha_rep][kb_coo.row] * kb_coo.data,
                kbq,
            ]
        )

    # exact Lagrangian Hessian, lower triangle; ftdt duplicates merged as in cra_problem
    kb_pair = (kb_coo.row // 3).astype(np.int64) * nq + kb_coo.col
    kb_unique, kb_inverse = np.unique(kb_pair, return_inverse=True)
    w_nz = np.flatnonzero(weights)

    h_rows = np.concatenate(
        [
            of + w_nz,  # objective f diagonal
            oa + np.arange(V),  # objective alpha diagonal
            oq + bn_coo.col,  # contact cross (q_j, fn+_i)
            fnm_idx,  # fn+ fn- cross (fn-_i, fn+_i)
            oa + (kb_unique // nq).astype(np.int64),  # ftdt cross (alpha_i, q_j)
        ]
    )
    h_cols = np.concatenate(
        [
            of + w_nz,
            oa + np.arange(V),
            fnp_idx[bn_coo.row],
            fnp_idx,
            oq + (kb_unique % nq).astype(np.int64),
        ]
    )

    def hessian(x, sigma, lam):
        lam_ct = lam[o_ct:o_np]
        lam_pm = lam[o_pm:o_ft]
        lam_ft = lam[o_ft:]
        ftdt = np.bincount(kb_inverse, weights=lam_ft[kb_coo.row] * kb_coo.data, minlength=kb_unique.shape[0])
        return np.concatenate(
            [
                2.0 * sigma * weights[w_nz],
                np.full(V, 2.0 * sigma),
                lam_ct[bn_coo.row] * bn_coo.data,
                lam_pm,
                ftdt,
            ]
        )

    problem = NLPProblem(
        n=n,
        x_l=x_l,
        x_u=x_u,
        g_l=g_l,
        g_u=g_u,
        x0=x0,
        objective=objective,
        gradient=gradient,
        constraints=constraints,
        jac_rows=jac_rows,
        jac_cols=jac_cols,
        jacobian=jacobian,
        hess_rows=h_rows,
        hess_cols=h_cols,
        hessian=hessian,
    )
    layout = {
        "nv": V,
        "nfree": F,
        "f": slice(of, oq),
        "q": slice(oq, oa),
        "alpha": slice(oa, n),
    }
    return problem, layout
