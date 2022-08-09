#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Some functions to help building pyomo optimisation problems
"""

import pyomo.environ as pyo

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['bounds',
           'initialisations',
           'objectives']


def initialisations(
    variable: str = 'f_tilde',
):
    """variable initialisations for pyomo

        Parameters
        ----------
        variable : str, optional
            * f_tilde: force, [fn+, fn-, fu, fv]

        Returns
        -------
        initialisations function for pyomo

    """
    def f_tilde_init(model, i):
        """initialise f ̃ with [1, 0, 1, 1]"""
        if i % 4 == 1:
            return 0.0
        return 1.0

    if variable == 'f_tilde':
        return f_tilde_init


def bounds(
    variable: str = 'd',
    d_bnd: float = 1e-3
):
    """variable bounds for pyomo

        Parameters
        ----------
        variable : str, optional
            * d: displacement
            * f: force, :math:`(fn, fu, fv)`
            * f_tilde: force, :math:`(fn^+, fn^-, fu, fv)`
        d_bnd : float, optional
            displacement bounds, -d_bnd <= d <= d_bnd

        Returns
        -------
        bounds function for pyomo

    """
    def f_tilde_bnds(model, i):
        """bounds of f ̃, f ̃ include [fn+, fn-, fu, fv]"""
        if i % 4 == 0 or i % 4 == 1:
            return pyo.NonNegativeReals
        return pyo.Reals

    def f_bnds(model, i):
        """bounds of f, f include [fn, fu, fv]"""
        if i % 3 == 0:
            return pyo.NonNegativeReals
        return pyo.Reals

    def d_bnds(model, i):
        return (-d_bnd, model.d[i], d_bnd)

    if variable == 'f':
        return f_bnds
    if variable == 'f_tilde':
        return f_tilde_bnds
    if variable == 'd':
        return d_bnds


def objectives(
    solver: str = 'cra',
    weights: tuple = (1e+0, 1e+0, 1e+6)
):
    """objective functions for pyomo

        Parameters
        ----------
        solver : str, optional
            * cra: CRA objective, :math:`W_{compression} * ||fn||_2^2 + W_{alpha} * ||alpha||_2^2`
            * cra_penalty: CRA penalty objective, :math:`W_{compression} * ||fn^+||_2^2 + W_{tension} * ||fn^-||_2^2 + W_{alpha} * ||alpha||_2^2`
            * rbe: RBE objective, :math:`W_{compression} * ||fn^+||_2^2 + W_{tension} * ||fn^-||_2^2`
        weights : tuple, optional
            weighting factors, :math:`(W_{alpha}, W_{compression}, W_{tension})`

        Returns
        -------
        objective function for pyomo

    """

    def obj_rbe(model):
        """RBE objective function"""
        return _obj_weights(model)

    def obj_cra(model):
        """CRA objective function"""
        alpha_sum = pyo.dot_product(model.alpha, model.alpha)
        f_sum = 0
        for i in model.fid:
            if i % 3 == 0:
                f_sum = f_sum + (model.f[i] * model.f[i])
        return f_sum + alpha_sum

    def obj_cra_penalty(model):
        """CRA penalty objective function"""
        alpha_sum = pyo.dot_product(model.alpha, model.alpha) * weights[0]  # alpha
        f_sum = _obj_weights(model)
        return alpha_sum + f_sum

    def _obj_weights(model):
        f_sum = 0
        for i in model.fid:
            if i % 4 == 1:
                f_sum = f_sum + (model.f[i] * model.f[i] * weights[2])  # tension
            elif i % 4 == 0:
                f_sum = f_sum + (model.f[i] * model.f[i] * weights[1])  # compression
        return f_sum

    if solver == 'cra':
        return obj_cra
    if solver == 'cra_penalty':
        return obj_cra_penalty
    if solver == 'rbe':
        return obj_rbe


def constraints(
    solver: str = 'cra',
    eps: float = 1e-4
):
    """constraint functions for pyomo"""

    def contact_con(model, i):
        """contact constraint"""
        dn = model.d[i * 3]
        fn = model.f[i * 3]
        return ((dn + eps) * fn, 0)

    def contact_con_penalty(model, i):
        """penalty contact constraint"""
        dn = model.d[i * 3]
        fn = model.f[i * 4]
        return ((dn + eps) * fn, 0)

    def fnp_con(model, i):
        """fn+ and fn- cannot coexist constraints"""
        return (model.f[i * 4] * model.f[i * 4 + 1], 0)

    def nonpen_con(m, t):
        """non penetration constraint"""
        return (0, m.d[t * 3] + eps, None)
