#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Some pyomo functions.
"""

import pyomo.environ as pyo

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"

__all__ = ['f_bnds',
           'f_tilde_bnds',
           'f_tilde_init',
           'obj_rbe',
           'obj_cra',
           'obj_cra_penalty']


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


def f_tilde_init(model, i):
    """initialise f ̃ with [1, 0, 1, 1]"""
    if i % 4 == 1:
        return 0.0
    return 1.0


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
    alpha_sum = pyo.dot_product(model.alpha, model.alpha) * 1e+0  # alpha
    f_sum = _obj_weights(model)
    return alpha_sum + f_sum


def _obj_weights(model):
    f_sum = 0
    for i in model.fid:
        if i % 4 == 1:
            f_sum = f_sum + (model.f[i] * model.f[i] * 1e+6)  # tension
        elif i % 4 == 0:
            f_sum = f_sum + (model.f[i] * model.f[i] * 0)  # compression
    return f_sum
