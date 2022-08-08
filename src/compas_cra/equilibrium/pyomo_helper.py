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
           'f_tilde_init']


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


def obj_rbe(m):
    f_sum = 0
    for i in m.fid:
        if i % 4 == 1:
            f_sum = f_sum + (m.f[i] * m.f[i] * 1e+6)  # tension
        elif i % 4 == 0:
            f_sum = f_sum + (m.f[i] * m.f[i] * 0)  # compression
    return f_sum


def obj_cra(m):
    alpha_sum = pyo.dot_product(m.alpha, m.alpha)
    f_sum = 0
    for i in m.fid:
        if i % 3 == 0:
            f_sum = f_sum + (m.f[i] * m.f[i])
    return f_sum + alpha_sum


def obj_cra_penalty(m):
    alpha_sum = pyo.dot_product(m.alpha, m.alpha) * 1e+0  # alpha
    f_sum = 0
    for i in m.fid:
        if i % 4 == 1:
            f_sum = f_sum + (m.f[i] * m.f[i] * 1e+6)  # tension
        elif i % 4 == 0:
            f_sum = f_sum + (m.f[i] * m.f[i] * 0)  # compression
    return alpha_sum + f_sum
