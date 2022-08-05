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


def f_tilde_bnds(m, i):
    # bounds of f ̃, f ̃ include [fn+, fn-, fu, fv]
    if i % 4 == 0 or i % 4 == 1:
        return pyo.NonNegativeReals
    else:
        return pyo.Reals


def f_bnds(m, i):
    # bounds of f̃, f̃ include [fn, fu, fv]
    if i % 3 == 0:
        return pyo.NonNegativeReals
    else:
        return pyo.Reals


def f_tilde_init(m, i):
    # initialise f ̃ with [1, 0, 1, 1]
    if i % 4 == 1:
        return 0.0
    else:
        return 1.0

